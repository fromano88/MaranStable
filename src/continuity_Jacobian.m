function[b]=continuity_Jacobian(b, flowopt)
    ax=flowopt.ax;
    energy=flowopt.energy;
    stability=flowopt.stability;
    
    rho=str2func(b.rho);
    
    drho=str2func(b.drho);
    
    if (max(strncmp('o', b.bc.z(:),1)) || max(strncmp('o', b.bc.r(:),1)))~=1
        ri=find(abs(b.r.w(2:end-1)-b.rp)==min(min(abs(b.r.w(2:end-1)-b.rp))));
        zi=find(abs(b.z.u(2:end-1)-b.zp)==min(min(abs(b.z.u(2:end-1)-b.zp))));
            
    end
        
    if energy==1
        % [J*r^ax*u0*drho0/dT0]_delta.xi*delta.eta + [J*r^ax*w0*drho0/dT0]_delta.eta*delta.xi
            F.u=b.JA.u.*b.R.u.^ax.*drho(b.T.u).*b.u.u;
            F.w=b.JA.w.*b.R.w.^ax.*drho(b.T.w).*b.w.w;

            T_S=zeros(b.K+2,b.J+2); T_W=T_S; T_E=T_S; T_N=T_S; T_C=T_S;

            T_S(1:end-2,2:end-1) = -b.S.u(:,2:end-1).*F.u(1:end-1,2:end-1).*b.DETA.c;
            T_W(2:end-1,1:end-2) = -b.W.w(2:end-1,:).*F.w(2:end-1,1:end-1).*b.DXI.c;
            T_E(2:end-1,3:end)   = b.E.w(2:end-1,:).*F.w(2:end-1,2:end).*b.DXI.c;
            T_N(3:end,2:end-1)   = b.N.u(:,2:end-1).*F.u(2:end,2:end-1).*b.DETA.c;
            T_C(3:end-1,2:end-1) = T_C(3:end-1,2:end-1)-b.N.u(1:end-1,2:end-1).*F.u(2:end-1,2:end-1).*b.DETA.c(2:end,:);
            T_C(2:end-1,3:end-1) = T_C(2:end-1,3:end-1)-b.E.w(2:end-1,1:end-1).*F.w(2:end-1,2:end-1).*b.DXI.c(:,2:end);
            T_C(2:end-1,2:end-2) = T_C(2:end-1,2:end-2)+b.W.w(2:end-1,2:end).*F.w(2:end-1,2:end-1).*b.DXI.c(:,1:end-1);
            T_C(2:end-2,2:end-1) = T_C(2:end-2,2:end-1)+b.S.u(2:end,2:end-1).*F.u(2:end-1,2:end-1).*b.DETA.c(1:end-1,:);

            T_S=T_S.'; T_S=T_S(:);
            T_W=T_W.'; T_W=T_W(:);
            T_C=T_C.'; T_C=T_C(:);
            T_E=T_E.'; T_E=T_E(:);
            T_N=T_N.'; T_N=T_N(:);

            b.Ja.C.T=spdiags([T_S T_W T_C T_E T_N],[-b.J-2 -1 0 1 b.J+2],(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
            b.Ja.C.T(1:b.J+2,:)=[];
            b.Ja.C.T(1:b.J+2:end,:)=[];
            b.Ja.C.T(b.J+1:b.J+1:end,:)=[];
            b.Ja.C.T((b.K)*(b.J)+1:end,:)=[];
            
            if (max(strncmp('o', b.bc.z(:),1)) || max(strncmp('o', b.bc.r(:),1)))~=1 && (stability~=1 || (stability==1 && flowopt.m==0))
                b.Ja.C.T((ri(1)-1)*b.J+zi(1),:)=0;
                
            end

    end
        
    if stability~=1
        index_bc_z=length(b.bc.z)/2;
    
        for o=1:2
            if max(strncmp('ss', {b.bc.z{(o-1)*index_bc_z+1:o*index_bc_z}},2))
                if o==1
                    % lower boundary
                        JA_rlb.u=b.JA_rlb1.u;
                        JA_rlb.w=b.JA_rlb1.w;
                        R_rlb.u=b.R_rlb1.u;
                        R_rlb.w=b.R_rlb1.w;
                        Ja.C.R=b.Ja.C.R1;
                        Ja.counter.C.R=b.Ja.counter.C.R1;

                elseif o==2
                    % upper boundary
                        JA_rlb.u=b.JA_rlbend.u;
                        JA_rlb.w=b.JA_rlbend.w;
                        R_rlb.u=b.R_rlbend.u;
                        R_rlb.w=b.R_rlbend.w;
                        Ja.C.R=b.Ja.C.Rend;
                        Ja.counter.C.R=b.Ja.counter.C.Rend;

                end

                F.u = (JA_rlb.u.*b.R.u.^ax+ax*b.JA.u.*R_rlb.u).*rho(b.T.u).*b.u.u;
                F.w = (JA_rlb.w.*b.R.w.^ax+ax*b.JA.w.*R_rlb.w).*rho(b.T.w).*b.w.w;

                Rs_W    = -b.W.w(2:end-1,:).*F.w(2:end-1,1:end-1).*b.DXI.c;
                Rs_E    = b.E.w(2:end-1,:).*F.w(2:end-1,2:end).*b.DXI.c;
                Rs_C    = (F.u(2:end,2:end-1)-F.u(1:end-1,2:end-1)).*b.DETA.c;
                Rs_ad_W = [zeros(b.K,1) -b.E.w(2:end-1,1:end-1).*F.w(2:end-1,2:end-1).*b.DXI.c(:,2:end)];
                Rs_ad_E = [b.W.w(2:end-1,2:end).*F.w(2:end-1,2:end-1).*b.DXI.c(:,1:end-1) zeros(b.K,1)];
                Rs_C    = Rs_C + Rs_ad_W;
                Rs_C    = Rs_C + Rs_ad_E;

                %Rs_W=Rs_W(:);
                %Rs_C=Rs_C(:);
                %Rs_E=Rs_E(:);

                index=kron(ones(1,b.J),1:b.J:(b.K)*(b.J))+ kron(0:(b.K)*(b.J)+1:(b.K)*(b.J)*(b.J),ones(1,(b.K)));
                Rs1=sparse((b.K)*(b.J),b.J); Rs2=Rs1; Rs3=Rs1;
                Rs1(index)=Rs_W;
                Rs2(index)=Rs_C;
                Rs3(index)=Rs_E;

                Ja.C.Rs_all=[Rs1 sparse(b.K*b.J,2)]+[sparse(b.K*b.J,1) Rs2 sparse(b.K*b.J,1)]+[sparse(b.K*b.J,2) Rs3];
                Ja.C.Rs=[];
                Ja_counter.C.Rs=[];

                for n=1:length(b.bc.z)/2
                    if (o==1 && strcmp(b.bc.z{n}(1),'s')==1 && strcmp(b.bc.z{n}(2),'s')==1) || (o==2 && strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'s')==1 && strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'s')==1)
                        Ja.C.Rs_part=Ja.C.Rs_all(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                        if n>1 && strcmp(b.bc.z{(o-1)*index_bc_z+n-1}(1),'s')~=1 && strcmp(b.bc.z{(o-1)*index_bc_z+n-1}(2),'s')~=1
                            Ja.C.Rs_part(length(b.z.w(b.z.w<b.geom.z(n)))+1:b.J:end,1)=Ja.C.Rs_part(length(b.z.w(b.z.w<b.geom.z(n)))+1:b.J:end,1)-Rs_ad_W(:,length(b.z.w(b.z.w<b.geom.z(n)))+1);
                            Ja.C.Rs_part(length(b.z.w(b.z.w<b.geom.z(n))):b.J:end,1)=0;

                        end
                        if n<length(b.bc.z)/2 && strcmp(b.bc.z{(o-1)*index_bc_z+n+1}(1),'s')~=1 && strcmp(b.bc.z{(o-1)*index_bc_z+n+1}(2),'s')~=1
                            Ja.C.Rs_part(length(b.z.w(b.z.w<b.geom.z(n+1))):b.J:end,end)=Ja.C.Rs_part(length(b.z.w(b.z.w<b.geom.z(n+1))):b.J:end,end)-Rs_ad_E(:,length(b.z.w(b.z.w<b.geom.z(n+1))));
                            Ja.C.Rs_part(length(b.z.w(b.z.w<b.geom.z(n+1)))+1:b.J:end,end)=0;

                        end

                        if eval(['b.' b.bc.z{(o-1)*index_bc_z+n} '.o==1'])
                            Ja.C.Rs=[Ja.C.Rs Ja.C.Rs_part];

                        elseif eval(['b.' b.bc.z{(o-1)*index_bc_z+n} '.o==0'])
                            Ja_counter.C.Rs=[Ja_counter.C.Rs Ja.C.Rs_part];

                        end

                    end

                end
                if (max(strncmp('o', b.bc.z(:),1)) || max(strncmp('o', b.bc.r(:),1)))~=1
                    if max(strncmp('ss', b.bc.z(:),2))==1 && isfield(b, 'V_l')==1
                        if isempty(Ja.C.Rs)==0
                            Ja.C.Rs((ri(1)-1)*b.J+zi(1),:)=Ja.C.R;

                        end
                        if isempty(Ja_counter.C.Rs)==0
                            Ja_counter.C.Rs((ri(1)-1)*b.J+zi(1),:)=Ja.counter.C.R;

                        end

                    elseif isfield(b,'p_fluid')==1
                        if isempty(Ja.C.Rs)==0
                            Ja.C.Rs((ri(1)-1)*b.J+zi(1),:)=0;

                        end
                        if isempty(Ja_counter.C.Rs)==0
                            Ja_counter.C.Rs((ri(1)-1)*b.J+zi(1),:)=0;

                        end

                    end
                end

                if o==1
                    b.Ja.C.Rs1=Ja.C.Rs;
                    b.Ja_counter.C.Rs1=Ja_counter.C.Rs;

                else
                    b.Ja.C.Rsend=Ja.C.Rs;
                    b.Ja_counter.C.Rsend=Ja_counter.C.Rs;

                end
                    
            end
            
        end
        
    end