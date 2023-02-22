function[b]=continuity(b, flowopt)
    ax=flowopt.ax;
    stability=flowopt.stability;
    rho=str2func(b.rho);
    
    % [J*r^ax*rho0*u0]_delta.xi*delta.eta
        F.u=b.JA.u.*b.R.u.^ax.*rho(b.T.u);

        u_s = [sparse(b.K,1) -F.u(1:end-1,2:end-1).*b.DETA.c sparse(b.K,1)];
        u_n = [sparse(b.K,1) F.u(2:end,2:end-1).*b.DETA.c sparse(b.K,1)];

        u_s=u_s.'; u_s=u_s(:);
        u_n=u_n.'; u_n=u_n(:);

        b.C.u=spdiags([u_s u_n],[0 b.J+2],(b.K)*(b.J+2),(b.K+1)*(b.J+2));
        b.C.u(1:b.J+2:end,:)=[];
        b.C.u(b.J+1:b.J+1:end,:)=[];

    % [J*r^ax*rho0*w0]_delta.eta*delta.xi
        F.w=b.JA.w.*b.R.w.^ax.*rho(b.T.w);

        w_w=zeros(b.K+2,b.J+2); w_e=w_w;

        w_w(2:end-1,2:end-1) = -F.w(2:end-1,1:end-1).*b.DXI.c;
        w_e(2:end-1,2:end-1) = F.w(2:end-1,2:end).*b.DXI.c;

        w_w=w_w.'; w_w=w_w(:);
        w_e=w_e.'; w_e=w_e(:);

        b.C.w=spdiags([w_w w_e],[0 1],(b.K+1)*(b.J+2),(b.K+2)*(b.J+2));
        b.C.w(1:b.J+1,:)=[];
        b.C.w(1:b.J+2:end,:)=[];
        b.C.w(1:b.J+1:end,:)=[];
        b.C.w(:,1:b.J+2:end)=[];
        
    b.C.p=sparse(b.K*b.J, b.K*b.J);
        
    if ((max(strncmp('o', b.bc.z(:),1)) || max(strncmp('o', b.bc.r(:),1)))~=1) && (stability~=1 || (stability==1 && flowopt.m==0))
         ri=find(abs(b.r.w(2:end-1)-b.rp)==min(min(abs(b.r.w(2:end-1)-b.rp))));
         zi=find(abs(b.z.u(2:end-1)-b.zp)==min(min(abs(b.z.u(2:end-1)-b.zp))));
         b.C.u((ri(1)-1)*b.J+zi(1),:)=0;
         b.C.w((ri(1)-1)*b.J+zi(1),:)=0;
         b.C.p((ri(1)-1)*b.J+zi(1),(ri(1)-1)*b.J+zi(1))=1; 
    end
    
    index_bc_z=length(b.bc.z)/2;
    
    if stability~=1
        for o=1:2
            if max(strncmp('ss', {b.bc.z{(o-1)*index_bc_z+1:o*index_bc_z}},2))
                if o==1
                    % lower boundary
                        C.R=b.C.R1;
                        counter.C.R=b.counter.C.R1;

                elseif o==2
                    % upper boundary
                        C.R=b.C.Rend;
                        counter.C.R=b.counter.C.Rend;

                end
                C.Rs=[];
                counter.C.Rs=[];
                for n=1:length(b.bc.z)/2
                    if (o==1 && strcmp(b.bc.z{n}(1),'s')==1 && strcmp(b.bc.z{n}(2),'s')==1) || (o==2 && strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'s')==1 && strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'s')==1)
                        if eval(['b.' b.bc.z{(o-1)*index_bc_z+n} '.o==1'])
                            C.Rs=[C.Rs sparse(b.K*b.J,length(b.z.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1))))];

                        elseif eval(['b.' b.bc.z{(o-1)*index_bc_z+n} '.o==0'])
                            counter.C.Rs=[counter.C.Rs sparse(b.K*b.J,length(b.z.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1))))];

                        end

                    end

                end
                %}
                if (max(strncmp('o', b.bc.z(:),1)) || max(strncmp('o', b.bc.r(:),1)))~=1
                    if max(strncmp('ss', b.bc.z(:),2))==1 && isfield(b, 'V_l')==1
                        if isempty(C.R)==0
                            C.Rs((ri(1)-1)*b.J+zi(1),:)=C.R;

                        end
                        if isempty(counter.C.R)==0
                            counter.C.Rs((ri(1)-1)*b.J+zi(1),:)=counter.C.R;

                        end
                        b.C.p((ri(1)-1)*b.J+zi(1),(ri(1)-1)*b.J+zi(1))=0;

                    end
                end

                if o==1
                    b.C.Rs1=C.Rs;
                    b.counter.C.Rs1=counter.C.Rs;

                else
                    b.C.Rsend=C.Rs;
                    b.counter.C.Rsend=counter.C.Rs;

                end

            end

        end
        
    end