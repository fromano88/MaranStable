function[b]=thermal_energy_Jacobian_V3d2(b, flowopt)
    ax=flowopt.ax;
    stability=flowopt.stability;
    
    rho=str2func(b.rho);
    lambda=str2func(b.lambda);
    cp=str2func(b.cp);
    
    drho=str2func(b.drho);
    dlambda=str2func(b.dlambda);
    dcp=str2func(b.dcp);
    ddcp=str2func(b.ddcp);
    
    % [J*r^ax*rho0*cp0*T0]_delta.xi*delta.eta
        dT_dxi.u = derivative_xi_eta(b,b.T,'xi','u');
        dT_dxi.u = [dT_dxi.u(:,1), dT_dxi.u, dT_dxi.u(:,end)];
        F.u=b.JA.u.*b.R.u.^ax.*rho(b.T.u).*b.T.u.*cp(b.T.u);
        T4_r = 0.5*dcp(b.T.u).*dT_dxi.u.*b.JA.u.*b.R.u.^ax.*rho(b.T.u).*b.T.u; % missing term 4

        u_s=zeros(b.K+2,b.J+2); u_n=u_s; w_w=u_s; w_e=u_s;

        u_s(1:end-2,2:end-1) = -F.u(1:end-1,2:end-1).*b.DETA.c;
        u_s(1:end-2,2:end-1) = u_s(1:end-2,2:end-1) + T4_r(1:end-1,2:end-1).*b.DETA.c.*b.DXI.c;
        u_n(2:end-1,2:end-1) = F.u(2:end,2:end-1).*b.DETA.c;
        u_n(2:end-1,2:end-1) = u_n(2:end-1,2:end-1) - T4_r(2:end,2:end-1).*b.DETA.c.*b.DXI.c;

        u_s=u_s.'; u_s=u_s(:);
        u_n=u_n.'; u_n=u_n(:);

        b.Ja.th_e.u=spdiags([u_s u_n],[-b.J-2 0],(b.K+2)*(b.J+2),(b.K+1)*(b.J+2));

    % [J*r^ax*rho0*cp0*dT0]_delta.eta*delta.xi
        dT_deta.w = derivative_xi_eta(b,b.T,'eta','w');
        dT_deta.w = [dT_deta.w(1,:); dT_deta.w; dT_deta.w(end,:)];
        F.w=b.JA.w.*b.R.w.^ax.*rho(b.T.w).*b.T.w.*cp(b.T.w);
        T4_z = 0.5*dcp(b.T.w).*dT_deta.w.*b.JA.w.*b.R.w.^ax.*rho(b.T.w).*b.T.w; % missing term 4

        w_w(2:end-1,2:end-1) = -F.w(2:end-1,1:end-1).*b.DXI.c;
        w_w(2:end-1,1:end-2) = w_w(2:end-1,1:end-2) + T4_z(2:end-1,1:end-1).*b.DETA.c.*b.DXI.c;
        w_e(2:end-1,3:end)   = F.w(2:end-1,2:end).*b.DXI.c;
        w_e(2:end-1,3:end)   = w_e(2:end-1,3:end) - T4_z(2:end-1,2:end).*b.DETA.c.*b.DXI.c;

        w_w=w_w.'; w_w=w_w(:);
        w_e=w_e.'; w_e=w_e(:);

        b.Ja.th_e.w=spdiags([w_w w_e],[0 1],(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
        b.Ja.th_e.w(:,1:b.J+2:end)=[];

    % [J*r^ax*u0*d(rho0*cp0*T0)/dT0]_delta.xi*delta.eta + [J*r^ax*w0*d(rho0*cp0*dT0)/dT0]_delta.eta*delta.xi
        F.u=b.JA.u.*b.R.u.^ax.*b.T.u.*(rho(b.T.u).*dcp(b.T.u)+drho(b.T.u).*cp(b.T.u)).*b.u.u;
        F.w=b.JA.w.*b.R.w.^ax.*b.T.w.*(rho(b.T.w).*dcp(b.T.w)+drho(b.T.w).*cp(b.T.w)).*b.w.w;
        T236_r = 0.5*b.JA.u.*b.R.u.^ax.*(-(drho(b.T.u).*b.T.u+rho(b.T.u)).*dcp(b.T.u).*dT_dxi.u - rho(b.T.u).*ddcp(b.T.u).*b.T.u.*dT_dxi.u).*b.u.u;     % missing terms T2, T3 and T6
        T236_z = 0.5*b.JA.w.*b.R.w.^ax.*(-(drho(b.T.w).*b.T.w+rho(b.T.w)).*dcp(b.T.w).*dT_deta.w - rho(b.T.w).*ddcp(b.T.w).*b.T.w.*dT_deta.w).*b.w.w;   % missing terms T2, T3 and T6

        % conductive terms
        % sign changed because of conductive terms!!!
            F1.u = -b.JA.u.*b.R.u.^ax.*(b.XI_r.u.^2+b.XI_z.u.^2).*dlambda(b.T.u);
            F2.u = -b.JA.u.*b.R.u.^ax.*b.ETA_z.u.*b.XI_z.u.*dlambda(b.T.u);
            F1.w = -b.JA.w.*b.R.w.^ax.*b.ETA_z.w.*b.XI_z.w.*dlambda(b.T.w);
            F2.w = -b.JA.w.*b.R.w.^ax.*b.ETA_z.w.^2.*dlambda(b.T.w);

            F_s.u = F1.u(1:end-1,2:end-1).*(b.T.T(2:end-1,2:end-1) - b.T.T(1:end-2,2:end-1))./b.DXI.rm(1:end-1,:) + F2.u(1:end-1,2:end-1).*(b.T.v(1:end-1,2:end) - b.T.v(1:end-1,1:end-1))./b.DETA.c;
            F_n.u = F1.u(2:end,2:end-1).*(b.T.T(3:end,2:end-1) - b.T.T(2:end-1,2:end-1))./b.DXI.rm(2:end,:) + F2.u(2:end,2:end-1).*(b.T.v(2:end,2:end) - b.T.v(2:end,1:end-1))./b.DETA.c;
            F_w.w = F1.w(2:end-1,1:end-1).*(b.T.v(2:end,1:end-1) - b.T.v(1:end-1,1:end-1))./b.DXI.c + F2.w(2:end-1,1:end-1).*(b.T.T(2:end-1,2:end-1) - b.T.T(2:end-1,1:end-2))./b.DETA.zm(:,1:end-1);
            F_e.w = F1.w(2:end-1,2:end).*(b.T.v(2:end,2:end) - b.T.v(1:end-1,2:end))./b.DXI.c + F2.w(2:end-1,2:end).*(b.T.T(2:end-1,3:end) - b.T.T(2:end-1,2:end-1))./b.DETA.zm(:,2:end);

        T_S=zeros(b.K+2,b.J+2); T_W=T_S; T_E=T_S; T_N=T_S; T_C=T_S;

        T_S(1:end-2,2:end-1) = -b.S.u(:,2:end-1).*(F.u(1:end-1,2:end-1) + F_s.u(1:end,:)).*b.DETA.c;
        T_S(1:end-2,2:end-1) = T_S(1:end-2,2:end-1) + T236_r(1:end-1,2:end-1).*b.DETA.c.*b.DXI.c;
        T_W(2:end-1,1:end-2) = -b.W.w(2:end-1,:).*(F.w(2:end-1,1:end-1) + F_w.w(:,1:end)).*b.DXI.c;
        T_W(2:end-1,1:end-2) = T_W(2:end-1,1:end-2) + T236_z(2:end-1,1:end-1).*b.DETA.c.*b.DXI.c;
        T_E(2:end-1,3:end)   = b.E.w(2:end-1,:).*(F.w(2:end-1,2:end) + F_e.w(:,1:end)).*b.DXI.c;
        T_E(2:end-1,3:end)   = T_E(2:end-1,3:end) - T236_z(2:end-1,2:end).*b.DETA.c.*b.DXI.c;
        T_N(3:end,2:end-1)   = b.N.u(:,2:end-1).*(F.u(2:end,2:end-1) + F_n.u(1:end,:)).*b.DETA.c;
        T_N(3:end,2:end-1)   = T_N(3:end,2:end-1) - T236_r(2:end,2:end-1).*b.DETA.c.*b.DXI.c;
        T_C(3:end-1,2:end-1) = T_C(3:end-1,2:end-1)-b.N.u(1:end-1,2:end-1).*(F.u(2:end-1,2:end-1) + F_s.u(2:end,:)).*b.DETA.c(2:end,:);
        T_C(2:end-1,3:end-1) = T_C(2:end-1,3:end-1)-b.E.w(2:end-1,1:end-1).*(F.w(2:end-1,2:end-1) + F_w.w(:,2:end)).*b.DXI.c(:,2:end);
        T_C(2:end-1,2:end-2) = T_C(2:end-1,2:end-2)+b.W.w(2:end-1,2:end).*(F.w(2:end-1,2:end-1) + F_e.w(:,1:end-1)).*b.DXI.c(:,1:end-1);
        T_C(2:end-2,2:end-1) = T_C(2:end-2,2:end-1)+b.S.u(2:end,2:end-1).*(F.u(2:end-1,2:end-1) + F_n.u(1:end-1,:)).*b.DETA.c(1:end-1,:);

        T_S=T_S.'; T_S=T_S(:);
        T_W=T_W.'; T_W=T_W(:);
        T_E=T_E.'; T_E=T_E(:);
        T_N=T_N.'; T_N=T_N(:);
        T_C=T_C.'; T_C=T_C(:);

        b.Ja.th_e.T=b.th_e.T+spdiags([T_S T_W T_C T_E T_N],[-b.J-2 -1 0 1 b.J+2],(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
        
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
                        XI_r_rlb.u=b.XI_r_rlb1.u;
                        XI_z_rlb_C.u=b.XI_z_rlb1_C.u;
                        XI_z_rlb_C.w=b.XI_z_rlb1_C.w;
                        XI_z_rlb_E.u=b.XI_z_rlb1_E.u;
                        XI_z_rlb_E.w=b.XI_z_rlb1_E.w;

                elseif o==2
                    % upper boundary
                        JA_rlb.u=b.JA_rlbend.u;
                        JA_rlb.w=b.JA_rlbend.w;
                        R_rlb.u=b.R_rlbend.u;
                        R_rlb.w=b.R_rlbend.w;
                        XI_r_rlb.u=b.XI_r_rlbend.u;
                        XI_z_rlb_C.u=b.XI_z_rlbend_C.u;
                        XI_z_rlb_C.w=b.XI_z_rlbend_C.w;
                        XI_z_rlb_E.u=b.XI_z_rlbend_E.u;
                        XI_z_rlb_E.w=b.XI_z_rlbend_E.w;

                end

                % convective terms
                    % [rho0*cp0*u0*T0*d(J*r^ax)/dRs]_delta.xi*delta.eta + [rho0*cp0*w0*dT0*d(J*r^ax)/dRs]_delta.eta*delta.xi
                        F.u=(JA_rlb.u.*b.R.u.^ax+ax*b.JA.u.*R_rlb.u).*rho(b.T.u).*cp(b.T.u).*b.T.u.*b.u.u;
                        FC.u=F.u(:,2:end-1);
                        F.w=(JA_rlb.w.*b.R.w.^ax+ax*b.JA.w.*R_rlb.w).*rho(b.T.w).*cp(b.T.w).*b.T.w.*b.w.w;
                        FC.w=F.w(2:end-1,:);

                % conductive terms with changed signs!
                    % -[(d(J*r^ax)/dRs*(xi_r^2+xi_z^2)+J*r^ax*2(xi_r*d(xi_r)/dRs+xi_z*d(xi_z)/dRs))*lambda0*dT0/xi]_delta.xi*delta.eta
                        F.u=-((JA_rlb.u.*b.R.u.^ax+ax*b.JA.u.*R_rlb.u).*(b.XI_r.u.^2+b.XI_z.u.^2)+b.JA.u.*b.R.u.^ax*2.*(b.XI_r.u.*XI_r_rlb.u+b.XI_z.u.*XI_z_rlb_C.u)).*lambda(b.T.u).*(b.T.T(2:end,:)-b.T.T(1:end-1,:));
                        FC.u=FC.u+F.u(:,2:end-1)./b.DXI.rm;
                        F.u=-b.JA.u.*b.R.u.^ax*2.*b.XI_z.u.*XI_z_rlb_E.u.*lambda(b.T.u).*(b.T.T(2:end,:)-b.T.T(1:end-1,:));
                        FE.u=F.u(:,2:end-1)./b.DXI.rm;

                    % -[(d(J*r^ax)/dRs*eta_z*xi_z+J*r^ax*eta_z*d(xi_z)/dRs))*lambda0*dT0/eta]_delta.xi*delta.eta
                        F.u=-((JA_rlb.u.*b.R.u.^ax+ax*b.JA.u.*R_rlb.u).*b.ETA_z.u.*b.XI_z.u+b.JA.u.*b.R.u.^ax.*b.ETA_z.u.*XI_z_rlb_C.u).*lambda(b.T.u);
                        FC.u=FC.u+F.u(:,2:end-1).*(b.T.v(:,2:end)-b.T.v(:,1:end-1))./b.DETA.rm;
                        F.u=-b.JA.u.*b.R.u.^ax.*b.ETA_z.u.*XI_z_rlb_E.u.*lambda(b.T.u);
                        FE.u=FE.u+F.u(:,2:end-1).*(b.T.v(:,2:end)-b.T.v(:,1:end-1))./b.DETA.rm;

                    % -[(d(J*r^ax)/dRs*eta_z*xi_z+J*r^ax*eta_z*d(xi_z)/dRs))*lambda0*dT0/xi]_delta.eta*delta.xi
                        F.w=-((JA_rlb.w.*b.R.w.^ax+ax*b.JA.w.*R_rlb.w).*b.ETA_z.w.*b.XI_z.w+b.JA.w.*b.R.w.^ax.*b.ETA_z.w.*XI_z_rlb_C.w).*lambda(b.T.w);
                        FC.w=FC.w+F.w(2:end-1,:).*(b.T.v(2:end,:)-b.T.v(1:end-1,:))./b.DXI.zm;
                        F.w=-b.JA.w.*b.R.w.^ax.*b.ETA_z.w.*XI_z_rlb_E.w.*lambda(b.T.w);
                        FE.w=F.w(2:end-1,:).*(b.T.v(2:end,:)-b.T.v(1:end-1,:))./b.DXI.zm;

                    % -[(d(J*r^ax)/dRs*eta_z.^2*lambda0*dT0/eta]_delta.eta*delta.xi
                        F.w=-(JA_rlb.w.*b.R.w.^ax+ax*b.JA.w.*R_rlb.w).*b.ETA_z.w.^2.*lambda(b.T.w).*(b.T.T(:,2:end)-b.T.T(:,1:end-1));
                        FC.w=FC.w+F.w(2:end-1,:)./b.DETA.zm;

                    % missing term
                        dT_dxi.w = derivative_xi_eta(b,b.T,'xi','w');
                        dT_dxi.w = [dT_dxi.w(1,:); dT_dxi.w; dT_dxi.w(end,:)];
                        dT_deta.u = derivative_xi_eta(b,b.T,'eta','u');
                        dT_deta.u = [dT_deta.u(:,1), dT_deta.u, dT_deta.u(:,end)];
                        F.u  = -rho(b.T.u).*b.T.u.*dcp(b.T.u).*(b.u.u.*dT_dxi.u + b.w.u.*dT_deta.u).*(JA_rlb.u.*b.R.u.^ax+ax*b.JA.u.*R_rlb.u);
                        FC.u = FC.u + F.u(:,2:end-1);
                        F.w  = -rho(b.T.w).*b.T.w.*dcp(b.T.w).*(b.u.w.*dT_dxi.w + b.w.w.*dT_deta.w).*(JA_rlb.w.*b.R.w.^ax+ax*b.JA.w.*R_rlb.w);
                        FC.w = FC.w + F.w(2:end-1,:);

                Rs_W=zeros(b.K+2,b.J); Rs_C=Rs_W; Rs_E=Rs_W; Rs_ad_W=Rs_W; Rs_ad_E=Rs_W;

                Rs_W(2:end-1,:)    = (-b.W.w(2:end-1,:).*FC.w(:,1:end-1)+FE.w(:,1:end-1)).*b.DXI.c-b.W.w(2:end-1,:).*(FE.u(2:end,:)-FE.u(1:end-1,:)).*b.DETA.c;
                Rs_E(2:end-1,:)    = (b.E.w(2:end-1,:).*FC.w(:,2:end)+FE.w(:,2:end)).*b.DXI.c+b.E.w(2:end-1,:).*(FE.u(2:end,:)-FE.u(1:end-1,:)).*b.DETA.c;
                Rs_C(2:end-1,:)    = (FC.u(2:end,:)-FC.u(1:end-1,:)).*b.DETA.c;
                Rs_C(2:end-1,:)    = Rs_C(2:end-1,:)-(FE.w(:,1:end-1)+FE.w(:,2:end)).*b.DXI.c;
                Rs_ad_W(2:end-1,:) = [zeros(b.K,1) -b.E.w(2:end-1,1:end-1).*(FC.w(:,2:end-1).*b.DXI.c(:,2:end)+(FE.u(2:end,2:end)-FE.u(1:end-1,2:end)).*b.DETA.c(:,2:end))];
                Rs_ad_E(2:end-1,:) = [b.W.w(2:end-1,2:end).*(FC.w(:,2:end-1).*b.DXI.c(:,1:end-1)+(FE.u(2:end,1:end-1)-FE.u(1:end-1,1:end-1)).*b.DETA.c(:,1:end-1)) zeros(b.K,1)];
                Rs_C               = Rs_C + Rs_ad_W;
                Rs_C               = Rs_C + Rs_ad_E;

                %Rs_W=Rs_W(:);
                %Rs_C=Rs_C(:);
                %Rs_E=Rs_E(:);

                index=kron(ones(1,b.J),2:b.J+2:(b.K+2)*(b.J+2))+ kron(0:(b.K+2)*(b.J+2)+1:(b.K+2)*(b.J+2)*(b.J),ones(1,(b.K+2)));
                Rs1=sparse((b.K+2)*(b.J+2),b.J); Rs2=Rs1; Rs3=Rs1;
                Rs1(index)=Rs_W;
                Rs2(index)=Rs_C;
                Rs3(index)=Rs_E;

                Ja.th_e.Rs_all=[Rs1 sparse((b.K+2)*(b.J+2),2)]+[sparse((b.K+2)*(b.J+2),1) Rs2 sparse((b.K+2)*(b.J+2),1)]+[sparse((b.K+2)*(b.J+2),2) Rs3];
                Ja.th_e.Rs=[];
                Ja_counter.th_e.Rs=[];

                for n=1:length(b.bc.z)/2
                    if (o==1 && strcmp(b.bc.z{n}(1),'s')==1 && strcmp(b.bc.z{n}(2),'s')==1) || (o==2 && strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'s')==1 && strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'s')==1)
                        Ja.th_e.Rs_part=Ja.th_e.Rs_all(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                        if (n>1 && strcmp(b.bc.z{(o-1)*index_bc_z+n-1}(1),'s')~=1 && strcmp(b.bc.z{(o-1)*index_bc_z+n-1}(2),'s')~=1)
                            Ja.th_e.Rs_part(find(b.z.w==b.geom.z(n))+1:b.J+2:end,1)=Ja.th_e.Rs_part(find(b.z.w==b.geom.z(n))+1:b.J+2:end,1)-Rs_ad_W(:,b.z.w==b.geom.z(n));
                            Ja.th_e.Rs_part(find(b.z.w==b.geom.z(n))+b.J+2:b.J+2:end-b.J-2,1)=FE.w(:,b.z.w==b.geom.z(n)).*b.DXI.c(:,1);

                        end
                        if (n<length(b.bc.z)/2 && strcmp(b.bc.z{(o-1)*index_bc_z+n+1}(1),'s')~=1 && strcmp(b.bc.z{(o-1)*index_bc_z+n+1}(2),'s')~=1)
                            Ja.th_e.Rs_part(find(b.z.w==b.geom.z(n+1)):b.J+2:end,end)=Ja.th_e.Rs_part(find(b.z.w==b.geom.z(n+1)):b.J+2:end,end)-Rs_ad_E(:,find(b.z.w==b.geom.z(n+1))-1);
                            Ja.th_e.Rs_part(find(b.z.w==b.geom.z(n+1))+b.J+3:b.J+2:end-b.J-2,end)=FE.w(:,b.z.w==b.geom.z(n+1)).*b.DXI.c(:,1);

                        end

                        if eval(['b.' b.bc.z{(o-1)*index_bc_z+n} '.o==1'])
                            Ja.th_e.Rs=[Ja.th_e.Rs Ja.th_e.Rs_part];

                        elseif eval(['b.' b.bc.z{(o-1)*index_bc_z+n} '.o==0'])
                            Ja_counter.th_e.Rs=[Ja_counter.th_e.Rs Ja.th_e.Rs_part];

                        end

                    end

                end

                if o==1
                    b.Ja.th_e.Rs1=Ja.th_e.Rs;
                    b.Ja_counter.th_e.Rs1=Ja_counter.th_e.Rs;

                else
                    b.Ja.th_e.Rsend=Ja.th_e.Rs;
                    b.Ja_counter.th_e.Rsend=Ja_counter.th_e.Rs;

                end

            end
            
        end
        
    end
