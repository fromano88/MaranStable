if o==1
    % lower boundary
        index_bc_z=0;
        index_bc_r=1;
        K=1;
        R_z.u=b.R1_z.u;
        R_z.w=b.R1_z.w;
        R_zz.u=b.R1_zz.u;

elseif o==2
    % upper boundary
        index_bc_z=length(b.bc.z)/2;
        index_bc_r=length(b.bc.r)/2;
        K=b.K+1;
        R_z.u=b.Rend_z.u;
        R_z.w=b.Rend_z.w;
        R_zz.u=b.Rend_zz.u;

end

N.u=sqrt(1+R_z.u(1,:).^2);

% equal velocity
    equal_u=1./N.u./b.XI_r.u(K,:);    
    
if stability~=1 && max(strncmp('ss', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},2))
    % normal stresses
        nabla_dot_n=-R_zz.u(1,:)./N.u.^3+ax./(b.R.u(K,:).^ax.*N.u);

        F.u=mu(b.T.u(K,2:end-1)).*2.*b.XI_r.u(K,2:end-1).*b.dz.c.*N.u(1,2:end-1);

        u_Wsw=zeros(1,b.J+2); u_s=u_Wsw; u_Ese=u_Wsw; u_n=u_Wsw;

        u_Wsw(1,2:end-1) = 2*mu(b.T.u(K,2:end-1))./N.u(1,2:end-1).*R_z.u(1,2:end-1).*b.ETA_z.u(K,2:end-1)./b.XI_r.w(K+(o-1),1:end-1)./b.DETA.c(1,:)./2.*b.dz.c;
        u_Wsw(1,2)       = 2*u_Wsw(1,2);
        u_Ese(1,2:end-1) = -2*mu(b.T.u(K,2:end-1))./N.u(1,2:end-1).*R_z.u(1,2:end-1).*b.ETA_z.u(K,2:end-1)./b.XI_r.w(K+(o-1),2:end)./b.DETA.c(1,:)./2.*b.dz.c;
        u_Ese(1,end-1)   = 2*u_Ese(1,end-1);
        u_s(1,3:end-1)   = u_s(1,3:end-1) + u_Wsw(1,3:end-1);
        u_s(1,2:end-2)   = u_s(1,2:end-2) + u_Ese(1,2:end-2);
        u_s(1,2:end-1)   = u_s(1,2:end-1) + (-1)^o*(F.u./b.XI_r.u(K,2:end-1)./b.DXI.c(K-(o-1),:) + 2/3*mu(b.T.u(K,2:end-1))./b.DXI.c(K-(o-1),:).*b.dz.c.*N.u(1,2:end-1));
        u_n(1,2:end-1)   = -(-1)^o*(F.u./b.XI_r.u(K-(-1)^o,2:end-1)./b.DXI.c(K-(o-1),:) - 2/3./b.JA.u(K,2:end-1).*mu(b.T.u(K,2:end-1))./b.R.u(K,2:end-1).^ax.*b.JA.u(K-(-1)^o,2:end-1).*b.R.u(K-(-1)^o,2:end-1).^ax./b.DXI.c(K-(o-1),:).*b.dz.c.*N.u(1,2:end-1));

        F.u=-mu(b.T.u(K,2:end-1)).*2.*b.XI_r.u(K,2:end-1).*b.dz.c.*N.u(1,2:end-1);
        F2.u=-mu(b.T.u(K,2:end-1))./N.u(1,2:end-1).*2.*R_z.u(1,2:end-1).*(R_z.w(1,2:end)-R_z.w(1,1:end-1))./b.DETA.c(1,:).*b.dz.c;
        F3.u=-2/3./b.JA.u(K,2:end-1).*mu(b.T.u(K,2:end-1))./b.R.u(K,2:end-1).^ax.*b.dz.c.*N.u(1,2:end-1);

        w_sw=zeros(1,b.J+2); w_se=w_sw; w_nw=w_sw; w_ne=w_sw;

        w_sw(1,2:end-1) = (-1)^o*F.u.*(b.JA.T(K+(o-1),2:end-1).*b.XI_z.T(K+(o-1),2:end-1)+R_z.u(1,2:end-1)./b.ETA_z.T(K+(o-1),2:end-1))./b.DXI.rm(K,:)./2 + 0.5*F2.u - F3.u.*b.JA.w(K+(o-1),1:end-1).*b.R.w(K+(o-1),1:end-1).^ax./b.DETA.c(1,:);
        w_se(1,2:end-1) = (-1)^o*F.u.*(b.JA.T(K+(o-1),2:end-1).*b.XI_z.T(K+(o-1),2:end-1)+R_z.u(1,2:end-1)./b.ETA_z.T(K+(o-1),2:end-1))./b.DXI.rm(K,:)./2 + 0.5*F2.u + F3.u.*b.JA.w(K+(o-1),2:end).*b.R.w(K+(o-1),2:end).^ax./b.DETA.c(1,:);
        w_nw(1,2:end-1) = -(-1)^o*F.u.*(b.JA.T(K-(o-2),2:end-1).*b.XI_z.T(K-(o-2),2:end-1)+R_z.u(1,2:end-1)./b.ETA_z.T(K-(o-2),2:end-1))./b.DXI.rm(K,:)./2;
        w_ne(1,2:end-1) = -(-1)^o*F.u.*(b.JA.T(K-(o-2),2:end-1).*b.XI_z.T(K-(o-2),2:end-1)+R_z.u(1,2:end-1)./b.ETA_z.T(K-(o-2),2:end-1))./b.DXI.rm(K,:)./2;

        p_C=zeros(1,b.J+2); p_N=p_C;

        %{
        p_C(1,2:end-1) = -(b.R_cyl.T(K+(o-1),2:end-1)-b.R_cyl.p(K-(-1)^o*o,:))./(b.R_cyl.p(K-(o-1),:)-b.R_cyl.p(K-(-1)^o*o,:));
        p_N(1,2:end-1) = -(b.R_cyl.p(K-(o-1),:)-b.R_cyl.T(K+(o-1),2:end-1))./(b.R_cyl.p(K-(o-1),:)-b.R_cyl.p(K-(-1)^o*o,:));
        %}
        %%{
        p_C(1,2:end-1) = -1.5.*b.dz.c.*N.u(1,2:end-1);
        p_N(1,2:end-1) = 0.5.*b.dz.c.*N.u(1,2:end-1);
        %}
        
    % Laplace pressure
        if isfield(b, 'sigma')==1
            LP=(-1)^o*sigma(b.T.u(K,:)).*nabla_dot_n.*[b.dz.zm(1) b.dz.c b.dz.zm(end)].*N.u;

        end

    % hydrostatic pressure
    % hp was implemented twice: once inside the field quantity p and once
    % again here --> therefore HP is set to 0 here --> don't change!!!
        HP=-0*(-1)^o*g*rho(b.T.u(K,:)).*b.z.u(1,:).*[b.dz.zm(1) b.dz.c b.dz.zm(end)].*N.u;

    if energy==1
        % derivative of normal stresses with respect to T
            % dmu/dT*2*xi_r*d(u0/xi_r)/delta.xi
                T_C=dmu(b.T.u(K,:)).*2.*b.XI_r.u(K,:).*(b.u.u(K+2-o,:)./b.XI_r.u(K+2-o,:)-b.u.u(K-(o-1),:)./b.XI_r.u(K-(o-1),:))./[b.DXI.c(K-(o-1),1) b.DXI.c(K-(o-1),:) b.DXI.c(K-(o-1),end)];

            % -2/3*dmu/dT/(J*Rs)*d(J*R*u0))/delta.xi
                T_C=T_C-2/3*dmu(b.T.u(K,:))./(b.JA.u(K,:).*b.R.u(K,:).^ax).*(b.JA.u(K+2-o,:).*b.R.u(K+2-o,:).^ax.*b.u.u(K+2-o,:)-b.JA.u(K-(o-1),:).*b.R.u(K-(o-1),:).^ax.*b.u.u(K-(o-1),:))./[b.DXI.c(K-(o-1),1) b.DXI.c(K-(o-1),:) b.DXI.c(K-(o-1),end)];

            % -dmu/dT*2*xi_r*d(J*xi_z*w0)/delta.xi
                T_C=T_C-dmu(b.T.u(K,:)).*2.*b.XI_r.u(K,:).*(b.JA.T(K+1,:).*b.XI_z.T(K+1,:).*b.w.T(K+1,:)-b.JA.T(K,:).*b.XI_z.T(K,:).*b.w.T(K,:))./[b.DXI.rm(K,1) b.DXI.rm(K,:) b.DXI.rm(K,end)];

            % -dmu/dT*r_z*2*xi_r*d(w0/eta_z)/delta.xi
                T_C=T_C-dmu(b.T.u(K,:)).*R_z.u(K,:).*2.*b.XI_r.u(K,:).*(b.w.T(K+1,:)./b.ETA_z.T(K+1,:)-b.w.T(K,:)./b.ETA_z.T(K,:))./[b.DXI.rm(K,1) b.DXI.rm(K,:) b.DXI.rm(K,end)];

            T_C=T_C(1,2:end-1);

            % -dmu/dT/N^2*2*r_z*eta_z*d(u0/xi_r)/delta.eta
                T_C=T_C-dmu(b.T.u(K,2:end-1))./N.u(1,2:end-1).^2.*2.*R_z.u(K,2:end-1).*b.ETA_z.u(K,2:end-1).*(b.u.w(K+(o-1),2:end)./b.XI_r.w(K+(o-1),2:end)-b.u.w(K+(o-1),1:end-1)./b.XI_r.w(K+(o-1),1:end-1))./b.DETA.c(1,:);

            % -dmu/N^2*2*r_z*w0*d(r_z)/delta.eta
                T_C=T_C-dmu(b.T.u(K,2:end-1))./N.u(1,2:end-1).^2.*2.*R_z.u(K,2:end-1).*b.w.u(K,2:end-1).*(R_z.w(K+(o-1),2:end)-R_z.w(K+(o-1),1:end-1))./b.DETA.c(1,:);

            % -2/3*dmu/dT/(J*Rs)*d(J*R*w0))/delta.eta
                T_C=T_C-2/3*dmu(b.T.u(K,2:end-1))./(b.JA.u(K,2:end-1).*b.R.u(K,2:end-1).^ax).*(b.JA.w(K+(o-1),2:end).*b.R.w(K+(o-1),2:end).^ax.*b.w.w(K+(o-1),2:end)-b.JA.w(K+(o-1),1:end-1).*b.R.w(K+(o-1),1:end-1).^ax.*b.w.w(K+(o-1),1:end-1))./b.DETA.c(1,:);

            % Laplace pressure
                if isfield(b, 'sigma')==1 && (sigma(30)-sigma(10))~=0
                    T_C=T_C+(-1)^o*dsigma(b.T.u(K,2:end-1)).*nabla_dot_n(1,2:end-1);

                end

            % hydrostatic pressure
                T_C=T_C-(-1)^o*g*drho(b.T.u(K,2:end-1)).*b.z.u(1,2:end-1);
                
            T_C=T_C.*b.dz.c.*N.u(1,2:end-1);

    end

    for p=1:2
        if max(strncmp('ss', {b.bc.z{(p-1)*length(b.bc.z)/2+1:(p-1)*length(b.bc.z)/2+length(b.bc.z)/2}},2))
            if p==1
                % derivative with respect to lower boundary
                    R_rlb.u=b.R_rlb1.u;
                    R_rlb.w=b.R_rlb1.w;
                    R_z_rlb.u=b.R1_z_rlb1.u;
                    R_z_rlb.w=b.R1_z_rlb1.w;
                    XI_r_rlb.u=b.XI_r_rlb1.u;
                    XI_r_rlb.w=b.XI_r_rlb1.w;
                    XI_r_rlb.T=b.XI_r_rlb1.T;
                    XI_r_inv_rlb.u=b.XI_r_inv_rlb1.u;
                    XI_r_inv_rlb.w=b.XI_r_inv_rlb1.w;
                    XI_z_rlb_C.u=b.XI_z_rlb1_C.u;
                    XI_z_rlb_C.w=b.XI_z_rlb1_C.w;
                    XI_z_rlb_C.T=b.XI_z_rlb1_C.T;
                    XI_z_rlb_E.u=b.XI_z_rlb1_E.u;
                    XI_z_rlb_E.w=b.XI_z_rlb1_E.w;
                    XI_z_rlb_E.T=b.XI_z_rlb1_E.T;
                    JA_rlb.u=b.JA_rlb1.u;
                    JA_rlb.w=b.JA_rlb1.w;
                    JA_rlb.T=b.JA_rlb1.T;

            elseif p==2
                % derivative with respect to upper boundary
                    R_rlb.u=b.R_rlbend.u;
                    R_rlb.w=b.R_rlbend.w;
                    R_z_rlb.u=b.Rend_z_rlbend.u;
                    R_z_rlb.w=b.Rend_z_rlbend.w;
                    XI_r_rlb.u=b.XI_r_rlbend.u;
                    XI_r_rlb.w=b.XI_r_rlbend.w;
                    XI_r_rlb.T=b.XI_r_rlbend.T;
                    XI_r_inv_rlb.u=b.XI_r_inv_rlbend.u;
                    XI_r_inv_rlb.w=b.XI_r_inv_rlbend.w;
                    XI_z_rlb_C.u=b.XI_z_rlbend_C.u;
                    XI_z_rlb_C.w=b.XI_z_rlbend_C.w;
                    XI_z_rlb_C.T=b.XI_z_rlbend_C.T;
                    XI_z_rlb_E.u=b.XI_z_rlbend_E.u;
                    XI_z_rlb_E.w=b.XI_z_rlbend_E.w;
                    XI_z_rlb_E.T=b.XI_z_rlbend_E.T;
                    JA_rlb.u=b.JA_rlbend.u;
                    JA_rlb.w=b.JA_rlbend.w;
                    JA_rlb.T=b.JA_rlbend.T;

            end

            % derivative of normal stresses with respect to Rs
                if o==p
                    %N_inv_sq_rlb.u=-1./N.u.^4.*2.*R_z.u(1,:).*R_z_rlb.u(1,:);
                    N_inv_rlb.u=-0.5./N.u.^3.*2.*R_z.u(1,:).*R_z_rlb.u(1,:);
                    N_rlb.u=0.5./N.u.*2.*R_z.u(1,:).*R_z_rlb.u(1,:);
                    
                    % mu*2*N*xi_r*d(u0/xi_r)/delta.xi
                        % mu*2*dxi_r/dRs*d(u0/xi_r)/delta.xi
                            FC.u=mu(b.T.u(K,:)).*2.*N.u.*XI_r_rlb.u(K,:).*(b.u.u(K+2-o,:)./b.XI_r.u(K+2-o,:)-b.u.u(K-(o-1),:)./b.XI_r.u(K-(o-1),:))./[b.DXI.c(K-(o-1),1) b.DXI.c(K-(o-1),:) b.DXI.c(K-(o-1),end)];

                        % mu*2*N*xi_r*d(u0*d(1/xi_r)/dRs)/delta.xi
                            FC.u=FC.u+mu(b.T.u(K,:)).*2.*N.u.*b.XI_r.u(K,:).*(b.u.u(K+2-o,:).*XI_r_inv_rlb.u(K+2-o,:)-b.u.u(K-(o-1),:).*XI_r_inv_rlb.u(K-(o-1),:))./[b.DXI.c(K-(o-1),1) b.DXI.c(K-(o-1),:) b.DXI.c(K-(o-1),end)];

                        % mu*2*dN/dRs*xi_r*d(u0/xi_r)/delta.xi
                            FE.u=mu(b.T.u(K,:)).*2.*N_rlb.u.*b.XI_r.u(K,:).*(b.u.u(K+2-o,:)./b.XI_r.u(K+2-o,:)-b.u.u(K-(o-1),:)./b.XI_r.u(K-(o-1),:))./[b.DXI.c(K-(o-1),1) b.DXI.c(K-(o-1),:) b.DXI.c(K-(o-1),end)];
                        
                    % -2/3*N*mu/(J*Rs)*d(J*R*u0))/delta.xi
                        % -2/3*N*mu*(-(dJ/dRs*Rs+J)/(J*Rs)^2*d(J*R*u0))/delta.xi
                            FC.u=FC.u-2/3.*N.u.*mu(b.T.u(K,:)).*(-(JA_rlb.u(K,:).*b.R.u(K,:).^ax+ax*b.JA.u(K,:).*R_rlb.u(K,:))./(b.JA.u(K,:).*b.R.u(K,:).^ax).^2).*(b.JA.u(K+2-o,:).*b.R.u(K+2-o,:).^ax.*b.u.u(K+2-o,:)-b.JA.u(K-(o-1),:).*b.R.u(K-(o-1),:).^ax.*b.u.u(K-(o-1),:))./[b.DXI.c(K-(o-1),1) b.DXI.c(K-(o-1),:) b.DXI.c(K-(o-1),end)];

                        % -2/3*N*mu/(J*Rs)*d((dJ/dRs*R+J*dR/dRs)*u0)/delta.xi
                            FC.u=FC.u-2/3.*N.u.*mu(b.T.u(K,:))./(b.JA.u(K,:).*b.R.u(K,:).^ax).*((JA_rlb.u(K+2-o,:).*b.R.u(K+2-o,:).^ax+ax*b.JA.u(K+2-o,:).*R_rlb.u(K+2-o,:)).*b.u.u(K+2-o,:)-(JA_rlb.u(K-(o-1),:).*b.R.u(K-(o-1),:).^ax+ax*b.JA.u(K-(o-1),:).*R_rlb.u(K-(o-1),:)).*b.u.u(K-(o-1),:))./[b.DXI.c(K-(o-1),1) b.DXI.c(K-(o-1),:) b.DXI.c(K-(o-1),end)];

                        % -2/3*dN/dRs*mu*/(J*Rs)*d(J*R*u0))/delta.xi
                            FE.u=FE.u-2/3.*N_rlb.u.*mu(b.T.u(K,:))./(b.JA.u(K,:).*b.R.u(K,:).^ax).*(b.JA.u(K+2-o,:).*b.R.u(K+2-o,:).^ax.*b.u.u(K+2-o,:)-b.JA.u(K-(o-1),:).*b.R.u(K-(o-1),:).^ax.*b.u.u(K-(o-1),:))./[b.DXI.c(K-(o-1),1) b.DXI.c(K-(o-1),:) b.DXI.c(K-(o-1),end)];

                    % -mu*2*N*xi_r*d(J*xi_z*w0)/delta.xi
                        % -mu*2*N*dxi_r/dRs*d(J*xi_z*w0)/delta.xi
                            FC.u=FC.u-mu(b.T.u(K,:)).*2.*N.u.*XI_r_rlb.u(K,:).*(b.JA.T(K+1,:).*b.XI_z.T(K+1,:).*b.w.T(K+1,:)-b.JA.T(K,:).*b.XI_z.T(K,:).*b.w.T(K,:))./[b.DXI.rm(K,1) b.DXI.rm(K,:) b.DXI.rm(K,end)];

                        % -mu*2*N*xi_r*d((dJ/dRs*xi_z+J*dxi_z/dRs)*w0)/delta.xi
                            FC.u=FC.u-mu(b.T.u(K,:)).*2.*N.u.*b.XI_r.u(K,:).*((JA_rlb.T(K+1,:).*b.XI_z.T(K+1,:)+b.JA.T(K+1,:).*XI_z_rlb_C.T(K+1,:)).*b.w.T(K+1,:)-(JA_rlb.T(K,:).*b.XI_z.T(K,:)+b.JA.T(K,:).*XI_z_rlb_C.T(K,:)).*b.w.T(K,:))./[b.DXI.rm(K,1) b.DXI.rm(K,:) b.DXI.rm(K,end)];
                            FE.u=FE.u-mu(b.T.u(K,:)).*2.*N.u.*b.XI_r.u(K,:).*(b.JA.T(K+1,:).*XI_z_rlb_E.T(K+1,:).*b.w.T(K+1,:)-b.JA.T(K,:).*XI_z_rlb_E.T(K,:).*b.w.T(K,:))./[b.DXI.rm(K,1) b.DXI.rm(K,:) b.DXI.rm(K,end)];
                            
                        % -mu*2*dN/dRs*xi_r*d(J*xi_z*w0)/delta.xi
                            FE.u=FE.u-mu(b.T.u(K,:)).*2.*N_rlb.u.*b.XI_r.u(K,:).*(b.JA.T(K+1,:).*b.XI_z.T(K+1,:).*b.w.T(K+1,:)-b.JA.T(K,:).*b.XI_z.T(K,:).*b.w.T(K,:))./[b.DXI.rm(K,1) b.DXI.rm(K,:) b.DXI.rm(K,end)];
                        
                    % -mu*r_z*2*N*xi_r*d(w0/eta_z)/delta.xi*delta.z
                        % -mu*r_z*2*N*dxi_r/dRs*d(w0/eta_z)/delta.xi
                            F.T=(b.w.T(K+1,:)./b.ETA_z.T(K+1,:)-b.w.T(K,:)./b.ETA_z.T(K,:))./[b.DXI.rm(K,1) b.DXI.rm(K,:) b.DXI.rm(K,end)];
                            FC.u=FC.u-mu(b.T.u(K,:)).*R_z.u(1,:).*2.*N.u.*XI_r_rlb.u(K,:).*F.T;

                        % -mu*d(r_z*N)/dRs*2*xi_r*d(w0/eta_z)/delta.xi
                            FE.u=FE.u-mu(b.T.u(K,:)).*(R_z_rlb.u(1,:).*N.u + R_z.u(1,:).*N_rlb.u).*2.*b.XI_r.u(K,:).*F.T;

                    FC.u=FC.u(1,2:end-1);
                    FE.u=FE.u(1,2:end-1);

                    % -mu/N*2*r_z*eta_z*d(u0/xi_r)/delta.eta
                        % -mu/N^2*2*r_z*eta_z*d(u0*d(1/xi_r)/dRs)/delta.eta
                            F.u=-mu(b.T.u(K,2:end-1))./N.u(1,2:end-1).*2.*R_z.u(1,2:end-1).*b.ETA_z.u(K,2:end-1);
                            F_W.u=F.u.*(b.u.w(K+(o-1),1:end-1).*XI_r_inv_rlb.w(K+(o-1),1:end-1))./b.DETA.c(1,:);
                            F_E.u=F.u.*(b.u.w(K+(o-1),2:end).*XI_r_inv_rlb.w(K+(o-1),2:end))./b.DETA.c(1,:);

                        % -mu*(d(1/N)/dRs*r_z+1/N*dr_z/dRs)*2*eta_z*d(u0/xi_r)/delta.eta
                            F.u=-mu(b.T.u(K,2:end-1)).*(N_inv_rlb.u(1,2:end-1).*R_z.u(1,2:end-1)+1./N.u(1,2:end-1).*R_z_rlb.u(1,2:end-1)).*2.*b.ETA_z.u(K,2:end-1).*(b.u.w(K+(o-1),2:end).*b.XI_r.w(K+(o-1),2:end)-b.u.w(K+(o-1),1:end-1).*b.XI_r.w(K+(o-1),1:end-1))./b.DETA.c(1,:);
                            FE.u=FE.u+F.u;

                    % -mu/N*2*r_z*w0*d(r_z)/delta.eta
                        % -mu/N*2*r_z*w0*d(dr_z/dRs)/delta.eta
                            F.u=-mu(b.T.u(K,2:end-1))./N.u(1,2:end-1).*2.*R_z.u(1,2:end-1).*b.w.u(K,2:end-1);
                            FE_W.u=F.u.*R_z_rlb.w(K+(o-1),1:end-1)./b.DETA.c(1,:);
                            FE_E.u=F.u.*R_z_rlb.w(K+(o-1),2:end)./b.DETA.c(1,:);

                        % -mu(d(1/N)/dRs*r_z+1/N*dr_z/dRs)*2*w0*d(r_z)/delta.eta
                            FE.u=FE.u-mu(b.T.u(K,2:end-1)).*(N_inv_rlb.u(1,2:end-1).*R_z.u(1,2:end-1)+1./N.u(1,2:end-1).*R_z_rlb.u(1,2:end-1)).*2.*b.w.u(K,2:end-1).*(R_z.w(K+(o-1),2:end)-R_z.w(K+(o-1),1:end-1))./b.DETA.c(1,:);
                                                    
                    % -2/3*mu*N*/(J*Rs)*d(J*R*w0))/delta.eta
                        % -2/3*mu*N*(-(dJ/dRs*Rs+J)/(J*Rs)^2*d(J*R*w0))/delta.eta
                            FC.u=FC.u-2/3*mu(b.T.u(K,2:end-1)).*N.u(1,2:end-1).*(-(JA_rlb.u(K,2:end-1).*b.R.u(K,2:end-1).^ax+ax*b.JA.u(K,2:end-1).*R_rlb.u(K,2:end-1))./(b.JA.u(K,2:end-1).*b.R.u(K,2:end-1).^ax).^2).*(b.JA.w(K+(o-1),2:end).*b.R.w(K+(o-1),2:end).^ax.*b.w.w(K+(o-1),2:end)-b.JA.w(K+(o-1),1:end-1).*b.R.w(K+(o-1),1:end-1).^ax.*b.w.w(K+(o-1),1:end-1))./b.DETA.c(1,:);

                        % -2/3*mu*N/(J*Rs)*d((dJ/dRs*R+J*dR/dRs)*w0)/delta.eta
                            F.u=-2/3*mu(b.T.u(K,2:end-1)).*N.u(1,2:end-1)./(b.JA.u(K,2:end-1).*b.R.u(K,2:end-1).^ax);
                            F_W.u=F_W.u+F.u.*(JA_rlb.w(K+(o-1),1:end-1).*b.R.w(K+(o-1),1:end-1).^ax+ax*b.JA.w(K+(o-1),1:end-1).*R_rlb.w(K+(o-1),1:end-1)).*b.w.w(K+(o-1),1:end-1)./b.DETA.c(1,:);
                            F_E.u=F_E.u+F.u.*(JA_rlb.w(K+(o-1),2:end).*b.R.w(K+(o-1),2:end).^ax+ax*b.JA.w(K+(o-1),2:end).*R_rlb.w(K+(o-1),2:end)).*b.w.w(K+(o-1),2:end)./b.DETA.c(1,:);

                        % -2/3*mu*dN/dRs/(J*Rs)*d(J*R*w0))/delta.eta
                            FE.u=FE.u-2/3*mu(b.T.u(K,2:end-1)).*N_rlb.u(1,2:end-1)./(b.JA.u(K,2:end-1).*b.R.u(K,2:end-1).^ax).*(b.JA.w(K+(o-1),2:end).*b.R.w(K+(o-1),2:end).^ax.*b.w.w(K+(o-1),2:end)-b.JA.w(K+(o-1),1:end-1).*b.R.w(K+(o-1),1:end-1).^ax.*b.w.w(K+(o-1),1:end-1))./b.DETA.c(1,:);

                    % Laplace pressure
                        if isfield(b, 'sigma')==1
                            FC.u=FC.u-(-1)^o*sigma(b.T.u(K,2:end-1)).*ax.*(b.R.u(K,2:end-1).^ax).^(-2).*R_rlb.u(K,2:end-1);
                            %%{
                            FE_W.u=FE_W.u-(-1)^o*sigma(b.T.u(K,2:end-1))./b.dz.zm(1:end-1)./b.dz.c./N.u(1,2:end-1).^2;
                            FE_E.u=FE_E.u-(-1)^o*sigma(b.T.u(K,2:end-1))./b.dz.zm(2:end)./b.dz.c./N.u(1,2:end-1).^2;
                            %}
                            %{
                            FE_W.u=FE_W.u-sigma(b.T.u(K,2:end-1)).*b.eta_z.w(1:end-1)./b.deta.zm(1:end-1).*b.eta_z.u(2:end-1)./b.deta.c./N.u(1,2:end-1).^3;
                            FE_E.u=FE_E.u-sigma(b.T.u(K,2:end-1)).*b.eta_z.w(2:end)./b.deta.zm(2:end).*b.eta_z.u(2:end-1)./b.deta.c./N.u(1,2:end-1).^3;
                            %}
                            FE.u=FE.u+(-1)^o*sigma(b.T.u(K,2:end-1)).*(R_zz.u(1,2:end-1).*1./N.u(1,2:end-1).^4.*2.*R_z.u(1,2:end-1).*R_z_rlb.u(1,2:end-1));
                            
                        end
                        
                     % hydrostatic pressure
                        FE.u=FE.u-(-1)^o*g*rho(b.T.u(K,2:end-1)).*b.z.u(1,2:end-1).*N_rlb.u(1,2:end-1);
                        
                    % -dN/dRs*p
                        FE.u=FE.u+(-1.5.*b.p.p(K-(o-1),:)+0.5.*b.p.p(K-o*(-1)^o,:)).*N_rlb.u(1,2:end-1);

                else
                    %{
                    FC.u=-(-1)^o*mu(b.T.u(K,:))./N.u.^2.*2.*(b.XI_r.u(K,:)-R_z.u(1,:).*b.XI_z.u(K,:)).*JA_rlb.u(K-(-1)^o,:).*b.ETA_z.u(K-(-1)^o,:).*b.u.u(K-(-1)^o,:)./[b.DXI.c(K-(o-1),1) b.DXI.c(K-(o-1),:) b.DXI.c(K-(o-1),end)];
                    FC.u=FC.u+(-1)^o*2/3*mu(b.T.u(K,:))./(b.JA.u(K,:).*b.R.u(K,:).^ax).*(JA_rlb.u(K-(-1)^o,:).*b.R.u(K-(-1)^o,:).^ax+ax*b.JA.u(K-(-1)^o,:).*R_rlb.u(K-(-1)^o,:)).*b.u.u(K-(-1)^o,:)./[b.DXI.c(K-(o-1),1) b.DXI.c(K-(o-1),:) b.DXI.c(K-(o-1),end)];
                    FC.u=FC.u+(-1)^o*mu(b.T.u(K,:))./N.u.^2.*2.*(b.XI_r.u(K,:)-R_z.u(1,:).*b.XI_z.u(K,:)).*(JA_rlb.T(K+(2-o),:).*b.XI_z.T(K+(2-o),:)+b.JA.T(K+(2-o),:).*XI_z_rlb_C.T(K+(2-o),:)).*b.w.T(K+(2-o),:)./[b.DXI.rm(K,1) b.DXI.rm(K,:) b.DXI.rm(K,end)];
                    FE.u=(-1)^o*mu(b.T.u(K,:))./N.u.^2.*2.*(b.XI_r.u(K,:)-R_z.u(1,:).*b.XI_z.u(K,:)).*b.JA.T(K+(2-o),:).*XI_z_rlb_E.T(K+(2-o),:).*b.w.T(K+(2-o),:)./[b.DXI.rm(K,1) b.DXI.rm(K,:) b.DXI.rm(K,end)];
                    FC.u=FC.u+(-1)^o*mu(b.T.u(K,:))./N.u.^2.*R_z.u(1,:).*2.*(b.XI_r.u(K,:)-R_z.u(1,:).*b.XI_z.u(K,:)).*(JA_rlb.T(K+(2-o),:).*b.XI_r.T(K+(2-o),:)+b.JA.T(K+(2-o),:).*XI_r_rlb.T(K+(2-o),:)).*b.w.T(K+(2-o),:)./[b.DXI.rm(K,1) b.DXI.rm(K,:) b.DXI.rm(K,end)];

                    FC.u=FC.u(1,2:end-1);
                    FE.u=FE.u(1,2:end-1);
                    %}
                    FC.u=sparse(1,b.J);
                    FE.u=sparse(1,b.J);
                    F_W.u=sparse(1,b.J);
                    F_E.u=sparse(1,b.J);
                    FE_W.u=sparse(1,b.J);
                    FE_E.u=sparse(1,b.J);

                end

            FC.u=FC.u.*b.dz.c;
            FE.u=FE.u.*b.dz.c;
            F_W.u=F_W.u.*b.dz.c;
            F_E.u=F_E.u.*b.dz.c;
            FE_W.u=FE_W.u.*b.dz.c;
            FE_E.u=FE_E.u.*b.dz.c;
                    
            Rs_W = -[0 b.W.w(1,2:end)].*(FE.u + F_W.u) + [0 FE_W.u(1,2:end)];
            Rs_E = [b.E.w(1,1:end-1) 0].*(FE.u + F_E.u) + [FE_E.u(1,1:end-1) 0];
            Rs_ad_W = - [0 b.E.w(1,1:end-1)].*(FE.u + F_W.u);
            Rs_ad_E = [b.W.w(1,2:end) 0].*(FE.u + F_E.u);
            Rs_C = FC.u + Rs_ad_W + Rs_ad_E - FE_W.u - FE_E.u;

            Ja.ns.Rs=sparse(b.J+2,b.J+2);
            Ja.ns.Rs(2:end-1,:)=spdiags([Rs_W.' Rs_C.' Rs_E.'],[0 1 2],b.J,b.J+2);

            if (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'s')==1
                Ja.ns.Rs(1,1)=XI_z_rlb_E.u(K,1).*(b.u.u(K+2-o,1)-b.u.u(K-(o-1),1))./b.DXI.c(K-(o-1),1);

            end
            if (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'s')==1
                Ja.ns.Rs(end,end)=XI_z_rlb_E.u(K,end).*(b.u.u(K+2-o,end)-b.u.u(K-(o-1),end))./b.DXI.c(K-(o-1),end);

            end

            if o==1 && o==p
                b.ns1.u=spdiags([u_Wsw.' u_s.' u_Ese.' u_n.'],[-1 0 1 b.J+2],(b.J+2),(b.K+1)*(b.J+2));
                b.ns1.w=spdiags([w_sw.' w_se.' w_nw.' w_ne.'],[-1 0 b.J b.J+1],(b.J+2),(b.K+2)*(b.J+1));
                b.ns1.p=spdiags([p_C.' p_N.'],[-1 b.J-1],(b.J+2),b.K*b.J);
                if isfield(b, 'sigma')==1
                    LP1=LP;

                end
                HP1=HP;
                % Jacobi of normal stresses with respect to T
                    if energy==1
                        Ja.ns.T1=sparse(b.J+2,(b.K+2)*(b.J+2));
                        Ja.ns.T1(2:end-1,:)=spdiags(T_C.',1,b.J,(b.K+2)*(b.J+2));

                    end
                % Jacobi of normal stresses with respect to Rs
                    Ja.ns.Rs1=Ja.ns.Rs;
                    Ja.ns.Rs1_ad_W=Rs_ad_W;
                    Ja.ns.Rs1_ad_E=Rs_ad_E;

            elseif o==2 && o==p
                b.nsend.u=spdiags([u_n.' u_Wsw.' u_s.' u_Ese.'],(b.K-1)*(b.J+2)+[0 b.J+1 b.J+2 b.J+3],(b.J+2),(b.K+1)*(b.J+2));
                b.nsend.w=spdiags([w_nw.' w_ne.' w_sw.' w_se.'],(b.K)*(b.J+1)+[-1 0 b.J b.J+1],(b.J+2),(b.K+2)*(b.J+1));
                b.nsend.p=spdiags([p_N.' p_C.'],(b.K-2)*(b.J)+[-1 b.J-1],(b.J+2),b.K*b.J);
                if isfield(b, 'sigma')==1
                    LPend=LP;

                end
                HPend=HP;
                % Jacobi of normal stresses with respect to T
                    if energy==1
                        Ja.ns.Tend=sparse(b.J+2,(b.K+2)*(b.J+2));
                        Ja.ns.Tend(2:end-1,:)=spdiags(T_C.',(b.K+1)*(b.J+2)+1,b.J,(b.K+2)*(b.J+2));

                    end
                % Jacobi of normal stresses with respect to Rs
                    Ja.ns.Rsend=Ja.ns.Rs;
                    Ja.ns.Rsend_ad_W=Rs_ad_W;
                    Ja.ns.Rsend_ad_E=Rs_ad_E;

            end

        end

    end

end
if o==1
    % equal velocity
        b.e1.u=spdiags(equal_u.',0,b.J+2,(b.K+1)*(b.J+2));

else
    % equal velocity
        b.eend.u=[sparse(b.J+2, b.K*(b.J+2)) spdiags(equal_u.',0,b.J+2,b.J+2)];

end