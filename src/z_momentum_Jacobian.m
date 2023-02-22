function[b]=z_momentum_Jacobian(b, flowopt)
    ax=flowopt.ax;
    g=flowopt.g;
    energy=flowopt.energy;
    stability=flowopt.stability;
    inviscid=flowopt.inviscid;
    creeping=flowopt.creeping;
    
    rho=str2func(b.rho);
    mu=str2func(b.mu);
    
    drho=str2func(b.drho);
    dmu=str2func(b.dmu);
    
    b.Ja.z_m.u=sparse((b.K+2)*(b.J+1),(b.K+1)*(b.J+2));
    b.Ja.z_m.w=sparse((b.K+2)*(b.J+1),(b.K+2)*(b.J+1));
    
    if creeping~=1
    %  [J^2*r^ax*rho0*xi_r*w0]_delta.xi*delta.eta
        F.v=b.JA.v.^2.*b.R.v.^ax.*rho(b.T.v).*b.XI_r.v.*b.w.v;

        u_sw=zeros(b.K+2,b.J+2); u_se=u_sw; u_nw=u_sw; u_ne=u_sw;

        u_sw(1:end-2,2:end-2)   = -F.v(1:end-1,2:end-1)/2.*b.DETA.zm(:,2:end-1);
        u_se(1:end-2,3:end-1)   = -F.v(1:end-1,2:end-1)/2.*b.DETA.zm(:,2:end-1);
        u_nw(2:end-1,2:end-2)   = F.v(2:end,2:end-1)/2.*b.DETA.zm(:,2:end-1);
        u_ne(2:end-1,3:end-1)   = F.v(2:end,2:end-1)/2.*b.DETA.zm(:,2:end-1);

        u_sw=u_sw.'; u_sw=u_sw(:);
        u_se=u_se.'; u_se=u_se(:);
        u_nw=u_nw.'; u_nw=u_nw(:);
        u_ne=u_ne.'; u_ne=u_ne(:);

        b.Ja.z_m.u=spdiags([u_sw u_se u_nw u_ne],[-(b.J+3) -(b.J+2) -1 0],(b.K+2)*(b.J+2),(b.K+1)*(b.J+2));
        b.Ja.z_m.u(1:b.J+2:end,:)=[];

    %  [J^2*r^ax*rho0*xi_r*u0]_delta.xi*delta.eta +[J^2*r^ax*rho0*xi_r*w0*2]_delta.eta*delta.xi
        F.p=b.JA.p.^2.*b.R.p.^ax.*rho(b.T.p).*b.XI_r.p.*b.w.p*2;
        F.v=b.JA.v.^2.*b.R.v.^ax.*rho(b.T.v).*b.XI_r.v.*b.u.v;

        w_S=zeros(b.K+2,b.J+1); w_W=w_S; w_C=w_S; w_E=w_S; w_N=w_S;

        w_S(1,2:end-1)       = -F.v(1,2:end-1).*b.DETA.zm(1,2:end-1);
        w_S(2:end-2,2:end-1) = -F.v(2:end-1,2:end-1)/2.*b.DETA.zm(2:end,2:end-1);
        w_W(2:end-1,1:end-2) = -F.p(:,1:end-1)/2.*b.DXI.zm(:,2:end-1);
        w_E(2:end-1,3:end)   = F.p(:,2:end)/2.*b.DXI.zm(:,2:end-1);
        w_N(end,2:end-1)     = F.v(end,2:end-1).*b.DETA.zm(end,2:end-1);
        w_N(3:end-1,2:end-1) = F.v(2:end-1,2:end-1)/2.*b.DETA.zm(1:end-1,2:end-1);
        w_C(3:end-1,2:end-1) = w_C(3:end-1,2:end-1)+w_S(2:end-2,2:end-1);
        w_C(2:end-1,2:end-1) = w_C(2:end-1,2:end-1)+w_W(2:end-1,1:end-2);
        w_C(2:end-1,2:end-1) = w_C(2:end-1,2:end-1)+w_E(2:end-1,3:end);
        w_C(2:end-2,2:end-1) = w_C(2:end-2,2:end-1)+w_N(3:end-1,2:end-1);

        w_S=w_S.'; w_S=w_S(:);
        w_W=w_W.'; w_W=w_W(:);
        w_C=w_C.'; w_C=w_C(:);
        w_E=w_E.'; w_E=w_E(:);
        w_N=w_N.'; w_N=w_N(:);

        b.Ja.z_m.w=spdiags([w_S w_W w_C w_E w_N],[-b.J-1 -1 0 1 b.J+1],(b.K+2)*(b.J+1),(b.K+2)*(b.J+1));
        
    end
        
    if inviscid~=1
        b.Ja.z_m.u=b.Ja.z_m.u + b.z_m_visc.u;
        b.Ja.z_m.w=b.Ja.z_m.w + b.z_m_visc.w;
        
    end
        
    if energy==1
        F.p=zeros(size(b.T.p)); F.v=zeros(size(b.T.v)); F1.v=F.v; F2.v=F.v; F3.v=F.v; F4.v=F.v; F5.v=F.v; F6.v=F.v; F7.v=F.v; F8.v=F.v; F1.p=F.p; F2.p=F.p; F3.p=F.p;
        
        if creeping~=1
        %  [J^2*r^ax*drho0/dT0*xi_r*u0*w0]_delta.xi*delta.eta +[J^2*r^ax*drho0/dT0*xi_r*w0^2]_delta.eta*delta.xi
            F.p=b.JA.p.^2.*b.R.p.^ax.*drho(b.T.p).*b.XI_r.p.*b.w.p.^2;
            F.v=b.JA.v.^2.*b.R.v.^ax.*drho(b.T.v).*b.XI_r.v.*b.u.v.*b.w.v;
            
        end

        if inviscid~=1
        % viscous terms for T
        % sign changed because of viscous terms!!!
            F1.v=-dmu(b.T.v).*b.JA.v.*b.XI_r.v.*b.XI_z.v.*b.R.v.^ax;
            F2.v=2/3*dmu(b.T.v).*b.XI_z.v;
            F3.v=-dmu(b.T.v).*b.JA.v.*b.XI_r.v.*b.ETA_z.v.*b.R.v.^ax;
            F4.v=-dmu(b.T.v).*b.JA.v.*(b.XI_r.v.^2+2*b.XI_z.v.^2).*b.R.v.^ax;
            F5.v=dmu(b.T.v).*b.JA.v.*b.XI_r.v.*b.XI_z.v.*b.R.v.^ax;
            F6.v=2/3*dmu(b.T.v).*b.XI_z.v;
            F7.v=-2*dmu(b.T.v).*b.JA.v.*b.XI_z.v.*b.ETA_z.v.*b.R.v.^ax;
            F8.v=dmu(b.T.v).*b.JA.v.*b.XI_r.v.*b.ETA_z.v.*b.R.v.^ax;
            
            F1.p=2/3*dmu(b.T.p).*b.ETA_z.p;
            F2.p=-2*dmu(b.T.p).*b.JA.p.*b.XI_z.p.*b.ETA_z.p.*b.R.p.^ax;
            F3.p=-2*dmu(b.T.p).*b.JA.p.*b.ETA_z.p.^2.*b.R.p.^ax;
            
        end
        
        if isfield(flowopt,'boussinesq')==1 && flowopt.boussinesq==1
            drho = @(theta) -rho(b.T_0).*b.beta;
        end
            
        F_s.v=F2.v(1:end-1,2:end-1).*(b.JA.w(2:end-1,2:end-1).*b.R.w(2:end-1,2:end-1).^ax.*b.u.w(2:end-1,2:end-1) - b.JA.w(1:end-2,2:end-1).*b.R.w(1:end-2,2:end-1).^ax.*b.u.w(1:end-2,2:end-1))./b.DXI.rm(1:end-1,1:end-1);
        F_s.v=F_s.v + F1.v(1:end-1,2:end-1).*(b.JA.w(2:end-1,2:end-1).*b.ETA_z.w(2:end-1,2:end-1).*b.u.w(2:end-1,2:end-1) - b.JA.w(1:end-2,2:end-1).*b.ETA_z.w(1:end-2,2:end-1).*b.u.w(1:end-2,2:end-1))./b.DXI.rm(1:end-1,1:end-1);
        F_s.v=F_s.v + F3.v(1:end-1,2:end-1).*(b.JA.u(1:end-1,3:end-1).*b.ETA_z.u(1:end-1,3:end-1).*b.u.u(1:end-1,3:end-1) - b.JA.u(1:end-1,2:end-2).*b.ETA_z.u(1:end-1,2:end-2).*b.u.u(1:end-1,2:end-2))./b.DETA.zm(:,2:end-1);
        F_s.v=F_s.v + F4.v(1:end-1,2:end-1).*(b.JA.w(2:end-1,2:end-1).*b.XI_r.w(2:end-1,2:end-1).*b.w.w(2:end-1,2:end-1) - b.JA.w(1:end-2,2:end-1).*b.XI_r.w(1:end-2,2:end-1).*b.w.w(1:end-2,2:end-1))./b.DXI.rm(1:end-1,1:end-1);
        F_s.v=F_s.v + F5.v(1:end-1,2:end-1).*(b.JA.w(2:end-1,2:end-1).*b.XI_z.w(2:end-1,2:end-1).*b.w.w(2:end-1,2:end-1) - b.JA.w(1:end-2,2:end-1).*b.XI_z.w(1:end-2,2:end-1).*b.w.w(1:end-2,2:end-1))./b.DXI.rm(1:end-1,1:end-1);
        F_s.v=F_s.v + F6.v(1:end-1,2:end-1).*(b.JA.u(1:end-1,3:end-1).*b.R.u(1:end-1,3:end-1).^ax.*b.w.u(1:end-1,3:end-1) - b.JA.u(1:end-1,2:end-2).*b.R.u(1:end-1,2:end-2).^ax.*b.w.u(1:end-1,2:end-2))./b.DETA.zm(:,2:end-1);
        F_s.v=F_s.v + F7.v(1:end-1,2:end-1).*(b.JA.u(1:end-1,3:end-1).*b.XI_r.u(1:end-1,3:end-1).*b.w.u(1:end-1,3:end-1) - b.JA.u(1:end-1,2:end-2).*b.XI_r.u(1:end-1,2:end-2).*b.w.u(1:end-1,2:end-2))./b.DETA.zm(:,2:end-1);
        F_s.v=F_s.v + F8.v(1:end-1,2:end-1).*(b.JA.u(1:end-1,3:end-1).*b.XI_z.u(1:end-1,3:end-1).*b.w.u(1:end-1,3:end-1) - b.JA.u(1:end-1,2:end-2).*b.XI_z.u(1:end-1,2:end-2).*b.w.u(1:end-1,2:end-2))./b.DETA.zm(:,2:end-1);

        F_n.v=F2.v(2:end,2:end-1).*(b.JA.w(3:end,2:end-1).*b.R.w(3:end,2:end-1).^ax.*b.u.w(3:end,2:end-1) - b.JA.w(2:end-1,2:end-1).*b.R.w(2:end-1,2:end-1).^ax.*b.u.w(2:end-1,2:end-1))./b.DXI.rm(2:end,1:end-1);
        F_n.v=F_n.v + F1.v(2:end,2:end-1).*(b.JA.w(3:end,2:end-1).*b.ETA_z.w(3:end,2:end-1).*b.u.w(3:end,2:end-1) - b.JA.w(2:end-1,2:end-1).*b.ETA_z.w(2:end-1,2:end-1).*b.u.w(2:end-1,2:end-1))./b.DXI.rm(2:end,1:end-1);
        F_n.v=F_n.v + F3.v(2:end,2:end-1).*(b.JA.u(2:end,3:end-1).*b.ETA_z.u(2:end,3:end-1).*b.u.u(2:end,3:end-1) - b.JA.u(2:end,2:end-2).*b.ETA_z.u(2:end,2:end-2).*b.u.u(2:end,2:end-2))./b.DETA.zm(:,2:end-1);
        F_n.v=F_n.v + F4.v(2:end,2:end-1).*(b.JA.w(3:end,2:end-1).*b.XI_r.w(3:end,2:end-1).*b.w.w(3:end,2:end-1) - b.JA.w(2:end-1,2:end-1).*b.XI_r.w(2:end-1,2:end-1).*b.w.w(2:end-1,2:end-1))./b.DXI.rm(2:end,1:end-1);
        F_n.v=F_n.v + F5.v(2:end,2:end-1).*(b.JA.w(3:end,2:end-1).*b.XI_z.w(3:end,2:end-1).*b.w.w(3:end,2:end-1) - b.JA.w(2:end-1,2:end-1).*b.XI_z.w(2:end-1,2:end-1).*b.w.w(2:end-1,2:end-1))./b.DXI.rm(2:end,1:end-1);
        F_n.v=F_n.v + F6.v(2:end,2:end-1).*(b.JA.u(2:end,3:end-1).*b.R.u(2:end,3:end-1).^ax.*b.w.u(2:end,3:end-1) - b.JA.u(2:end,2:end-2).*b.R.u(2:end,2:end-2).^ax.*b.w.u(2:end,2:end-2))./b.DETA.zm(:,2:end-1);
        F_n.v=F_n.v + F7.v(2:end,2:end-1).*(b.JA.u(2:end,3:end-1).*b.XI_r.u(2:end,3:end-1).*b.w.u(2:end,3:end-1) - b.JA.u(2:end,2:end-2).*b.XI_r.u(2:end,2:end-2).*b.w.u(2:end,2:end-2))./b.DETA.zm(:,2:end-1);
        F_n.v=F_n.v + F8.v(2:end,2:end-1).*(b.JA.u(2:end,3:end-1).*b.XI_z.u(2:end,3:end-1).*b.w.u(2:end,3:end-1) - b.JA.u(2:end,2:end-2).*b.XI_z.u(2:end,2:end-2).*b.w.u(2:end,2:end-2))./b.DETA.zm(:,2:end-1);

        F_w.p=F1.p(:,1:end-1).*(b.JA.u(2:end,2:end-2).*b.R.u(2:end,2:end-2).^ax.*b.u.u(2:end,2:end-2) - b.JA.u(1:end-1,2:end-2).*b.R.u(1:end-1,2:end-2).^ax.*b.u.u(1:end-1,2:end-2))./b.DXI.zm(:,2:end-1);
        F_w.p=F_w.p + F2.p(:,1:end-1).*(b.JA.u(2:end,2:end-2).*b.XI_r.u(2:end,2:end-2).*b.w.u(2:end,2:end-2) - b.JA.u(1:end-1,2:end-2).*b.XI_r.u(1:end-1,2:end-2).*b.w.u(1:end-1,2:end-2))./b.DXI.zm(:,2:end-1);
        F_w.p=F_w.p + F3.p(:,1:end-1).*(b.JA.w(2:end-1,2:end-1).*b.XI_r.w(2:end-1,2:end-1).*b.w.w(2:end-1,2:end-1) - b.JA.w(2:end-1,1:end-2).*b.XI_r.w(2:end-1,1:end-2).*b.w.w(2:end-1,1:end-2))./b.DETA.c(:,1:end-1);
        F_w.p=F_w.p + F1.p(:,1:end-1).*(b.JA.w(2:end-1,2:end-1).*b.R.w(2:end-1,2:end-1).^ax.*b.w.w(2:end-1,2:end-1) - b.JA.w(2:end-1,1:end-2).*b.R.w(2:end-1,1:end-2).^ax.*b.w.w(2:end-1,1:end-2))./b.DETA.c(:,1:end-1);
        F_w.p=F_w.p + 0.5*g*drho(b.T.w(2:end-1,2:end-1)).*b.R.w(2:end-1,2:end-1).^ax.*b.JA.w(2:end-1,2:end-1).*b.DETA.zm(:,2:end-1);

        F_e.p=F1.p(:,2:end).*(b.JA.u(2:end,3:end-1).*b.R.u(2:end,3:end-1).^ax.*b.u.u(2:end,3:end-1) - b.JA.u(1:end-1,3:end-1).*b.R.u(1:end-1,3:end-1).^ax.*b.u.u(1:end-1,3:end-1))./b.DXI.zm(:,2:end-1);
        F_e.p=F_e.p + F2.p(:,2:end).*(b.JA.u(2:end,3:end-1).*b.XI_r.u(2:end,3:end-1).*b.w.u(2:end,3:end-1) - b.JA.u(1:end-1,3:end-1).*b.XI_r.u(1:end-1,3:end-1).*b.w.u(1:end-1,3:end-1))./b.DXI.zm(:,2:end-1);
        F_e.p=F_e.p + F3.p(:,2:end).*(b.JA.w(2:end-1,3:end).*b.XI_r.w(2:end-1,3:end).*b.w.w(2:end-1,3:end) - b.JA.w(2:end-1,2:end-1).*b.XI_r.w(2:end-1,2:end-1).*b.w.w(2:end-1,2:end-1))./b.DETA.c(:,2:end);
        F_e.p=F_e.p + F1.p(:,2:end).*(b.JA.w(2:end-1,3:end).*b.R.w(2:end-1,3:end).^ax.*b.w.w(2:end-1,3:end) - b.JA.w(2:end-1,2:end-1).*b.R.w(2:end-1,2:end-1).^ax.*b.w.w(2:end-1,2:end-1))./b.DETA.c(:,2:end);
        F_e.p=F_e.p - 0.5*g*drho(b.T.w(2:end-1,2:end-1)).*b.R.w(2:end-1,2:end-1).^ax.*b.JA.w(2:end-1,2:end-1).*b.DETA.zm(:,2:end-1);

        T_Sse=zeros(b.K+2,b.J+2); T_Ssw=T_Sse; T_e=T_Sse; T_w=T_Sse; T_Nne=T_Sse; T_Nnw=T_Sse;

        T_Ssw(2:end-2,2:end-2)   = -(F.v(2:end-1,2:end-1) + F_s.v(2:end,:))/4.*b.DETA.zm(2:end,2:end-1);
        T_Sse(2:end-2,3:end-1)   = -(F.v(2:end-1,2:end-1) + F_s.v(2:end,:))/4.*b.DETA.zm(2:end,2:end-1);
        T_Nnw(3:end-1,2:end-2)   = (F.v(2:end-1,2:end-1) + F_n.v(1:end-1,:))/4.*b.DETA.zm(1:end-1,2:end-1);
        T_Nne(3:end-1,3:end-1)   = (F.v(2:end-1,2:end-1) + F_n.v(1:end-1,:))/4.*b.DETA.zm(1:end-1,2:end-1);
        T_w(2:end-1,2:end-2)     = T_Ssw(1:end-2,2:end-2)+T_Nnw(3:end,2:end-2)-(F.p(:,1:end-1) + F_w.p).*b.DXI.zm(:,2:end-1);
        T_e(2:end-1,3:end-1)     = T_Sse(1:end-2,3:end-1)+T_Nne(3:end,3:end-1)+(F.p(:,2:end) + F_e.p).*b.DXI.zm(:,2:end-1);

        T_Ssw(1,2:end-2)         = -(F.v(1,2:end-1) + F_s.v(1,:))/2.*b.DETA.zm(1,2:end-1);
        T_Sse(1,3:end-1)         = -(F.v(1,2:end-1) + F_s.v(1,:))/2.*b.DETA.zm(1,2:end-1);
        T_Nnw(end,2:end-2)       = (F.v(end,2:end-1) + F_n.v(end,:))/2.*b.DETA.zm(end,2:end-1);
        T_Nne(end,3:end-1)       = (F.v(end,2:end-1) + F_n.v(end,:))/2.*b.DETA.zm(end,2:end-1);

        T_Ssw=T_Ssw.'; T_Ssw=T_Ssw(:);
        T_Sse=T_Sse.'; T_Sse=T_Sse(:);
        T_w=T_w.'; T_w=T_w(:);
        T_e=T_e.'; T_e=T_e(:);
        T_Nnw=T_Nnw.'; T_Nnw=T_Nnw(:);
        T_Nne=T_Nne.'; T_Nne=T_Nne(:);

        b.Ja.z_m.T=spdiags([T_Ssw T_Sse T_w T_e T_Nnw T_Nne],[-(b.J+3) -(b.J+2) -1 0 b.J+1 b.J+2],(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
        b.Ja.z_m.T(1:b.J+2:end,:)=[];
        
        drho = str2func(b.drho);
            
    end
    
    if stability~=1
        index_bc_z=length(b.bc.z)/2;
    
        for o=1:2
            if max(strncmp('ss', {b.bc.z{(o-1)*index_bc_z+1:o*index_bc_z}},2))
                if o==1
                    % lower boundary
                        JA_rlb.u=b.JA_rlb1.u;
                        JA_rlb.w=b.JA_rlb1.w;
                        JA_rlb.p=b.JA_rlb1.p;
                        JA_rlb.v=b.JA_rlb1.v;

                        R_rlb.u=b.R_rlb1.u;
                        R_rlb.w=b.R_rlb1.w;
                        R_rlb.p=b.R_rlb1.p;
                        R_rlb.v=b.R_rlb1.v;

                        XI_r_rlb.u=b.XI_r_rlb1.u;
                        XI_r_rlb.w=b.XI_r_rlb1.w;
                        XI_r_rlb.p=b.XI_r_rlb1.p;
                        XI_r_rlb.v=b.XI_r_rlb1.v;

                        XI_z_rlb_C.u=b.XI_z_rlb1_C.u;
                        XI_z_rlb_C.w=b.XI_z_rlb1_C.w;
                        XI_z_rlb_C.p=b.XI_z_rlb1_C.p;
                        XI_z_rlb_C.v=b.XI_z_rlb1_C.v;

                        XI_z_rlb_E.u=b.XI_z_rlb1_E.u;
                        XI_z_rlb_E.w=b.XI_z_rlb1_E.w;
                        XI_z_rlb_E.p=b.XI_z_rlb1_E.p;
                        XI_z_rlb_E.v=b.XI_z_rlb1_E.v;

                elseif o==2
                    % upper boundary
                        JA_rlb.u=b.JA_rlbend.u;
                        JA_rlb.w=b.JA_rlbend.w;
                        JA_rlb.p=b.JA_rlbend.p;
                        JA_rlb.v=b.JA_rlbend.v;

                        R_rlb.u=b.R_rlbend.u;
                        R_rlb.w=b.R_rlbend.w;
                        R_rlb.p=b.R_rlbend.p;
                        R_rlb.v=b.R_rlbend.v;

                        XI_r_rlb.u=b.XI_r_rlbend.u;
                        XI_r_rlb.w=b.XI_r_rlbend.w;
                        XI_r_rlb.p=b.XI_r_rlbend.p;
                        XI_r_rlb.v=b.XI_r_rlbend.v;

                        XI_z_rlb_C.u=b.XI_z_rlbend_C.u;
                        XI_z_rlb_C.w=b.XI_z_rlbend_C.w;
                        XI_z_rlb_C.p=b.XI_z_rlbend_C.p;
                        XI_z_rlb_C.v=b.XI_z_rlbend_C.v;

                        XI_z_rlb_E.u=b.XI_z_rlbend_E.u;
                        XI_z_rlb_E.w=b.XI_z_rlbend_E.w;
                        XI_z_rlb_E.p=b.XI_z_rlbend_E.p;
                        XI_z_rlb_E.v=b.XI_z_rlbend_E.v;

                end

                % inertia terms
                if creeping~=1
                    %  [d(J^2*r^ax)/dRs*rho0*xi_r*u0*w0]_delta.xi*delta.eta +[d(J^2*r^ax)/dRs*rho0*xi_r*w0^2]_delta.eta*delta.xi
                        F.p = (b.JA.p.^2.*(ax*R_rlb.p.*b.XI_r.p+XI_r_rlb.p.*b.R.p.^ax)+2*b.JA.p.*JA_rlb.p.*b.XI_r.p.*b.R.p.^ax).*rho(b.T.p).*b.w.p.^2;
                        FC.v = (b.JA.v.^2.*(ax*R_rlb.v.*b.XI_r.v+XI_r_rlb.v.*b.R.v.^ax)+2*b.JA.v.*JA_rlb.v.*b.XI_r.v.*b.R.v.^ax).*rho(b.T.v).*b.u.v.*b.w.v;

                else
                    F.p  = zeros(size(b.JA.p));
                    FC.v = zeros(size(b.JA.v));

                end

                % viscous terms with changed signs!
                if inviscid~=1
                     % -[mu*d(J*r^ax*xi_r*xi_z)/dRs*d(J*eta_z*u0)/dxi]_delta.xi*delta.eta
                        FC.v=FC.v - mu(b.T.v).*((JA_rlb.v.*b.R.v.^ax+ax*b.JA.v.*R_rlb.v).*b.XI_r.v.*b.XI_z.v+b.JA.v.*b.R.v.^ax.*(XI_r_rlb.v.*b.XI_z.v+b.XI_r.v.*XI_z_rlb_C.v)).*(b.JA.w(2:end,:).*b.ETA_z.w(2:end,:).*b.u.w(2:end,:)-b.JA.w(1:end-1,:).*b.ETA_z.w(1:end-1,:).*b.u.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];
                        FE.v=-mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.XI_r.v.*XI_z_rlb_E.v.*(b.JA.w(2:end,:).*b.ETA_z.w(2:end,:).*b.u.w(2:end,:)-b.JA.w(1:end-1,:).*b.ETA_z.w(1:end-1,:).*b.u.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % -[mu*J*r^ax*xi_r*xi_z*d(dJ/dRs*eta_z*u0)/dxi]_delta.xi*delta.eta
                        FC.v=FC.v - mu(b.T.v).*b.JA.v.*b.XI_r.v.*b.XI_z.v.*b.R.v.^ax.*(JA_rlb.w(2:end,:).*b.ETA_z.w(2:end,:).*b.u.w(2:end,:)-JA_rlb.w(1:end-1,:).*b.ETA_z.w(1:end-1,:).*b.u.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % [2/3*mu*d(xi_z)/dRs*d(J*r^ax*u0)/dxi]_delta.xi*delta.eta
                        FC.v=FC.v + 2/3*mu(b.T.v).*XI_z_rlb_C.v.*(b.JA.w(2:end,:).*b.R.w(2:end,:).^ax.*b.u.w(2:end,:)-b.JA.w(1:end-1,:).*b.R.w(1:end-1,:).^ax.*b.u.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];
                        FE.v=FE.v + 2/3*mu(b.T.v).*XI_z_rlb_E.v.*(b.JA.w(2:end,:).*b.R.w(2:end,:).^ax.*b.u.w(2:end,:)-b.JA.w(1:end-1,:).*b.R.w(1:end-1,:).^ax.*b.u.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % [2/3*mu*xi_z*d(d(J*r^ax)/dRs*u0)/dxi]_delta.xi*delta.eta
                        FC.v=FC.v + 2/3*mu(b.T.v).*b.XI_z.v.*((JA_rlb.w(2:end,:).*b.R.w(2:end,:).^ax+ax*b.JA.w(2:end,:).*R_rlb.w(2:end,:)).*b.u.w(2:end,:)-(JA_rlb.w(1:end-1,:).*b.R.w(1:end-1,:).^ax+ax*b.JA.w(1:end-1,:).*R_rlb.w(1:end-1,:)).*b.u.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % -[mu*d(J*r^ax*xi_r)/dRs*eta_z*d(J*eta_z*u0)/dxi]_delta.xi*delta.eta
                        FC.v=FC.v - mu(b.T.v).*((JA_rlb.v.*b.R.v.^ax+ax*b.JA.v.*R_rlb.v).*b.XI_r.v+b.JA.v.*b.R.v.^ax.*XI_r_rlb.v).*b.ETA_z.v.*(b.JA.u(:,2:end).*b.ETA_z.u(:,2:end).*b.u.u(:,2:end)-b.JA.u(:,1:end-1).*b.ETA_z.u(:,1:end-1).*b.u.u(:,1:end-1))./[b.DETA.zm(1,:); b.DETA.zm];

                    % -[mu*J*r^ax*xi_r*eta_z*d(dJ/dRs*eta_z*u0)/dxi]_delta.xi*delta.eta
                        F_w.v= - mu(b.T.v).*b.JA.v.*b.XI_r.v.*b.ETA_z.v.*b.R.v.^ax.*JA_rlb.u(:,1:end-1).*b.ETA_z.u(:,1:end-1).*b.u.u(:,1:end-1);
                        F_e.v= - mu(b.T.v).*b.JA.v.*b.XI_r.v.*b.ETA_z.v.*b.R.v.^ax.*JA_rlb.u(:,2:end).*b.ETA_z.u(:,2:end).*b.u.u(:,2:end);

                    % [2/3*mu*eta_z*d(J*r^ax*u0)/dxi]_delta.eta*delta.xi
                        F.p=F.p + 2/3*mu(b.T.p).*b.ETA_z.p.*((JA_rlb.u(2:end,2:end-1).*b.R.u(2:end,2:end-1).^ax+ax*b.JA.u(2:end,2:end-1).*R_rlb.u(2:end,2:end-1)).*b.u.u(2:end,2:end-1)-(JA_rlb.u(1:end-1,2:end-1).*b.R.u(1:end-1,2:end-1).^ax+ax*b.JA.u(1:end-1,2:end-1).*R_rlb.u(1:end-1,2:end-1)).*b.u.u(1:end-1,2:end-1))./b.DXI.c;

                    % -[mu*d(J*r^ax*(xi_r^2+2*xi_z^2)/dRs*d(J*xi_r*w0)/dxi]_delta.xi*delta.eta
                        FC.v=FC.v - mu(b.T.v).*((JA_rlb.v.*b.R.v.^ax+ax*b.JA.v.*R_rlb.v).*(b.XI_r.v.^2+2*b.XI_z.v.^2)+b.JA.v.*b.R.v.^ax.*(2*b.XI_r.v.*XI_r_rlb.v+4*b.XI_z.v.*XI_z_rlb_C.v)).*(b.JA.w(2:end,:).*b.XI_r.w(2:end,:).*b.w.w(2:end,:)-b.JA.w(1:end-1,:).*b.XI_r.w(1:end-1,:).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];
                        FE.v=FE.v - mu(b.T.v).*b.JA.v.*b.R.v.^ax.*4.*b.XI_z.v.*XI_z_rlb_E.v.*(b.JA.w(2:end,:).*b.XI_r.w(2:end,:).*b.w.w(2:end,:)-b.JA.w(1:end-1,:).*b.XI_r.w(1:end-1,:).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % -[mu*J*r^ax*(xi_r^2+2*xi_z^2*d(d(J*xi_r)/dRs*w0)/dxi]_delta.xi*delta.eta
                        FC.v=FC.v - mu(b.T.v).*b.JA.v.*b.R.v.^ax.*(b.XI_r.v.^2+2*b.XI_z.v.^2).*((JA_rlb.w(2:end,:).*b.XI_r.w(2:end,:)+b.JA.w(2:end,:).*XI_r_rlb.w(2:end,:)).*b.w.w(2:end,:)-(JA_rlb.w(1:end-1,:).*b.XI_r.w(1:end-1,:)+b.JA.w(1:end-1,:).*XI_r_rlb.w(1:end-1,:)).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % [mu*d(J*r^ax*xi_r*xi_z)/dRs*d(J*xi_z*w0)/dxi]_delta.xi*delta.eta
                        FC.v=FC.v + mu(b.T.v).*((JA_rlb.v.*b.R.v.^ax+ax*b.JA.v.*R_rlb.v).*b.XI_r.v.*b.XI_z.v+b.JA.v.*b.R.v.^ax.*(b.XI_z.v.*XI_r_rlb.v+b.XI_r.v.*XI_z_rlb_C.v)).*(b.JA.w(2:end,:).*b.XI_z.w(2:end,:).*b.w.w(2:end,:)-b.JA.w(1:end-1,:).*b.XI_z.w(1:end-1,:).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];
                        FE.v=FE.v + mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.XI_r.v.*XI_z_rlb_E.v.*(b.JA.w(2:end,:).*b.XI_z.w(2:end,:).*b.w.w(2:end,:)-b.JA.w(1:end-1,:).*b.XI_z.w(1:end-1,:).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % [mu*J*r^ax*xi_r*xi_z*d(d(J*xi_z)/dRs*w0)/dxi]_delta.xi*delta.eta
                        FC.v=FC.v + mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.XI_r.v.*b.XI_z.v.*((JA_rlb.w(2:end,:).*b.XI_z.w(2:end,:)+b.JA.w(2:end,:).*XI_z_rlb_C.w(2:end,:)).*b.w.w(2:end,:)-(JA_rlb.w(1:end-1,:).*b.XI_z.w(1:end-1,:)+b.JA.w(1:end-1,:).*XI_z_rlb_C.w(1:end-1,:)).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];
                        FE.v=FE.v + mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.XI_r.v.*b.XI_z.v.*(b.JA.w(2:end,:).*XI_z_rlb_E.w(2:end,:).*b.w.w(2:end,:)-b.JA.w(1:end-1,:).*XI_z_rlb_E.w(1:end-1,:).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % [2/3*mu*d(xi_z)/dRs*d(J*r^ax*w0)/deta]_delta.xi*delta.eta
                        FC.v=FC.v + 2/3*mu(b.T.v).*XI_z_rlb_C.v.*(b.JA.u(:,2:end).*b.R.u(:,2:end).^ax.*b.w.u(:,2:end)-b.JA.u(:,1:end-1).*b.R.u(:,1:end-1).^ax.*b.w.u(:,1:end-1))./[b.DETA.zm(1,:); b.DETA.zm];
                        FE.v=FE.v + 2/3*mu(b.T.v).*XI_z_rlb_E.v.*(b.JA.u(:,2:end).*b.R.u(:,2:end).^ax.*b.w.u(:,2:end)-b.JA.u(:,1:end-1).*b.R.u(:,1:end-1).^ax.*b.w.u(:,1:end-1))./[b.DETA.zm(1,:); b.DETA.zm];

                    % [2/3*mu*xi_z*d(d(J*r^ax)/dRs*w0)/deta]_delta.xi*delta.eta
                        F_w.v=F_w.v + 2/3*mu(b.T.v).*b.XI_z.v.*(JA_rlb.u(:,1:end-1).*b.R.u(:,1:end-1).^ax+ax*b.JA.u(:,1:end-1).*R_rlb.u(:,1:end-1).^ax).*b.w.u(:,1:end-1);
                        F_e.v=F_e.v + 2/3*mu(b.T.v).*b.XI_z.v.*(JA_rlb.u(:,2:end).*b.R.u(:,2:end).^ax+ax*b.JA.u(:,2:end).*R_rlb.u(:,2:end).^ax).*b.w.u(:,2:end);

                    % -[2*mu*d(J*r^ax*xi_z)/dRs*eta_z*d(J*xi_r*w0)/deta]_delta.xi*delta.eta
                        FC.v=FC.v - 2*mu(b.T.v).*((JA_rlb.v.*b.R.v.^ax+ax*b.JA.v.*R_rlb.v).*b.XI_z.v+b.JA.v.*b.R.v.^ax.*XI_z_rlb_C.v).*b.ETA_z.v.*(b.JA.u(:,2:end).*b.XI_r.u(:,2:end).*b.w.u(:,2:end)-b.JA.u(:,1:end-1).*b.XI_r.u(:,1:end-1).*b.w.u(:,1:end-1))./[b.DETA.zm(1,:); b.DETA.zm];
                        FE.v=FE.v - 2*mu(b.T.v).*b.JA.v.*b.R.v.^ax.*XI_z_rlb_E.v.*b.ETA_z.v.*(b.JA.u(:,2:end).*b.XI_r.u(:,2:end).*b.w.u(:,2:end)-b.JA.u(:,1:end-1).*b.XI_r.u(:,1:end-1).*b.w.u(:,1:end-1))./[b.DETA.zm(1,:); b.DETA.zm];

                    % -[2*mu*J*r^ax*xi_z*eta_z*d(d(J*xi_r)/dRs*w0)/deta]_delta.xi*delta.eta
                        F_w.v=F_w.v - 2*mu(b.T.v).*b.JA.v.*b.XI_z.v.*b.ETA_z.v.*b.R.v.^ax.*(JA_rlb.u(:,1:end-1).*b.XI_r.u(:,1:end-1)+b.JA.u(:,1:end-1).*XI_r_rlb.u(:,1:end-1)).*b.w.u(:,1:end-1);
                        F_e.v=F_e.v - 2*mu(b.T.v).*b.JA.v.*b.XI_z.v.*b.ETA_z.v.*b.R.v.^ax.*(JA_rlb.u(:,2:end).*b.XI_r.u(:,2:end)+b.JA.u(:,2:end).*XI_r_rlb.u(:,2:end)).*b.w.u(:,2:end);

                    % [mu*d(J*r^ax*xi_r)/dRs*eta_z*d(J*xi_z*w0)/deta]_delta.xi*delta.eta
                        FC.v=FC.v + mu(b.T.v).*((JA_rlb.v.*b.R.v.^ax+ax*b.JA.v.*R_rlb.v).*b.XI_r.v+b.JA.v.*b.R.v.^ax.*XI_r_rlb.v).*b.ETA_z.v.*(b.JA.u(:,2:end).*b.XI_z.u(:,2:end).*b.w.u(:,2:end)-b.JA.u(:,1:end-1).*b.XI_z.u(:,1:end-1).*b.w.u(:,1:end-1))./[b.DETA.zm(1,:); b.DETA.zm];

                    % [mu*J*r^ax*xi_r*eta_z*d(d(J*xi_z)/dRs*w0)/deta]_delta.xi*delta.eta
                        F_w.v=F_w.v + mu(b.T.v).*b.JA.v.*b.XI_r.v.*b.ETA_z.v.*b.R.v.^ax.*(JA_rlb.u(:,1:end-1).*b.XI_z.u(:,1:end-1)+b.JA.u(:,1:end-1).*XI_z_rlb_C.u(:,1:end-1)).*b.w.u(:,1:end-1);
                        F_e.v=F_e.v + mu(b.T.v).*b.JA.v.*b.XI_r.v.*b.ETA_z.v.*b.R.v.^ax.*(JA_rlb.u(:,2:end).*b.XI_z.u(:,2:end)+b.JA.u(:,2:end).*XI_z_rlb_C.u(:,2:end)).*b.w.u(:,2:end);
                        F_W.v=mu(b.T.v).*b.JA.v.*b.XI_r.v.*b.ETA_z.v.*b.R.v.^ax.*b.JA.u(:,1:end-1).*XI_z_rlb_E.u(:,1:end-1).*b.w.u(:,1:end-1);
                        F_E.v=mu(b.T.v).*b.JA.v.*b.XI_r.v.*b.ETA_z.v.*b.R.v.^ax.*b.JA.u(:,2:end).*XI_z_rlb_E.u(:,2:end).*b.w.u(:,2:end);

                    % -[2*mu*d(J*r^ax*xi_z)/dRs*eta_z*d(J*xi_r*w0)/dxi]_delta.eta*delta.xi
                        F.p=F.p - 2*mu(b.T.p).*((JA_rlb.p.*b.R.p.^ax+ax*b.JA.p.*R_rlb.p).*b.XI_z.p+b.JA.p.*b.R.p.^ax.*XI_z_rlb_C.p).*b.ETA_z.p.*(b.JA.u(2:end,2:end-1).*b.XI_r.u(2:end,2:end-1).*b.w.u(2:end,2:end-1)-b.JA.u(1:end-1,2:end-1).*b.XI_r.u(1:end-1,2:end-1).*b.w.u(1:end-1,2:end-1))./b.DXI.c;
                        FE.p= - 2*mu(b.T.p).*b.JA.p.*b.R.p.^ax.*XI_z_rlb_E.p.*b.ETA_z.p.*(b.JA.u(2:end,2:end-1).*b.XI_r.u(2:end,2:end-1).*b.w.u(2:end,2:end-1)-b.JA.u(1:end-1,2:end-1).*b.XI_r.u(1:end-1,2:end-1).*b.w.u(1:end-1,2:end-1));

                    % -[2*mu*J*r^ax*xi_z*eta_z*d(d(J*xi_r)/dRs*w0)/dxi]_delta.eta*delta.xi
                        F.p=F.p - 2*mu(b.T.p).*b.JA.p.*b.R.p.^ax.*b.XI_z.p.*b.ETA_z.p.*((JA_rlb.u(2:end,2:end-1).*b.XI_r.u(2:end,2:end-1) + b.JA.u(2:end,2:end-1).*XI_r_rlb.u(2:end,2:end-1)).*b.w.u(2:end,2:end-1) - (JA_rlb.u(1:end-1,2:end-1).*b.XI_r.u(1:end-1,2:end-1) + b.JA.u(1:end-1,2:end-1).*XI_r_rlb.u(1:end-1,2:end-1)).*b.w.u(1:end-1,2:end-1))./b.DXI.c;

                    % -[2*mu*d(J*r^ax)/dRs*eta_z^2*d(J*xi_r*w0)/deta]_delta.eta*delta.xi
                        F.p=F.p - 2*mu(b.T.p).*(JA_rlb.p.*b.R.p.^ax+ax*b.JA.p.*R_rlb.p).*b.ETA_z.p.^2.*(b.JA.w(2:end-1,2:end).*b.XI_r.w(2:end-1,2:end).*b.w.w(2:end-1,2:end)-b.JA.w(2:end-1,1:end-1).*b.XI_r.w(2:end-1,1:end-1).*b.w.w(2:end-1,1:end-1))./b.DETA.c;

                    % -[2*mu*J*r^ax*eta_z^2*d(d(J*xi_r)/dRs*w0)/deta]_delta.eta*delta.xi
                        F_W.p= -2*mu(b.T.p).*b.JA.p.*b.R.p.^ax.*b.ETA_z.p.^2.*(JA_rlb.w(2:end-1,1:end-1).*b.XI_r.w(2:end-1,1:end-1) + b.JA.w(2:end-1,1:end-1).*XI_r_rlb.w(2:end-1,1:end-1)).*b.w.w(2:end-1,1:end-1)./b.DETA.c;
                        F_E.p= -2*mu(b.T.p).*b.JA.p.*b.R.p.^ax.*b.ETA_z.p.^2.*(JA_rlb.w(2:end-1,2:end).*b.XI_r.w(2:end-1,2:end) + b.JA.w(2:end-1,2:end).*XI_r_rlb.w(2:end-1,2:end)).*b.w.w(2:end-1,2:end)./b.DETA.c;

                    % [2/3*mu*eta_z*d(d(J*r^ax)/dRs*w0)/deta]_delta.eta*delta.xi
                        F_W.p=F_W.p + 2/3*mu(b.T.p).*b.ETA_z.p.*(JA_rlb.w(2:end-1,1:end-1).*b.R.w(2:end-1,1:end-1).^ax + ax*b.JA.w(2:end-1,1:end-1).*R_rlb.w(2:end-1,1:end-1)).*b.w.w(2:end-1,1:end-1)./b.DETA.c;
                        F_E.p=F_E.p + 2/3*mu(b.T.p).*b.ETA_z.p.*(JA_rlb.w(2:end-1,2:end).*b.R.w(2:end-1,2:end).^ax + ax*b.JA.w(2:end-1,2:end).*R_rlb.w(2:end-1,2:end)).*b.w.w(2:end-1,2:end)./b.DETA.c;

                else
                    F_w.v = zeros(size(b.JA.v));
                    F_e.v = F_w.v;
                    FE.v  = F_w.v;
                    F_W.v = F_w.v;
                    F_E.v = F_w.v;
                    FE.p  = zeros(size(b.JA.p));
                    F_W.p = FE.p;
                    F_E.p = FE.p;

                end

                % pressure terms with changed signs!
                    % [d(J*xi_z*r^ax)/dRs*p0]_delta.xi*delta.eta + [d(J*eta_z*r^ax)/dRs*p0]_delta.eta*delta.xi
                        F.p=F.p + (JA_rlb.p.*b.R.p.^ax+ax*b.JA.p.*R_rlb.p).*b.ETA_z.p.*b.p.p;
                        FC.v=FC.v + (JA_rlb.v.*b.XI_z.v.*b.R.v.^ax+b.JA.v.*(XI_z_rlb_C.v.*b.R.v.^ax+ax*b.XI_z.v.*R_rlb.v)).*b.p.v;
                        FE.v=FE.v + b.JA.v.*XI_z_rlb_E.v.*b.R.v.^ax.*b.p.v;

                % gravity term with changed signs!
                    % g*rho*d(J*r)*delta.xi*delta.eta
                        if isfield(flowopt,'boussinesq')==1 && flowopt.boussinesq==1
                            rho = @(theta) rho(b.T_0).*(1-b.beta*(b.T.w-b.T_0));
                        end
                        FC.w= - g*rho(b.T.w).*(JA_rlb.w.*b.R.w.^ax+ax*b.JA.w.*R_rlb.w);
                        rho = str2func(b.rho);

                Rs_w=zeros(b.K+2,b.J+1); Rs_e=Rs_w; Rs_W=Rs_w; Rs_E=Rs_w; Rs_ad_w=Rs_w; Rs_ad_e=Rs_w; Rs_ad_W=Rs_w; Rs_ad_E=Rs_w;

                Rs_W(2:end-1,2:end-1)    = -b.W.w(2:end-1,1:end-1).*(-F_W.v(2:end,2:end-1)+F_W.v(1:end-1,2:end-1) - FE.p(:,1:end-1) - F_W.p(:,1:end-1).*b.DXI.zm(:,2:end-1));
                Rs_ad_w(2:end-1,2:end-1) = b.W.w(2:end-1,2:end).*(FC.v(2:end,2:end-1)-FC.v(1:end-1,2:end-1)+FC.w(2:end-1,2:end-1).*b.DXI.zm(:,2:end-1)).*b.DETA.zm(:,2:end-1) + b.W.w(2:end-1,2:end).*(-F_W.v(2:end,2:end-1) + F_W.v(1:end-1,2:end-1) - F_E.v(2:end,2:end-1) + F_E.v(1:end-1,2:end-1) - FE.p(:,1:end-1) - FE.p(:,2:end) - (F_W.p(:,2:end) + F_E.p(:,1:end-1)).*b.DXI.zm(:,2:end-1));
                Rs_ad_e(2:end-1,2:end-1) = b.E.w(2:end-1,1:end-1).*(FC.v(2:end,2:end-1)-FC.v(1:end-1,2:end-1)+FC.w(2:end-1,2:end-1).*b.DXI.zm(:,2:end-1)).*b.DETA.zm(:,2:end-1) - b.E.w(2:end-1,1:end-1).*(F_E.v(2:end,2:end-1)-F_E.v(1:end-1,2:end-1) + F_W.v(2:end,2:end-1) - F_W.v(1:end-1,2:end-1) + FE.p(:,1:end-1) + FE.p(:,2:end) + (F_W.p(:,2:end) + F_E.p(:,1:end-1)).*b.DXI.zm(:,2:end-1));
                Rs_ad_W(2:end-1,2:end-1) = -[zeros(b.K,1) b.E.w(2:end-1,1:end-2)].*(-F_W.v(2:end,2:end-1)+F_W.v(1:end-1,2:end-1) - FE.p(:,1:end-1) - F_W.p(:,1:end-1).*b.DXI.zm(:,2:end-1));
                Rs_ad_E(2:end-1,2:end-1) = [b.W.w(2:end-1,3:end) zeros(b.K,1)].*(F_E.v(2:end,2:end-1)-F_E.v(1:end-1,2:end-1) + FE.p(:,2:end) + F_E.p(:,2:end).*b.DXI.zm(:,2:end-1));
                Rs_w(2:end-1,2:end-1)    = Rs_ad_w(2:end-1,2:end-1) + Rs_ad_W(2:end-1,2:end-1) - (FE.v(2:end,2:end-1)-FE.v(1:end-1,2:end-1)).*b.DETA.zm(:,2:end-1) + -F_w.v(2:end,2:end-1)+F_w.v(1:end-1,2:end-1) - F.p(:,1:end-1).*b.DXI.c(:,1:end-1);
                Rs_e(2:end-1,2:end-1)    = Rs_ad_e(2:end-1,2:end-1) + Rs_ad_E(2:end-1,2:end-1) + (FE.v(2:end,2:end-1)-FE.v(1:end-1,2:end-1)).*b.DETA.zm(:,2:end-1) + F_e.v(2:end,2:end-1)-F_e.v(1:end-1,2:end-1) + F.p(:,2:end).*b.DXI.c(:,2:end);
                Rs_E(2:end-1,2:end-1)    = b.E.w(2:end-1,2:end).*(F_E.v(2:end,2:end-1)-F_E.v(1:end-1,2:end-1) + FE.p(:,2:end) + F_E.p(:,2:end).*b.DXI.zm(:,2:end-1));

                %Rs_W=Rs_W(:);
                %Rs_w=Rs_w(:);
                %Rs_e=Rs_e(:);
                %Rs_E=Rs_E(:);

                index=kron(ones(1,b.J+1),1:b.J+1:(b.K+2)*(b.J+1))+ kron(0:(b.K+2)*(b.J+1)+1:(b.K+2)*(b.J+1)*(b.J+1),ones(1,(b.K+2)));
                Rs1=sparse((b.K+2)*(b.J+1),b.J+1); Rs2=Rs1; Rs3=Rs1; Rs4=Rs1;
                Rs1(index)=Rs_w;
                Rs2(index)=Rs_e;
                Rs3(index)=Rs_W;
                Rs4(index)=Rs_E;

                Ja.z_m.Rs_all=[sparse((b.K+2)*(b.J+1),1) Rs1 sparse((b.K+2)*(b.J+1),2)]+[sparse((b.K+2)*(b.J+1),2) Rs2 sparse((b.K+2)*(b.J+1),1)]+[Rs3 sparse((b.K+2)*(b.J+1),3)]+[sparse((b.K+2)*(b.J+1),3) Rs4];
                Ja.z_m.Rs_all=Ja.z_m.Rs_all(:,2:end-1);
                Ja.z_m.Rs=[];
                Ja_counter.z_m.Rs=[];

                for n=1:length(b.bc.z)/2
                    if (o==1 && strcmp(b.bc.z{n}(1),'s')==1 && strcmp(b.bc.z{n}(2),'s')==1) || (o==2 && strcmp(b.bc.z{index_bc_z+n}(1),'s')==1 && strcmp(b.bc.z{index_bc_z+n}(2),'s')==1)
                        Ja.z_m.Rs_part=Ja.z_m.Rs_all(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                        if (n>1 && strcmp(b.bc.z{(o-1)*index_bc_z+n-1}(1),'s')~=1 && strcmp(b.bc.z{(o-1)*index_bc_z+n-1}(2),'s')~=1)
                            Ja.z_m.Rs_part(find(b.z.w==b.geom.z(n))+b.J:b.J+1:end-b.J-2,1)=0;
                            Ja.z_m.Rs_part(find(b.z.w==b.geom.z(n))+b.J+1:b.J+1:end-b.J-1,1)=Ja.z_m.Rs_part(find(b.z.w==b.geom.z(n))+b.J+1:b.J+1:end-b.J-1,1)-Rs_ad_e(2:end-1,b.z.w==b.geom.z(n));
                            Ja.z_m.Rs_part(find(b.z.w==b.geom.z(n))+b.J+2:b.J+1:end-b.J,1)=Ja.z_m.Rs_part(find(b.z.w==b.geom.z(n))+b.J+2:b.J+1:end-b.J,1)-Rs_ad_W(2:end-1,find(b.z.w==b.geom.z(n))+1);

                        end
                        if (n<length(b.bc.z)/2 && strcmp(b.bc.z{(o-1)*index_bc_z+n+1}(1),'s')~=1 && strcmp(b.bc.z{(o-1)*index_bc_z+n+1}(2),'s')~=1)
                            Ja.z_m.Rs_part(find(b.z.w==b.geom.z(n+1))+b.J+2:b.J+1:end-b.J,end)=0;
                            Ja.z_m.Rs_part(find(b.z.w==b.geom.z(n+1))+b.J+1:b.J+1:end-b.J-1,end)=Ja.z_m.Rs_part(find(b.z.w==b.geom.z(n+1))+b.J+1:b.J+1:end-b.J-1,end)-Rs_ad_w(2:end-1,b.z.w==b.geom.z(n+1));
                            Ja.z_m.Rs_part(find(b.z.w==b.geom.z(n+1))+b.J:b.J+1:end-b.J-2,end)=Ja.z_m.Rs_part(find(b.z.w==b.geom.z(n+1))+b.J:b.J+1:end-b.J-2,end)-Rs_ad_E(2:end-1,find(b.z.w==b.geom.z(n+1))-1);

                        end

                        if eval(['b.' b.bc.z{(o-1)*index_bc_z+n} '.o==1'])
                            Ja.z_m.Rs=[Ja.z_m.Rs Ja.z_m.Rs_part];

                        elseif eval(['b.' b.bc.z{(o-1)*index_bc_z+n} '.o==0'])
                            Ja_counter.z_m.Rs=[Ja_counter.z_m.Rs Ja.z_m.Rs_part];

                        end

                    end

                end

                if o==1
                    b.Ja.z_m.Rs1=Ja.z_m.Rs;
                    b.Ja_counter.z_m.Rs1=Ja_counter.z_m.Rs;

                else
                    b.Ja.z_m.Rsend=Ja.z_m.Rs;
                    b.Ja_counter.z_m.Rsend=Ja_counter.z_m.Rs;

                end

            end
            
        end

    end