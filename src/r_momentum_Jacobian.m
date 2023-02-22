function[b]=r_momentum_Jacobian(b, flowopt)
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
    
    b.Ja.r_m.u=sparse((b.K+1)*(b.J+2),(b.K+1)*(b.J+2));
    b.Ja.r_m.w=sparse((b.K+1)*(b.J+2),(b.K+2)*(b.J+1));
    
    if creeping~=1
    % [J^2*r^ax*rho0*(eta_z*2*u0-xi_z*w0)]_delta.xi*delta.eta + [J^2*r^ax*rho0*eta_z*w0]_delta.eta*delta.xi
        F.p=b.JA.p.^2.*b.R.p.^ax.*rho(b.T.p).*(b.ETA_z.p.*b.u.p*2-b.XI_z.p.*b.w.p);
        F.v=b.JA.v.^2.*b.R.v.^ax.*rho(b.T.v).*b.ETA_z.v.*b.w.v;

        u_S=zeros(b.K+1,b.J+2); u_N=u_S; u_C=u_S; u_W=u_S; u_E=u_S;

        u_S(1:end-2,2:end-1) = -F.p(1:end-1,:)/2.*b.DETA.rm(2:end-1,:);
        u_W(2:end-1,1)       = -F.v(2:end-1,1).*b.DXI.rm(2:end-1,1);
        u_W(2:end-1,2:end-2) = -F.v(2:end-1,2:end-1)/2.*b.DXI.rm(2:end-1,2:end);
        u_E(2:end-1,end)     = F.v(2:end-1,end).*b.DXI.rm(2:end-1,end);
        u_E(2:end-1,3:end-1) = F.v(2:end-1,2:end-1)/2.*b.DXI.rm(2:end-1,1:end-1);
        u_N(3:end,2:end-1)   = F.p(2:end,:)/2.*b.DETA.rm(2:end-1,:);
        u_C(2:end-1,2:end-1) = u_C(2:end-1,2:end-1)+u_S(1:end-2,2:end-1);
        u_C(2:end-1,3:end-1) = u_C(2:end-1,3:end-1)+u_W(2:end-1,2:end-2);
        u_C(2:end-1,2:end-2) = u_C(2:end-1,2:end-2)+u_E(2:end-1,3:end-1);
        u_C(2:end-1,2:end-1) = u_C(2:end-1,2:end-1)+u_N(3:end,2:end-1);

        u_S=u_S.'; u_S=u_S(:);
        u_W=u_W.'; u_W=u_W(:);
        u_E=u_E.'; u_E=u_E(:);
        u_N=u_N.'; u_N=u_N(:);
        u_C=u_C.'; u_C=u_C(:);

        b.Ja.r_m.u=spdiags([u_S u_W u_C u_E u_N],[-b.J-2 -1 0 1 b.J+2],(b.K+1)*(b.J+2),(b.K+1)*(b.J+2));

    % [J^2*r^ax*rho0*(eta_z*u0-xi_z*2*w0)]_delta.eta*delta.xi + [J^2*r^ax*rho0*(-xi_z*u0)]_delta.xi*delta.eta
        F.v=b.JA.v.^2.*b.R.v.^ax.*rho(b.T.v).*(b.ETA_z.v.*b.u.v-b.XI_z.v.*b.w.v*2);
        F.p=b.JA.p.^2.*b.R.p.^ax.*rho(b.T.p).*(-b.XI_z.p.*b.u.p);

        w_sw=zeros(b.K+1,b.J+2); w_se=w_sw; w_nw=w_sw; w_ne=w_sw;

        w_sw(2:end-1,2:end-1) = -F.v(2:end-1,1:end-1)/2.*b.DXI.rm(2:end-1,:)-F.p(1:end-1,:)/2.*b.DETA.rm(2:end-1,:);
        w_se(2:end-1,3:end)   = F.v(2:end-1,2:end)/2.*b.DXI.rm(2:end-1,:)-F.p(1:end-1,:)/2.*b.DETA.rm(2:end-1,:);
        w_nw(3:end,2:end-1)   = -F.v(2:end-1,1:end-1)/2.*b.DXI.rm(2:end-1,:)+F.p(2:end,:)/2.*b.DETA.rm(2:end-1,:);
        w_ne(3:end,3:end)     = F.v(2:end-1,2:end)/2.*b.DXI.rm(2:end-1,:)+F.p(2:end,:)/2.*b.DETA.rm(2:end-1,:);

        w_sw=w_sw.'; w_sw=w_sw(:);
        w_se=w_se.'; w_se=w_se(:);
        w_nw=w_nw.'; w_nw=w_nw(:);
        w_ne=w_ne.'; w_ne=w_ne(:);

        b.Ja.r_m.w=spdiags([w_sw w_se w_nw w_ne],[0 1 b.J+2 b.J+3],(b.K+1)*(b.J+2),(b.K+1)*(b.J+2));
        b.Ja.r_m.w(:,1:b.J+2:end)=[];
        b.Ja.r_m.w=[b.Ja.r_m.w sparse((b.K+1)*(b.J+2),b.J+1)];
        
    end
    
    if inviscid~=1
        b.Ja.r_m.u=b.Ja.r_m.u + b.r_m_visc.u;
        b.Ja.r_m.w=b.Ja.r_m.w + b.r_m_visc.w;
        
    end

    if energy==1
        F.p=zeros(size(b.T.p)); F.v=zeros(size(b.T.v)); F1.p=F.p; F2.p=F.p; F3.p=F.p; F4.p=F.p; F5.p=F.p; F6.p=F.p; F7.p=F.p; F1.v=F.v; F2.v=F.v; F3.v=F.v; F4.v=F.v; F5.v=F.v; F1.u=zeros(size(b.T.u)); F2.u=F1.u; F3.u=F1.u; F4.u=F1.u;
        
        if creeping~=1
        % [J^2*r^ax*(eta_z*u0^2-xi_z*u0*w0)*drho0/dT0]_delta.xi*delta.eta + [J^2*r^ax*(eta_z*u0*w0-xi_z*w0^2)*drho0/dT0]_delta.eta*delta.xi
            F.p=b.JA.p.^2.*b.R.p.^ax.*drho(b.T.p).*(b.ETA_z.p.*b.u.p.^2-b.XI_z.p.*b.u.p.*b.w.p);
            F.v=b.JA.v.^2.*b.R.v.^ax.*drho(b.T.v).*(b.ETA_z.v.*b.u.v.*b.w.v-b.XI_z.v.*b.w.v.^2);

        end
        
        if inviscid~=1
        % viscous terms with changed signs!
            F1.p=2/3*dmu(b.T.p).*b.XI_r.p;
            F2.p=-dmu(b.T.p).*b.JA.p.*(2*b.XI_r.p.^2+b.XI_z.p.^2).*b.R.p.^ax;
            F3.p=-dmu(b.T.p).*b.JA.p.*b.XI_z.p.*b.ETA_z.p.*b.R.p.^ax;

            F4.p=2/3*dmu(b.T.p).*b.XI_r.p;
            F5.p=dmu(b.T.p).*b.JA.p.*(2*b.XI_r.p.^2+b.XI_z.p.^2).*b.R.p.^ax;
            F6.p=dmu(b.T.p).*b.JA.p.*b.XI_z.p.*b.ETA_z.p.*b.R.p.^ax;
            F7.p=-dmu(b.T.p).*b.JA.p.*b.XI_r.p.*b.XI_z.p.*b.R.p.^ax;

            F1.v=-dmu(b.T.v).*b.JA.v.*b.XI_z.v.*b.ETA_z.v.*b.R.v.^ax;
            F2.v=-dmu(b.T.v).*b.JA.v.*b.ETA_z.v.^2.*b.R.v.^ax;

            F3.v=dmu(b.T.v).*b.JA.v.*b.XI_z.v.*b.ETA_z.v.*b.R.v.^ax;
            F4.v=dmu(b.T.v).*b.JA.v.*b.ETA_z.v.^2.*b.R.v.^ax;
            F5.v=-dmu(b.T.v).*b.JA.v.*b.XI_r.v.*b.ETA_z.v.*b.R.v.^ax;

            F1.u=ax*2*dmu(b.T.u).*b.JA.u.^2.*b.ETA_z.u./b.R.u.^ax;
            F2.u=-ax*2/3*dmu(b.T.u)./b.R.u.^ax;
            
            F3.u=-ax*2*dmu(b.T.u).*b.JA.u.^2.*b.XI_z.u./b.R.u.^ax;
            F4.u=-ax*2/3*dmu(b.T.u)./b.R.u.^ax;
            
        end

            F_s.p=F2.p(1:end-1,:).*(b.JA.u(2:end-1,2:end-1).*b.ETA_z.u(2:end-1,2:end-1).*b.u.u(2:end-1,2:end-1) - b.JA.u(1:end-2,2:end-1).*b.ETA_z.u(1:end-2,2:end-1).*b.u.u(1:end-2,2:end-1))./b.DXI.c(1:end-1,:);
            F_s.p=F_s.p + F3.p(1:end-1,:).*(b.JA.w(2:end-2,2:end).*b.ETA_z.w(2:end-2,2:end).*b.u.w(2:end-2,2:end) - b.JA.w(2:end-2,1:end-1).*b.ETA_z.w(2:end-2,1:end-1).*b.u.w(2:end-2,1:end-1))./b.DETA.rm(2:end-1,:);
            F_s.p=F_s.p + F1.p(1:end-1,:).*(b.JA.u(2:end-1,2:end-1).*b.R.u(2:end-1,2:end-1).^ax.*b.u.u(2:end-1,2:end-1) - b.JA.u(1:end-2,2:end-1).*b.R.u(1:end-2,2:end-1).^ax.*b.u.u(1:end-2,2:end-1))./b.DXI.c(1:end-1,:);
            F_s.p=F_s.p + F5.p(1:end-1,:).*(b.JA.u(2:end-1,2:end-1).*b.XI_z.u(2:end-1,2:end-1).*b.w.u(2:end-1,2:end-1) - b.JA.u(1:end-2,2:end-1).*b.XI_z.u(1:end-2,2:end-1).*b.w.u(1:end-2,2:end-1))./b.DXI.c(1:end-1,:);
            F_s.p=F_s.p + F7.p(1:end-1,:).*(b.JA.u(2:end-1,2:end-1).*b.XI_r.u(2:end-1,2:end-1).*b.w.u(2:end-1,2:end-1) - b.JA.u(1:end-2,2:end-1).*b.XI_r.u(1:end-2,2:end-1).*b.w.u(1:end-2,2:end-1))./b.DXI.c(1:end-1,:);
            F_s.p=F_s.p + F6.p(1:end-1,:).*(b.JA.w(2:end-2,2:end).*b.XI_z.w(2:end-2,2:end).*b.w.w(2:end-2,2:end) - b.JA.w(2:end-2,1:end-1).*b.XI_z.w(2:end-2,1:end-1).*b.w.w(2:end-2,1:end-1))./b.DETA.rm(2:end-1,:);
            F_s.p=F_s.p + F4.p(1:end-1,:).*(b.JA.w(2:end-2,2:end).*b.R.w(2:end-2,2:end).^ax.*b.w.w(2:end-2,2:end) - b.JA.w(2:end-2,1:end-1).*b.R.w(2:end-2,1:end-1).^ax.*b.w.w(2:end-2,1:end-1))./b.DETA.rm(2:end-1,:);
            
            F_n.p=F2.p(2:end,:).*(b.JA.u(3:end,2:end-1).*b.ETA_z.u(3:end,2:end-1).*b.u.u(3:end,2:end-1) - b.JA.u(2:end-1,2:end-1).*b.ETA_z.u(2:end-1,2:end-1).*b.u.u(2:end-1,2:end-1))./b.DXI.c(2:end,:);
            F_n.p=F_n.p + F3.p(2:end,:).*(b.JA.w(3:end-1,2:end).*b.ETA_z.w(3:end-1,2:end).*b.u.w(3:end-1,2:end) - b.JA.w(3:end-1,1:end-1).*b.ETA_z.w(3:end-1,1:end-1).*b.u.w(3:end-1,1:end-1))./b.DETA.rm(2:end-1,:);
            F_n.p=F_n.p + F1.p(2:end,:).*(b.JA.u(3:end,2:end-1).*b.R.u(3:end,2:end-1).^ax.*b.u.u(3:end,2:end-1) - b.JA.u(2:end-1,2:end-1).*b.R.u(2:end-1,2:end-1).^ax.*b.u.u(2:end-1,2:end-1))./b.DXI.c(2:end,:);
            F_n.p=F_n.p + F5.p(2:end,:).*(b.JA.u(3:end,2:end-1).*b.XI_z.u(3:end,2:end-1).*b.w.u(3:end,2:end-1) - b.JA.u(2:end-1,2:end-1).*b.XI_z.u(2:end-1,2:end-1).*b.w.u(2:end-1,2:end-1))./b.DXI.c(2:end,:);
            F_n.p=F_n.p + F7.p(2:end,:).*(b.JA.u(3:end,2:end-1).*b.XI_r.u(3:end,2:end-1).*b.w.u(3:end,2:end-1) - b.JA.u(2:end-1,2:end-1).*b.XI_r.u(2:end-1,2:end-1).*b.w.u(2:end-1,2:end-1))./b.DXI.c(2:end,:);
            F_n.p=F_n.p + F6.p(2:end,:).*(b.JA.w(3:end-1,2:end).*b.XI_z.w(3:end-1,2:end).*b.w.w(3:end-1,2:end) - b.JA.w(3:end-1,1:end-1).*b.XI_z.w(3:end-1,1:end-1).*b.w.w(3:end-1,1:end-1))./b.DETA.rm(2:end-1,:);
            F_n.p=F_n.p + F4.p(2:end,:).*(b.JA.w(3:end-1,2:end).*b.R.w(3:end-1,2:end).^ax.*b.w.w(3:end-1,2:end) - b.JA.w(3:end-1,1:end-1).*b.R.w(3:end-1,1:end-1).^ax.*b.w.w(3:end-1,1:end-1))./b.DETA.rm(2:end-1,:);
            
            F_w.v=F1.v(2:end-1,1:end-1).*(b.JA.w(3:end-1,1:end-1).*b.ETA_z.w(3:end-1,1:end-1).*b.u.w(3:end-1,1:end-1) - b.JA.w(2:end-2,1:end-1).*b.ETA_z.w(2:end-2,1:end-1).*b.u.w(2:end-2,1:end-1))./b.DXI.rm(2:end-1,:);
            F_w.v=F_w.v + F2.v(2:end-1,1:end-1).*(b.JA.u(2:end-1,2:end-1).*b.ETA_z.u(2:end-1,2:end-1).*b.u.u(2:end-1,2:end-1) - b.JA.u(2:end-1,1:end-2).*b.ETA_z.u(2:end-1,1:end-2).*b.u.u(2:end-1,1:end-2))./b.DETA.zm(1:end-1,1:end-1);
            F_w.v=F_w.v + F5.v(2:end-1,1:end-1).*(b.JA.w(3:end-1,1:end-1).*b.XI_r.w(3:end-1,1:end-1).*b.w.w(3:end-1,1:end-1) - b.JA.w(2:end-2,1:end-1).*b.XI_r.w(2:end-2,1:end-1).*b.w.w(2:end-2,1:end-1))./b.DXI.rm(2:end-1,:);
            F_w.v=F_w.v + F3.v(2:end-1,1:end-1).*(b.JA.w(3:end-1,1:end-1).*b.XI_z.w(3:end-1,1:end-1).*b.w.w(3:end-1,1:end-1) - b.JA.w(2:end-2,1:end-1).*b.XI_z.w(2:end-2,1:end-1).*b.w.w(2:end-2,1:end-1))./b.DXI.rm(2:end-1,:);
            F_w.v=F_w.v + F4.v(2:end-1,1:end-1).*(b.JA.u(2:end-1,2:end-1).*b.XI_z.u(2:end-1,2:end-1).*b.w.u(2:end-1,2:end-1) - b.JA.u(2:end-1,1:end-2).*b.XI_z.u(2:end-1,1:end-2).*b.w.u(2:end-1,1:end-2))./b.DETA.zm(1:end-1,1:end-1);
            
            F_e.v=F1.v(2:end-1,2:end).*(b.JA.w(3:end-1,2:end).*b.ETA_z.w(3:end-1,2:end).*b.u.w(3:end-1,2:end) - b.JA.w(2:end-2,2:end).*b.ETA_z.w(2:end-2,2:end).*b.u.w(2:end-2,2:end))./b.DXI.rm(2:end-1,:);
            F_e.v=F_e.v + F2.v(2:end-1,2:end).*(b.JA.u(2:end-1,3:end).*b.ETA_z.u(2:end-1,3:end).*b.u.u(2:end-1,3:end) - b.JA.u(2:end-1,2:end-1).*b.ETA_z.u(2:end-1,2:end-1).*b.u.u(2:end-1,2:end-1))./b.DETA.zm(1:end-1,2:end);
            F_e.v=F_e.v + F5.v(2:end-1,2:end).*(b.JA.w(3:end-1,2:end).*b.XI_r.w(3:end-1,2:end).*b.w.w(3:end-1,2:end) - b.JA.w(2:end-2,2:end).*b.XI_r.w(2:end-2,2:end).*b.w.w(2:end-2,2:end))./b.DXI.rm(2:end-1,:);
            F_e.v=F_e.v + F3.v(2:end-1,2:end).*(b.JA.w(3:end-1,2:end).*b.XI_z.w(3:end-1,2:end).*b.w.w(3:end-1,2:end) - b.JA.w(2:end-2,2:end).*b.XI_z.w(2:end-2,2:end).*b.w.w(2:end-2,2:end))./b.DXI.rm(2:end-1,:);
            F_e.v=F_e.v + F4.v(2:end-1,2:end).*(b.JA.u(2:end-1,3:end).*b.XI_z.u(2:end-1,3:end).*b.w.u(2:end-1,3:end) - b.JA.u(2:end-1,2:end-1).*b.XI_z.u(2:end-1,2:end-1).*b.w.u(2:end-1,2:end-1))./b.DETA.zm(1:end-1,2:end);
            
            F_C.u=F1.u(2:end-1,2:end-1).*b.u.u(2:end-1,2:end-1).*b.DXI.rm(2:end-1,:).*b.DETA.rm(2:end-1,:);
            F_C.u=F_C.u + F2.u(2:end-1,2:end-1).*(b.JA.T(3:end-1,2:end-1).*b.R.T(3:end-1,2:end-1).^ax.*b.u.T(3:end-1,2:end-1) - b.JA.T(2:end-2,2:end-1).*b.R.T(2:end-2,2:end-1).^ax.*b.u.T(2:end-2,2:end-1)).*b.DETA.rm(2:end-1,:);
            F_C.u=F_C.u + F3.u(2:end-1,2:end-1).*b.w.u(2:end-1,2:end-1).*b.DXI.rm(2:end-1,:).*b.DETA.rm(2:end-1,:);
            F_C.u=F_C.u + F4.u(2:end-1,2:end-1).*(b.JA.v(2:end-1,2:end).*b.R.v(2:end-1,2:end).^ax.*b.w.v(2:end-1,2:end) - b.JA.v(2:end-1,1:end-1).*b.R.v(2:end-1,1:end-1).^ax.*b.w.v(2:end-1,1:end-1)).*b.DXI.rm(2:end-1,:);
            
            T_Wsw=zeros(b.K+1,b.J+2); T_s=T_Wsw; T_Ese=T_Wsw; T_Wnw=T_Wsw; T_n=T_Wsw; T_Ene=T_Wsw;

            %%{
            T_Wsw(2:end-1,1:end-2)   = -b.SW.v(2:end,:).*(F.v(2:end-1,1:end-1) + F_w.v).*b.DXI.rm(2:end-1,:);
            T_Ese(2:end-1,3:end)     = b.SE.v(2:end,:).*(F.v(2:end-1,2:end) + F_e.v).*b.DXI.rm(2:end-1,:);
            T_Wnw(3:end,1:end-2)     = -b.NW.v(1:end-1,:).*(F.v(2:end-1,1:end-1) + F_w.v).*b.DXI.rm(2:end-1,:);
            T_Ene(3:end,3:end)       = b.NE.v(1:end-1,:).*(F.v(2:end-1,2:end) + F_e.v).*b.DXI.rm(2:end-1,:);
            T_s(2:end-1,2:end-1)     = -(F.p(1:end-1,:) + F_s.p).*b.DETA.rm(2:end-1,:) + b.S.u(2:end,2:end-1).*F_C.u;
            T_n(3:end,2:end-1)       = (F.p(2:end,:) + F_n.p).*b.DETA.rm(2:end-1,:) + b.N.u(1:end-1,2:end-1).*F_C.u;
            T_s(2:end-1,3:end-1)     = T_s(2:end-1,3:end-1) - b.SE.v(2:end,1:end-1).*(F.v(2:end-1,2:end-1) + F_w.v(:,2:end)).*b.DXI.rm(2:end-1,2:end);
            T_s(2:end-1,2:end-2)     = T_s(2:end-1,2:end-2) + b.SW.v(2:end,2:end).*(F.v(2:end-1,2:end-1) + F_e.v(:,1:end-1)).*b.DXI.rm(2:end-1,1:end-1);
            T_n(3:end,3:end-1)       = T_n(3:end,3:end-1) - b.NE.v(1:end-1,1:end-1).*(F.v(2:end-1,2:end-1) + F_w.v(:,2:end)).*b.DXI.rm(2:end-1,2:end);
            T_n(3:end,2:end-2)       = T_n(3:end,2:end-2) + b.NW.v(1:end-1,2:end).*(F.v(2:end-1,2:end-1) + F_e.v(:,1:end-1)).*b.DXI.rm(2:end-1,1:end-1);
            %}
            
            %{
            T_Wsw(2:end-1,2:end-2)   = -(F.v(2:end-1,2:end-1) + F_w.v(:,2:end))/4.*b.DXI.rm(2:end-1,2:end);
            T_Ese(2:end-1,3:end-1)   = (F.v(2:end-1,2:end-1) + F_e.v(:,1:end-1))/4.*b.DXI.rm(2:end-1,1:end-1);
            T_Wnw(3:end,2:end-2)     = -(F.v(2:end-1,2:end-1) + F_w.v(:,2:end))/4.*b.DXI.rm(2:end-1,2:end);
            T_Ene(3:end,3:end-1)     = (F.v(2:end-1,2:end-1) + F_e.v(:,1:end-1))/4.*b.DXI.rm(2:end-1,1:end-1);
            T_s(2:end-1,2:end-1)     = T_Wsw(2:end-1,1:end-2) + T_Ese(2:end-1,3:end) - (F.p(1:end-1,:) + F_s.p).*b.DETA.rm(2:end-1,:) + F_C.u/2;
            T_n(3:end,2:end-1)       = T_Wnw(3:end,1:end-2) + T_Ene(3:end,3:end) + (F.p(2:end,:) + F_n.p).*b.DETA.rm(2:end-1,:) + F_C.u/2;

            T_Wsw(2:end-1,1)         = -(F.v(2:end-1,1) + F_w.v(:,1))/2.*b.DXI.rm(2:end-1,1);
            T_Ese(2:end-1,end)       = (F.v(2:end-1,end) + F_e.v(:,end))/2.*b.DXI.rm(2:end-1,end);
            T_Wnw(3:end,1)           = -(F.v(2:end-1,1) + F_w.v(:,1))/2.*b.DXI.rm(2:end-1,1);
            T_Ene(3:end,end)         = (F.v(2:end-1,end) + F_e.v(:,end))/2.*b.DXI.rm(2:end-1,end);
            %}
            
            T_Wsw=T_Wsw.'; T_Wsw=T_Wsw(:);
            T_s=T_s.'; T_s=T_s(:);
            T_Ese=T_Ese.'; T_Ese=T_Ese(:);
            T_Wnw=T_Wnw.'; T_Wnw=T_Wnw(:);
            T_n=T_n.'; T_n=T_n(:);
            T_Ene=T_Ene.'; T_Ene=T_Ene(:);

            b.Ja.r_m.T=spdiags([T_Wsw T_s T_Ese T_Wnw T_n T_Ene],[-1 0 1 b.J+1 b.J+2 b.J+3],(b.K+1)*(b.J+2),(b.K+1)*(b.J+2));
            b.Ja.r_m.T=[b.Ja.r_m.T sparse((b.K+1)*(b.J+2),b.J+2)];

    end
    
    b.Ja_counter.r_m.Rs1=[];
    b.Ja_counter.r_m.Rsend=[];
    
    if stability~=1
        index_bc_z=length(b.bc.z)/2;
    
        for o=1:2
            if max(strncmp('ss', {b.bc.z{(o-1)*index_bc_z+1:o*index_bc_z}},2))
                if o==1
                    % lower boundary
                        JA_rlb.u=b.JA_rlb1.u;
                        JA_rlb.w=b.JA_rlb1.w;
                        JA_rlb.p=b.JA_rlb1.p;
                        JA_rlb.T=b.JA_rlb1.T;
                        JA_rlb.v=b.JA_rlb1.v;

                        R_rlb.u=b.R_rlb1.u;
                        R_rlb.w=b.R_rlb1.w;
                        R_rlb.p=b.R_rlb1.p;
                        R_rlb.T=b.R_rlb1.T;
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
                        JA_rlb.T=b.JA_rlbend.T;
                        JA_rlb.v=b.JA_rlbend.v;

                        R_rlb.u=b.R_rlbend.u;
                        R_rlb.w=b.R_rlbend.w;
                        R_rlb.p=b.R_rlbend.p;
                        R_rlb.T=b.R_rlbend.T;
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
                    % [rho0*d((eta_z*u0^2-xi_z*u0*w0)*J^2*r^ax)/dRs]_delta.xi*delta.eta + [rho0*d((eta_z*u0*w0-xi_z*w0^2)*J^2*r^ax)/dRs]_delta.eta*delta.xi
                        FC.p = (ax*b.JA.p.^2.*R_rlb.p+2*b.JA.p.*JA_rlb.p.*b.R.p.^ax).*b.ETA_z.p.*rho(b.T.p).*b.u.p.^2 - (b.JA.p.^2.*(ax*R_rlb.p.*b.XI_z.p+XI_z_rlb_C.p.*b.R.p.^ax)+2*b.JA.p.*JA_rlb.p.*b.XI_z.p.*b.R.p.^ax).*rho(b.T.p).*b.u.p.*b.w.p;
                        FC.v = (ax*b.JA.v.^2.*R_rlb.v+2*b.JA.v.*JA_rlb.v.*b.R.v.^ax).*b.ETA_z.v.*rho(b.T.v).*b.u.v.*b.w.v - (b.JA.v.^2.*(ax*R_rlb.v.*b.XI_z.v+XI_z_rlb_C.v.*b.R.v.^ax)+2*b.JA.v.*JA_rlb.v.*b.XI_z.v.*b.R.v.^ax).*rho(b.T.v).*b.w.v.^2;

                    % [rho0*(-dxi_z/dRs*w0^2*J^2*r^ax]_delta.eta*delta.xi
                        FE.v = - b.JA.v.^2.*XI_z_rlb_E.v.*b.R.v.^ax.*rho(b.T.v).*b.w.v.^2;

                    % [rho0*(-dxi_z/dRs)*u0*w0*J^2*r^ax]_delta.xi*delta.eta
                        FE.p = - b.JA.p.^2.*XI_z_rlb_E.p.*b.R.p.^ax.*rho(b.T.p).*b.u.p.*b.w.p;

                else
                    FC.v = zeros(size(b.JA.v));
                    FE.v = FC.v;
                    FC.p = zeros(size(b.JA.p));
                    FE.p = FC.p;

                end

                % viscous terms with changed signs!
                if inviscid~=1
                    % -[mu*d(J*r^ax*(2*xi_r^2+xi_z^2)/dRs*d(J*eta_z*u0)/dxi]_delta.xi*delta.eta
                        FC.p = FC.p - mu(b.T.p).*((JA_rlb.p.*b.R.p.^ax + ax*b.JA.p.*R_rlb.p).*(2*b.XI_r.p.^2+b.XI_z.p.^2) + b.JA.p.*b.R.p.^ax.*(4*b.XI_r.p.*XI_r_rlb.p+2*b.XI_z.p.*XI_z_rlb_C.p)).*(b.JA.u(2:end,2:end-1).*b.ETA_z.u(2:end,2:end-1).*b.u.u(2:end,2:end-1)-b.JA.u(1:end-1,2:end-1).*b.ETA_z.u(1:end-1,2:end-1).*b.u.u(1:end-1,2:end-1))./b.DXI.c;
                        FE.p = FE.p - mu(b.T.p).*b.JA.p.*b.R.p.^ax.*2.*b.XI_z.p.*XI_z_rlb_E.p.*(b.JA.u(2:end,2:end-1).*b.ETA_z.u(2:end,2:end-1).*b.u.u(2:end,2:end-1)-b.JA.u(1:end-1,2:end-1).*b.ETA_z.u(1:end-1,2:end-1).*b.u.u(1:end-1,2:end-1))./b.DXI.c;

                    % -[mu*J*r^ax*(2*xi_r^2+xi_z^2)*d(dJ/dRs*eta_z*u0)/dxi]_delta.xi*delta.eta
                        FC.p = FC.p - mu(b.T.p).*b.JA.p.*b.R.p.^ax.*(2*b.XI_r.p.^2+b.XI_z.p.^2).*(JA_rlb.u(2:end,2:end-1).*b.ETA_z.u(2:end,2:end-1).*b.u.u(2:end,2:end-1)-JA_rlb.u(1:end-1,2:end-1).*b.ETA_z.u(1:end-1,2:end-1).*b.u.u(1:end-1,2:end-1))./b.DXI.c;

                    % [2/3*mu*dxi_r/dRs*d(J*r^ax*u0)/dxi]_delta.xi*delta.eta
                        FC.p = FC.p + 2/3*mu(b.T.p).*XI_r_rlb.p.*(b.JA.u(2:end,2:end-1).*b.R.u(2:end,2:end-1).^ax.*b.u.u(2:end,2:end-1)-b.JA.u(1:end-1,2:end-1).*b.R.u(1:end-1,2:end-1).^ax.*b.u.u(1:end-1,2:end-1))./b.DXI.c;

                    % [2/3*mu*xi_r*d(d(J*r^ax)/dRs*u0)/dxi]_delta.xi*delta.eta
                        FC.p = FC.p + 2/3*mu(b.T.p).*b.XI_r.p.*((JA_rlb.u(2:end,2:end-1).*b.R.u(2:end,2:end-1).^ax + ax*b.JA.u(2:end,2:end-1).*R_rlb.u(2:end,2:end-1)).*b.u.u(2:end,2:end-1)-(JA_rlb.u(1:end-1,2:end-1).*b.R.u(1:end-1,2:end-1).^ax + ax*b.JA.u(1:end-1,2:end-1).*R_rlb.u(1:end-1,2:end-1)).*b.u.u(1:end-1,2:end-1))./b.DXI.c;

                    % -[mu*d(J*r^ax*xi_z)/dRs*eta_z*d(J*eta_z*u0)/deta]_delta.xi*delta.eta
                        FC.p = FC.p - mu(b.T.p).*((JA_rlb.p.*b.R.p.^ax + ax*b.JA.p.*R_rlb.p).*b.XI_z.p + b.JA.p.*b.R.p.^ax.*XI_z_rlb_C.p).*b.ETA_z.p.*(b.JA.w(2:end-1,2:end).*b.ETA_z.w(2:end-1,2:end).*b.u.w(2:end-1,2:end)-b.JA.w(2:end-1,1:end-1).*b.ETA_z.w(2:end-1,1:end-1).*b.u.w(2:end-1,1:end-1))./b.DETA.c;
                        FE.p = FE.p - mu(b.T.p).*b.JA.p.*b.R.p.^ax.*XI_z_rlb_E.p.*b.ETA_z.p.*(b.JA.w(2:end-1,2:end).*b.ETA_z.w(2:end-1,2:end).*b.u.w(2:end-1,2:end)-b.JA.w(2:end-1,1:end-1).*b.ETA_z.w(2:end-1,1:end-1).*b.u.w(2:end-1,1:end-1))./b.DETA.c;

                    % -[mu*J*r^ax*xi_z*eta_z*d(dJ/dRs*eta_z*u0)/deta]_delta.xi*delta.eta
                        F_W.p = -mu(b.T.p).*b.JA.p.*b.R.p.^ax.*b.XI_z.p.*b.ETA_z.p.*JA_rlb.w(2:end-1,1:end-1).*b.ETA_z.w(2:end-1,1:end-1).*b.u.w(2:end-1,1:end-1)./b.DETA.c;
                        F_E.p = -mu(b.T.p).*b.JA.p.*b.R.p.^ax.*b.XI_z.p.*b.ETA_z.p.*JA_rlb.w(2:end-1,2:end).*b.ETA_z.w(2:end-1,2:end).*b.u.w(2:end-1,2:end)./b.DETA.c;

                    % -[mu*d(J*r^ax*xi_z)/dRs*eta_z*d(J*eta_z*u0)/dxi]_delta.eta*delta.xi
                        FC.v = FC.v - mu(b.T.v).*((JA_rlb.v.*b.R.v.^ax + ax*b.JA.v.*R_rlb.v).*b.XI_z.v + b.JA.v.*b.R.v.^ax.*XI_z_rlb_C.v).*b.ETA_z.v.*(b.JA.w(2:end,:).*b.ETA_z.w(2:end,:).*b.u.w(2:end,:)-b.JA.w(1:end-1,:).*b.ETA_z.w(1:end-1,:).*b.u.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];
                        FE.v = FE.v - mu(b.T.v).*b.JA.v.*b.R.v.^ax.*XI_z_rlb_E.v.*b.ETA_z.v.*(b.JA.w(2:end,:).*b.ETA_z.w(2:end,:).*b.u.w(2:end,:)-b.JA.w(1:end-1,:).*b.ETA_z.w(1:end-1,:).*b.u.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % -[mu*J*r^ax*xi_z*eta_z*d(dJ/dRs*eta_z*u0)/dxi]_delta.eta*delta.xi
                        FC.v = FC.v - mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.XI_z.v.*b.ETA_z.v.*(JA_rlb.w(2:end,:).*b.ETA_z.w(2:end,:).*b.u.w(2:end,:)-JA_rlb.w(1:end-1,:).*b.ETA_z.w(1:end-1,:).*b.u.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % -[mu*d(J*r^ax)/dRs*eta_z^2*d(J*eta_z*u0)/deta]_delta.eta*delta.xi
                        FC.v = FC.v - mu(b.T.v).*(JA_rlb.v.*b.R.v.^ax + ax*b.JA.v.*R_rlb.v).*b.ETA_z.v.^2.*(b.JA.u(:,2:end).*b.ETA_z.u(:,2:end).*b.u.u(:,2:end)-b.JA.u(:,1:end-1).*b.ETA_z.u(:,1:end-1).*b.u.u(:,1:end-1))./[b.DETA.zm(1,:); b.DETA.zm];

                    % -[mu*J*r^ax*eta_z^2*d(dJ/dRs*eta_z*u0)/deta]_delta.eta*delta.xi
                        F_W.v = - mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.ETA_z.v.^2.*JA_rlb.u(:,1:end-1).*b.ETA_z.u(:,1:end-1).*b.u.u(:,1:end-1)./[b.DETA.zm(1,:); b.DETA.zm];
                        F_E.v = - mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.ETA_z.v.^2.*JA_rlb.u(:,2:end).*b.ETA_z.u(:,2:end).*b.u.u(:,2:end)./[b.DETA.zm(1,:); b.DETA.zm];

                    % ax*2*mu*d(J^2/r)/dRs*eta_z*u0*delta.xi*delta.eta
                        F.u = ax*2*mu(b.T.u).*(2*b.JA.u.*JA_rlb.u./b.R.u.^ax-b.JA.u.^2.*b.R.u.^(-2).*R_rlb.u).*b.ETA_z.u.*b.u.u;
                        FC.u = F.u(2:end-1,2:end-1).*b.DXI.rm(2:end-1,:).*b.DETA.rm(2:end-1,:);

                    % -ax*2/3*mu*d(1/r)/dRs*d(J*r*u0)/dxi*delta.xi*delta.eta
                        F.u = -ax*2/3*mu(b.T.u).*(-b.R.u.^(-2).*R_rlb.u).*(b.JA.T(2:end,:).*b.R.T(2:end,:).*b.u.T(2:end,:)-b.JA.T(1:end-1,:).*b.R.T(1:end-1,:).*b.u.T(1:end-1,:));
                        FC.u = FC.u + F.u(2:end-1,2:end-1).*b.DETA.rm(2:end-1,:);

                    % -ax*2/3*mu/r*d(d(J*r)/dRs*u0)/dxi*delta.xi*delta.eta
                        F.u = -ax*2/3*mu(b.T.u)./b.R.u.^ax.*((JA_rlb.T(2:end,:).*b.R.T(2:end,:) + b.JA.T(2:end,:).*R_rlb.T(2:end,:)).*b.u.T(2:end,:)-(JA_rlb.T(1:end-1,:).*b.R.T(1:end-1,:) + b.JA.T(1:end-1,:).*R_rlb.T(1:end-1,:)).*b.u.T(1:end-1,:));
                        FC.u = FC.u + F.u(2:end-1,2:end-1).*b.DETA.rm(2:end-1,:);

                    % [2/3*mu*dxi_r/dRs*d(J*r^ax*w0)/deta]_delta.xi*delta.eta
                        FC.p = FC.p + 2/3*mu(b.T.p).*XI_r_rlb.p.*(b.JA.w(2:end-1,2:end).*b.R.w(2:end-1,2:end).^ax.*b.w.w(2:end-1,2:end)-b.JA.w(2:end-1,1:end-1).*b.R.w(2:end-1,1:end-1).^ax.*b.w.w(2:end-1,1:end-1))./b.DETA.c;

                    % [2/3*mu*xi_r*d(d(J*r^ax)/dRs*w0)/deta]_delta.xi*delta.eta
                        F_W.p = F_W.p + 2/3*mu(b.T.p).*b.XI_r.p.*((JA_rlb.w(2:end-1,1:end-1).*b.R.w(2:end-1,1:end-1).^ax + ax*b.JA.w(2:end-1,1:end-1).*R_rlb.w(2:end-1,1:end-1)).*b.w.w(2:end-1,1:end-1))./b.DETA.c;
                        F_E.p = F_E.p + 2/3*mu(b.T.p).*b.XI_r.p.*((JA_rlb.w(2:end-1,2:end).*b.R.w(2:end-1,2:end).^ax + ax*b.JA.w(2:end-1,2:end).*R_rlb.w(2:end-1,2:end)).*b.w.w(2:end-1,2:end))./b.DETA.c;

                    % [mu*d(J*r^ax*xi_z)/dRs*eta_z*d(J*xi_z*w0)/deta]_delta.xi*delta.eta
                        FC.p = FC.p + mu(b.T.p).*((JA_rlb.p.*b.R.p.^ax + ax*b.JA.p.*R_rlb.p).*b.XI_z.p + b.JA.p.*b.R.p.^ax.*XI_z_rlb_C.p).*b.ETA_z.p.*(b.JA.w(2:end-1,2:end).*b.XI_z.w(2:end-1,2:end).*b.w.w(2:end-1,2:end)-b.JA.w(2:end-1,1:end-1).*b.XI_z.w(2:end-1,1:end-1).*b.w.w(2:end-1,1:end-1))./b.DETA.c;
                        FE.p = FE.p + mu(b.T.p).*b.JA.p.*b.R.p.^ax.*XI_z_rlb_E.p.*b.ETA_z.p.*(b.JA.w(2:end-1,2:end).*b.XI_z.w(2:end-1,2:end).*b.w.w(2:end-1,2:end)-b.JA.w(2:end-1,1:end-1).*b.XI_z.w(2:end-1,1:end-1).*b.w.w(2:end-1,1:end-1))./b.DETA.c;

                    % [mu*J*r^ax*xi_z*eta_z*d(d(J*xi_z)/dRs*w0)/deta]_delta.xi*delta.eta
                        F_W.p = F_W.p + mu(b.T.p).*b.JA.p.*b.R.p.^ax.*b.XI_z.p.*b.ETA_z.p.*(JA_rlb.w(2:end-1,1:end-1).*b.XI_z.w(2:end-1,1:end-1) + b.JA.w(2:end-1,1:end-1).*XI_z_rlb_C.w(2:end-1,1:end-1)).*b.w.w(2:end-1,1:end-1)./b.DETA.c;
                        F_E.p = F_E.p + mu(b.T.p).*b.JA.p.*b.R.p.^ax.*b.XI_z.p.*b.ETA_z.p.*(JA_rlb.w(2:end-1,2:end).*b.XI_z.w(2:end-1,2:end) + b.JA.w(2:end-1,2:end).*XI_z_rlb_C.w(2:end-1,2:end)).*b.w.w(2:end-1,2:end)./b.DETA.c;
                        FE_W.p = mu(b.T.p).*b.JA.p.*b.R.p.^ax.*b.XI_z.p.*b.ETA_z.p.*b.JA.w(2:end-1,1:end-1).*XI_z_rlb_E.w(2:end-1,1:end-1).*b.w.w(2:end-1,1:end-1)./b.DETA.c;
                        FE_E.p = mu(b.T.p).*b.JA.p.*b.R.p.^ax.*b.XI_z.p.*b.ETA_z.p.*b.JA.w(2:end-1,2:end).*XI_z_rlb_E.w(2:end-1,2:end).*b.w.w(2:end-1,2:end)./b.DETA.c;

                    % [mu*d(J*r^ax*(2*xi_r^2+xi_z^2)/dRs*d(J*xi_z*w0)/dxi]_delta.xi*delta.eta
                        FC.p = FC.p + mu(b.T.p).*((JA_rlb.p.*b.R.p.^ax + ax*b.JA.p.*R_rlb.p).*(2*b.XI_r.p.^2+b.XI_z.p.^2) + b.JA.p.*b.R.p.^ax.*(4*b.XI_r.p.*XI_r_rlb.p+2*b.XI_z.p.*XI_z_rlb_C.p)).*(b.JA.u(2:end,2:end-1).*b.XI_z.u(2:end,2:end-1).*b.w.u(2:end,2:end-1)-b.JA.u(1:end-1,2:end-1).*b.XI_z.u(1:end-1,2:end-1).*b.w.u(1:end-1,2:end-1))./b.DXI.c;
                        FE.p = FE.p + mu(b.T.p).*b.JA.p.*b.R.p.^ax.*2.*b.XI_z.p.*XI_z_rlb_E.p.*(b.JA.u(2:end,2:end-1).*b.XI_z.u(2:end,2:end-1).*b.w.u(2:end,2:end-1)-b.JA.u(1:end-1,2:end-1).*b.XI_z.u(1:end-1,2:end-1).*b.w.u(1:end-1,2:end-1))./b.DXI.c;

                    % [mu*J*r^ax*(2*xi_r^2+xi_z^2)*d(d(J*xi_z)/dRs*w0)/dxi]_delta.xi*delta.eta
                        FC.p = FC.p + mu(b.T.p).*b.JA.p.*b.R.p.^ax.*(2*b.XI_r.p.^2+b.XI_z.p.^2).*((JA_rlb.u(2:end,2:end-1).*b.XI_z.u(2:end,2:end-1) + b.JA.u(2:end,2:end-1).*XI_z_rlb_C.u(2:end,2:end-1)).*b.w.u(2:end,2:end-1)-(JA_rlb.u(1:end-1,2:end-1).*b.XI_z.u(1:end-1,2:end-1) + b.JA.u(1:end-1,2:end-1).*XI_z_rlb_C.u(1:end-1,2:end-1)).*b.w.u(1:end-1,2:end-1))./b.DXI.c;
                        FE.p = FE.p + mu(b.T.p).*b.JA.p.*b.R.p.^ax.*(2*b.XI_r.p.^2+b.XI_z.p.^2).*(b.JA.u(2:end,2:end-1).*XI_z_rlb_E.u(2:end,2:end-1).*b.w.u(2:end,2:end-1)-b.JA.u(1:end-1,2:end-1).*XI_z_rlb_E.u(1:end-1,2:end-1).*b.w.u(1:end-1,2:end-1))./b.DXI.c;

                     % [mu*d(J*r^ax*xi_r*xi_z)/dRs*d(J*xi_r*w0)/dxi]_delta.xi*delta.eta
                        FC.p = FC.p - mu(b.T.p).*((JA_rlb.p.*b.R.p.^ax + ax*b.JA.p.*R_rlb.p).*b.XI_r.p.*b.XI_z.p + b.JA.p.*b.R.p.^ax.*(XI_r_rlb.p.*b.XI_z.p+b.XI_r.p.*XI_z_rlb_C.p)).*(b.JA.u(2:end,2:end-1).*b.XI_r.u(2:end,2:end-1).*b.w.u(2:end,2:end-1)-b.JA.u(1:end-1,2:end-1).*b.XI_r.u(1:end-1,2:end-1).*b.w.u(1:end-1,2:end-1))./b.DXI.c;
                        FE.p = FE.p - mu(b.T.p).*b.JA.p.*b.R.p.^ax.*b.XI_r.p.*XI_z_rlb_E.p.*(b.JA.u(2:end,2:end-1).*b.XI_r.u(2:end,2:end-1).*b.w.u(2:end,2:end-1)-b.JA.u(1:end-1,2:end-1).*b.XI_r.u(1:end-1,2:end-1).*b.w.u(1:end-1,2:end-1))./b.DXI.c;

                    % [mu*J*r^ax*xi_r*xi_z*d(d(J*xi_r)/dRs*w0)/dxi]_delta.xi*delta.eta
                        FC.p = FC.p - mu(b.T.p).*b.JA.p.*b.R.p.^ax.*b.XI_r.p.*b.XI_z.p.*((JA_rlb.u(2:end,2:end-1).*b.XI_r.u(2:end,2:end-1) + b.JA.u(2:end,2:end-1).*XI_r_rlb.u(2:end,2:end-1)).*b.w.u(2:end,2:end-1)-(JA_rlb.u(1:end-1,2:end-1).*b.XI_r.u(1:end-1,2:end-1) + b.JA.u(1:end-1,2:end-1).*XI_r_rlb.u(1:end-1,2:end-1)).*b.w.u(1:end-1,2:end-1))./b.DXI.c;

                    % -[mu*d(J*r^ax*xi_r)/dRs*eta_z*d(J*xi_r*w0)/dxi]_delta.eta*delta.xi
                        FC.v = FC.v - mu(b.T.v).*((JA_rlb.v.*b.R.v.^ax + ax*b.JA.v.*R_rlb.v).*b.XI_r.v + b.JA.v.*b.R.v.^ax.*XI_r_rlb.v).*b.ETA_z.v.*(b.JA.w(2:end,:).*b.XI_r.w(2:end,:).*b.w.w(2:end,:)-b.JA.w(1:end-1,:).*b.XI_r.w(1:end-1,:).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % -[mu*J*r^ax*xi_r*eta_z*d(dJ/dRs*xi_r*w0)/dxi]_delta.eta*delta.xi
                        FC.v = FC.v - mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.XI_r.v.*b.ETA_z.v.*((JA_rlb.w(2:end,:).*b.XI_r.w(2:end,:) + b.JA.w(2:end,:).*XI_r_rlb.w(2:end,:)).*b.w.w(2:end,:)-(JA_rlb.w(1:end-1,:).*b.XI_r.w(1:end-1,:) + b.JA.w(1:end-1,:).*XI_r_rlb.w(1:end-1,:)).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % -[mu*d(J*r^ax*xi_z)/dRs*eta_z*d(J*xi_z*w0)/dxi]_delta.eta*delta.xi
                        FC.v = FC.v + mu(b.T.v).*((JA_rlb.v.*b.R.v.^ax + ax*b.JA.v.*R_rlb.v).*b.XI_z.v + b.JA.v.*b.R.v.^ax.*XI_z_rlb_C.v).*b.ETA_z.v.*(b.JA.w(2:end,:).*b.XI_z.w(2:end,:).*b.w.w(2:end,:)-b.JA.w(1:end-1,:).*b.XI_z.w(1:end-1,:).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];
                        FE.v = FE.v + mu(b.T.v).*b.JA.v.*b.R.v.^ax.*XI_z_rlb_E.v.*b.ETA_z.v.*(b.JA.w(2:end,:).*b.XI_z.w(2:end,:).*b.w.w(2:end,:)-b.JA.w(1:end-1,:).*b.XI_z.w(1:end-1,:).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % -[mu*J*r^ax*xi_z*eta_z*d(d(J*xi_z)/dRs*w0)/dxi]_delta.eta*delta.xi
                        FC.v = FC.v + mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.XI_z.v.*b.ETA_z.v.*((JA_rlb.w(2:end,:).*b.XI_z.w(2:end,:) + b.JA.w(2:end,:).*XI_z_rlb_C.w(2:end,:)).*b.w.w(2:end,:)-(JA_rlb.w(1:end-1,:).*b.XI_z.w(1:end-1,:) + b.JA.w(1:end-1,:).*XI_z_rlb_C.w(1:end-1,:)).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];
                        FE.v = FE.v + mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.XI_z.v.*b.ETA_z.v.*(b.JA.w(2:end,:).*XI_z_rlb_E.w(2:end,:).*b.w.w(2:end,:)-b.JA.w(1:end-1,:).*XI_z_rlb_E.w(1:end-1,:).*b.w.w(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm];

                    % [mu*d(J*r^ax)/dRs*eta_z^2*d(J*xi_z*w0)/deta]_delta.eta*delta.xi
                        FC.v = FC.v + mu(b.T.v).*(JA_rlb.v.*b.R.v.^ax + ax*b.JA.v.*R_rlb.v).*b.ETA_z.v.^2.*(b.JA.u(:,2:end).*b.XI_z.u(:,2:end).*b.w.u(:,2:end)-b.JA.u(:,1:end-1).*b.XI_z.u(:,1:end-1).*b.w.u(:,1:end-1))./[b.DETA.zm(1,:); b.DETA.zm];

                    % [mu*J*r^ax*eta_z^2*d(dJ/dRs*xi_z*w0)/deta]_delta.eta*delta.xi
                        F_W.v = F_W.v + mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.ETA_z.v.^2.*(JA_rlb.u(:,1:end-1).*b.XI_z.u(:,1:end-1) + b.JA.u(:,1:end-1).*XI_z_rlb_C.u(:,1:end-1)).*b.w.u(:,1:end-1)./[b.DETA.zm(1,:); b.DETA.zm];
                        F_E.v = F_E.v + mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.ETA_z.v.^2.*(JA_rlb.u(:,2:end).*b.XI_z.u(:,2:end) + b.JA.u(:,2:end).*XI_z_rlb_C.u(:,2:end)).*b.w.u(:,2:end)./[b.DETA.zm(1,:); b.DETA.zm];
                        FE_W.v = mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.ETA_z.v.^2.*b.JA.u(:,1:end-1).*XI_z_rlb_E.u(:,1:end-1).*b.w.u(:,1:end-1)./[b.DETA.zm(1,:); b.DETA.zm];
                        FE_E.v = mu(b.T.v).*b.JA.v.*b.R.v.^ax.*b.ETA_z.v.^2.*b.JA.u(:,2:end).*XI_z_rlb_E.u(:,2:end).*b.w.u(:,2:end)./[b.DETA.zm(1,:); b.DETA.zm];

                    % -ax*2*mu*d(J^2/r*xi_z)/dRs*w0*delta.xi*delta.eta
                        F.u = -ax*2*mu(b.T.u).*((2*b.JA.u.*JA_rlb.u./b.R.u.^ax-b.JA.u.^2.*b.R.u.^(-2).*R_rlb.u).*b.XI_z.u + b.JA.u.^2./b.R.u.^ax.*XI_z_rlb_C.u).*b.w.u;
                        FC.u = FC.u + F.u(2:end-1,2:end-1).*b.DXI.rm(2:end-1,:).*b.DETA.rm(2:end-1,:);
                        F.u = -ax*2*mu(b.T.u).*b.JA.u.^2./b.R.u.^ax.*XI_z_rlb_E.u.*b.w.u;
                        FE.u = F.u(2:end-1,2:end-1).*b.DXI.rm(2:end-1,:).*b.DETA.rm(2:end-1,:);

                    % -ax*2/3*mu*d(1/r)/dRs*d(J*r*w0)/deta*delta.xi*delta.eta
                        F.u = -ax*2/3*mu(b.T.u).*(-b.R.u.^(-2).*R_rlb.u);
                        FC.u = FC.u + F.u(2:end-1,2:end-1).*(b.JA.v(2:end-1,2:end).*b.R.v(2:end-1,2:end).*b.w.v(2:end-1,2:end)-b.JA.v(2:end-1,1:end-1).*b.R.v(2:end-1,1:end-1).*b.w.v(2:end-1,1:end-1)).*b.DXI.rm(2:end-1,:);

                    % -ax*2/3*mu/r*d(d(J*r)/dRs*w0)/deta*delta.xi*delta.eta
                        F.u = -ax*2/3*mu(b.T.u)./b.R.u.^ax;
                        F_W.u = F.u(2:end-1,2:end-1).*(JA_rlb.v(2:end-1,1:end-1).*b.R.v(2:end-1,1:end-1) + b.JA.v(2:end-1,1:end-1).*R_rlb.v(2:end-1,1:end-1)).*b.w.v(2:end-1,1:end-1).*b.DXI.rm(2:end-1,:);
                        F_E.u = F.u(2:end-1,2:end-1).*(JA_rlb.v(2:end-1,2:end).*b.R.v(2:end-1,2:end) + b.JA.v(2:end-1,2:end).*R_rlb.v(2:end-1,2:end)).*b.w.v(2:end-1,2:end).*b.DXI.rm(2:end-1,:);


                else
                    F_W.p  = zeros(size(b.JA.p));
                    F_E.p  = F_W.p;
                    F_W.v  = zeros(size(b.JA.v));
                    F_E.v  = F_W.v;
                    FE_W.p = F_W.p;
                    FE_E.p = F_W.p;
                    FE_W.v = F_W.v;
                    FE_E.v = F_W.v;
                    FC.u   = zeros(size(b.JA.u(2:end-1,2:end-1)));
                    FE.u   = FC.u;
                    F_W.u  = FC.u;
                    F_E.u  = FC.u;

                end

                % pressure terms with changed signs!
                    % (d(d(J*xi_r)/dRs*p0)/delta.xi*r^ax+ax*d(J*xi_r*p0)/delta.xi*dr/dRs)*delta.xi*delta.eta
                        F.p  = b.JA.p.*b.XI_r.p.*b.p.p;
                        FC.u = FC.u + ax*(F.p(2:end,:)-F.p(1:end-1,:)).*R_rlb.u(2:end-1,2:end-1).*b.DETA.rm(2:end-1,:);
                        F.p  = (b.JA.p.*XI_r_rlb.p+JA_rlb.p.*b.XI_r.p).*b.p.p;
                        FC.u = FC.u + (F.p(2:end,:)-F.p(1:end-1,:)).*b.R.u(2:end-1,2:end-1).^ax.*b.DETA.rm(2:end-1,:);

                Rs_W=zeros(b.K+1,b.J); Rs_C=Rs_W; Rs_E=Rs_W; Rs_ad_W=Rs_W; Rs_ad_E=Rs_W; Rs_WW=Rs_W; Rs_EE=Rs_W; Rs_ad_WW=Rs_W; Rs_ad_EE=Rs_W;

                Rs_ad_WW(2:end-1,:) = -[zeros(b.K-1,2) b.E.w(3:end-1,1:end-2)].*FE_W.v(2:end-1,1:end-1).*b.DXI.rm(2:end-1,:);
                Rs_ad_EE(2:end-1,:) = [b.W.w(3:end-1,3:end) zeros(b.K-1,2)].*FE_E.v(2:end-1,2:end).*b.DXI.rm(2:end-1,:);
                Rs_WW(2:end-1,:)    = -[zeros(b.K-1,1) b.W.w(3:end-1,1:end-1)].*FE_W.v(2:end-1,1:end-1).*b.DXI.rm(2:end-1,:);
                Rs_W(2:end-1,:)     = Rs_ad_WW(2:end-1,:) + (-b.W.w(3:end-1,:).*FC.v(2:end-1,1:end-1) + FE.v(2:end-1,1:end-1) + F_W.v(2:end-1,1:end-1)).*b.DXI.rm(2:end-1,:) - b.W.w(3:end-1,:).*(FE.p(2:end,:)-FE.p(1:end-1,:) + F_W.p(2:end,:)-F_W.p(1:end-1,:)).*b.DETA.rm(2:end-1,:) - b.W.w(3:end-1,:).*(FE.u + F_W.u) + (FE_W.p(2:end,:) - FE_W.p(1:end-1,:)).*b.DETA.rm(2:end-1,:) + ([zeros(b.K-1,1) b.W.w(3:end-1,2:end)].*(FE_W.v(2:end-1,1:end-1) + FE_W.v(2:end-1,2:end) + FE_E.v(2:end-1,1:end-1))).*b.DXI.rm(2:end-1,:);
                Rs_E(2:end-1,:)     = Rs_ad_EE(2:end-1,:) + (b.E.w(3:end-1,:).*FC.v(2:end-1,2:end) + FE.v(2:end-1,2:end) + F_E.v(2:end-1,2:end)).*b.DXI.rm(2:end-1,:) + b.E.w(3:end-1,:).*(FE.p(2:end,:)-FE.p(1:end-1,:) + F_E.p(2:end,:)-F_E.p(1:end-1,:)).*b.DETA.rm(2:end-1,:) + b.E.w(3:end-1,:).*(FE.u + F_E.u) + (FE_E.p(2:end,:) - FE_E.p(1:end-1,:)).*b.DETA.rm(2:end-1,:) - ([b.E.w(3:end-1,1:end-1) zeros(b.K-1,1)].*(FE_E.v(2:end-1,1:end-1) + FE_E.v(2:end-1,2:end) + FE_W.v(2:end-1,2:end))).*b.DXI.rm(2:end-1,:);
                Rs_EE(2:end-1,:)    = [b.E.w(3:end-1,2:end) zeros(b.K-1,1)].*FE_E.v(2:end-1,2:end).*b.DXI.rm(2:end-1,:);
                Rs_ad_W(2:end-1,:)  = [zeros(b.K-1,1) -b.E.w(3:end-1,1:end-1)].*(FC.v(2:end-1,1:end-1).*b.DXI.rm(2:end-1,:) + (FE.p(2:end,:)-FE.p(1:end-1,:) + F_W.p(2:end,:)-F_W.p(1:end-1,:)).*b.DETA.rm(2:end-1,:) + FE.u + F_W.u) + ([zeros(b.K-1,1) b.E.w(3:end-1,1:end-1)].*(FE_W.v(2:end-1,2:end) + FE_E.v(2:end-1,1:end-1)) + [ones(b.K-1,1) b.E.w(3:end-1,1:end-1)].*FE_W.v(2:end-1,1:end-1)).*b.DXI.rm(2:end-1,:);
                Rs_ad_E(2:end-1,:)  = [b.W.w(3:end-1,2:end) zeros(b.K-1,1)].*(FC.v(2:end-1,2:end).*b.DXI.rm(2:end-1,:) + (FE.p(2:end,:)-FE.p(1:end-1,:) + F_E.p(2:end,:)-F_E.p(1:end-1,:)).*b.DETA.rm(2:end-1,:) + FE.u + F_E.u) - ([b.W.w(3:end-1,2:end) zeros(b.K-1,1)].*(FE_E.v(2:end-1,1:end-1) + FE_W.v(2:end-1,2:end)) + [b.W.w(3:end-1,2:end) ones(b.K-1,1)].*FE_E.v(2:end-1,2:end)).*b.DXI.rm(2:end-1,:);
                Rs_C(2:end-1,:)     = Rs_ad_W(2:end-1,:) + Rs_ad_E(2:end-1,:)  - (FE.v(2:end-1,1:end-1) + FE.v(2:end-1,2:end) + F_W.v(2:end-1,2:end) + F_E.v(2:end-1,1:end-1)).*b.DXI.rm(2:end-1,:) + (FC.p(2:end,:)-FC.p(1:end-1,:) - FE_W.p(2:end,:)+FE_W.p(1:end-1,:) - FE_E.p(2:end,:)+FE_E.p(1:end-1,:)).*b.DETA.rm(2:end-1,:) + FC.u;

                %Rs_W=Rs_W(:);
                %Rs_WW=Rs_WW(:);
                %Rs_C=Rs_C(:);
                %Rs_E=Rs_E(:);
                %Rs_EE=Rs_EE(:);

                index=kron(ones(1,b.J),2:b.J+2:(b.K+1)*(b.J+2))+ kron(0:(b.K+1)*(b.J+2)+1:(b.K+1)*(b.J+2)*(b.J),ones(1,(b.K+1)));
                Rs1=sparse((b.K+1)*(b.J+2),b.J); Rs2=Rs1; Rs3=Rs1; Rs4=Rs1; Rs5=Rs1;
                Rs1(index)=Rs_W;
                Rs2(index)=Rs_C;
                Rs3(index)=Rs_E;
                Rs4(index)=Rs_WW;
                Rs5(index)=Rs_EE;

                Ja.r_m.Rs_all=[sparse((b.K+1)*(b.J+2),1) Rs1 sparse((b.K+1)*(b.J+2),3)] + [sparse((b.K+1)*(b.J+2),2) Rs2 sparse((b.K+1)*(b.J+2),2)] + [sparse((b.K+1)*(b.J+2),3) Rs3 sparse((b.K+1)*(b.J+2),1)] + [Rs4 sparse((b.K+1)*(b.J+2),4)] + [sparse((b.K+1)*(b.J+2),4) Rs5];
                Ja.r_m.Rs_all=Ja.r_m.Rs_all(:,2:end-1);
                Ja.r_m.Rs=[];
                Ja_counter.r_m.Rs=[];

                for n=1:length(b.bc.z)/2
                    if (o==1 && strcmp(b.bc.z{n}(1),'s')==1 && strcmp(b.bc.z{n}(2),'s')==1) || (o==2 && strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'s')==1 && strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'s')==1)
                        Ja.r_m.Rs_part=Ja.r_m.Rs_all(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                        if (n>1 && strcmp(b.bc.z{(o-1)*index_bc_z+n-1}(1),'s')~=1 && strcmp(b.bc.z{(o-1)*index_bc_z+n-1}(2),'s')~=1)
                            Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n))+b.J+1:b.J+2:end-b.J-3,1)=0;
                            Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n))+b.J+2:b.J+2:end-b.J-2,1)=(FE.v(2:end-1,b.z.w==b.geom.z(n)) + F_E.v(2:end-1,b.z.w==b.geom.z(n))).*b.DXI.rm(2:end-1,1) + (FE_E.p(2:end,find(b.z.w==b.geom.z(n))-1) - FE_E.p(1:end-1,find(b.z.w==b.geom.z(n))-1)).*b.DETA.rm(2:end-1,find(b.z.w==b.geom.z(n))-1) + Rs_ad_EE(2:end-1,find(b.z.w==b.geom.z(n))-1);
                            Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n))+1:b.J+2:end,1)=Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n))+1:b.J+2:end,1)-Rs_ad_W(:,b.z.w==b.geom.z(n));
                            Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n))+b.J+4:b.J+2:end-b.J,1)=Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n))+b.J+4:b.J+2:end-b.J,1) - Rs_ad_WW(2:end-1,find(b.z.w==b.geom.z(n))+1);

                        end
                        %Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n))+b.J+2:b.J+2:end-b.J-2,1)=-b.w.u(2:end-1,1).*XI_z_rlb_E.u(2:end-1,1)./b.ETA_z.u(2:end-1,1);

                        if (n<length(b.bc.z)/2 && strcmp(b.bc.z{(o-1)*index_bc_z+n+1}(1),'s')~=1 && strcmp(b.bc.z{(o-1)*index_bc_z+n+1}(2),'s')~=1)
                            Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n+1))+b.J+3:b.J+2:end-b.J-1,end)=(FE.v(2:end-1,b.z.w==b.geom.z(n+1)) + F_W.v(2:end-1,b.z.w==b.geom.z(n+1))).*b.DXI.rm(2:end-1,1) + (FE_W.p(2:end,b.z.w==b.geom.z(n+1)) - FE_W.p(1:end-1,b.z.w==b.geom.z(n+1))).*b.DETA.rm(2:end-1,b.z.w==b.geom.z(n+1)) + Rs_ad_WW(2:end-1,b.z.w==b.geom.z(n+1));
                            Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n+1))+b.J+4:b.J+2:end-b.J,end)=0;
                            Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n+1)):b.J+2:end,end)=Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n+1)):b.J+2:end,end)-Rs_ad_E(:,find(b.z.w==b.geom.z(n+1))-1);
                            Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n+1))+b.J+1:b.J+2:end-b.J-3,end)=Ja.r_m.Rs_part(find(b.z.w==b.geom.z(n+1))+b.J+1:b.J+2:end-b.J-3,end) - Rs_ad_EE(2:end-1,find(b.z.w==b.geom.z(n+1))-2);

                        end

                        if eval(['b.' b.bc.z{(o-1)*index_bc_z+n} '.o==1'])
                            Ja.r_m.Rs=[Ja.r_m.Rs Ja.r_m.Rs_part];

                        elseif eval(['b.' b.bc.z{(o-1)*index_bc_z+n} '.o==0'])
                            Ja_counter.r_m.Rs=[Ja_counter.r_m.Rs Ja.r_m.Rs_part];

                        end

                    end

                end

                if o==1
                    b.Ja.r_m.Rs1=Ja.r_m.Rs;
                    b.Ja_counter.r_m.Rs1=Ja_counter.r_m.Rs;

                else
                    b.Ja.r_m.Rsend=Ja.r_m.Rs;
                    b.Ja_counter.r_m.Rsend=Ja_counter.r_m.Rs;

                end

            end
            
        end

    end