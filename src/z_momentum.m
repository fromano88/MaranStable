function[b]=z_momentum(b, flowopt)
    ax=flowopt.ax;
    g=flowopt.g;
    inviscid=flowopt.inviscid;
    creeping=flowopt.creeping;
    rho=str2func(b.rho);
    mu=str2func(b.mu);
    
    b.z_m_inert.w=sparse((b.K+2)*(b.J+1),(b.K+2)*(b.J+1));
    b.z_m_visc.u=sparse((b.K+2)*(b.J+1),(b.K+1)*(b.J+2));
    b.z_m_visc.w=sparse((b.K+2)*(b.J+1),(b.K+2)*(b.J+1));
        
    if creeping~=1
    % inertia terms for w
    %  [J^2*r^ax*rho0*xi_r*u0*w0]_delta.xi*delta.eta +[J^2*r^ax*rho0*xi_r*w0^2]_delta.eta*delta.xi
        F.p=b.JA.p.^2.*b.R.p.^ax.*rho(b.T.p).*b.XI_r.p.*b.w.p;
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

        b.z_m_inert.w=spdiags([w_S w_W w_C w_E w_N],[-b.J-1 -1 0 1 b.J+1],(b.K+2)*(b.J+1),(b.K+2)*(b.J+1));
        
    end
    
    if inviscid~=1
    % viscous terms for u
        F1.v=mu(b.T.v).*b.JA.v.*b.XI_r.v.*b.XI_z.v.*b.R.v.^ax;
        F2.v=-2/3*mu(b.T.v).*b.XI_z.v;
        F3.v=mu(b.T.v).*b.JA.v.*b.XI_r.v.*b.ETA_z.v.*b.R.v.^ax;
        
        F1.p=-2/3*mu(b.T.p).*b.ETA_z.p;
        
        u_sw=zeros(b.K+2,b.J+2); u_se=u_sw; u_nw=u_sw; u_ne=u_sw; u_SSw=u_sw; u_SSe=u_sw; u_NNw=u_sw; u_NNe=u_sw;

        C=-0.25*((F1.v(1:end-1,2:end-1)./b.DXI.rm(1:end-1,1:end-1) + F1.v(2:end,2:end-1)./b.DXI.rm(2:end,1:end-1)).*b.ETA_z.w(2:end-1,2:end-1) + (F2.v(1:end-1,2:end-1)./b.DXI.rm(1:end-1,1:end-1) + F2.v(2:end,2:end-1)./b.DXI.rm(2:end,1:end-1)).*b.R.w(2:end-1,2:end-1).^ax).*b.JA.w(2:end-1,2:end-1).*b.DETA.zm(:,2:end-1);
        
        u_sw(1:end-2,2:end-2) = C + (F3.v(1:end-1,2:end-1).*b.ETA_z.u(1:end-1,2:end-2) + F1.p(:,1:end-1).*b.R.u(1:end-1,2:end-2).^ax).*b.JA.u(1:end-1,2:end-2);
        u_se(1:end-2,3:end-1) = C + -(F3.v(1:end-1,2:end-1).*b.ETA_z.u(1:end-1,3:end-1) + F1.p(:,2:end).*b.R.u(1:end-1,3:end-1).^ax).*b.JA.u(1:end-1,3:end-1);
        u_nw(2:end-1,2:end-2) = C + -(F3.v(2:end,2:end-1).*b.ETA_z.u(2:end,2:end-2) + F1.p(:,1:end-1).*b.R.u(2:end,2:end-2).^ax).*b.JA.u(2:end,2:end-2);
        u_ne(2:end-1,3:end-1) = C + (F3.v(2:end,2:end-1).*b.ETA_z.u(2:end,3:end-1) + F1.p(:,2:end).*b.R.u(2:end,3:end-1).^ax).*b.JA.u(2:end,3:end-1);

        u_SSw(1:end-3,2:end-2) = 0.25*(F1.v(2:end-1,2:end-1).*b.ETA_z.w(2:end-2,2:end-1) + F2.v(2:end-1,2:end-1).*b.R.w(2:end-2,2:end-1).^ax).*b.JA.w(2:end-2,2:end-1)./b.DXI.rm(2:end-1,1:end-1).*b.DETA.zm(2:end,2:end-1);
        u_SSe(1:end-3,3:end-1) = u_SSw(1:end-3,2:end-2);
        u_NNw(3:end-1,2:end-2) = 0.25*(F1.v(2:end-1,2:end-1).*b.ETA_z.w(3:end-1,2:end-1) + F2.v(2:end-1,2:end-1).*b.R.w(3:end-1,2:end-1).^ax).*b.JA.w(3:end-1,2:end-1)./b.DXI.rm(2:end-1,1:end-1).*b.DETA.zm(1:end-1,2:end-1);
        u_NNe(3:end-1,3:end-1) = u_NNw(3:end-1,2:end-2);
        u_sw(1,2:end-2)        = u_sw(1,2:end-2) + 0.5*(F1.v(1,2:end-1).*b.ETA_z.w(1,2:end-1) + F2.v(1,2:end-1).*b.R.w(1,2:end-1).^ax).*b.JA.w(1,2:end-1)./b.DXI.rm(1,1:end-1).*b.DETA.zm(1,2:end-1);
        u_sw(2:end-2,2:end-2)  = u_sw(2:end-2,2:end-2) + u_SSw(1:end-3,2:end-2);
        u_se(1,3:end-1)        = u_se(1,3:end-1) + 0.5*(F1.v(1,2:end-1).*b.ETA_z.w(1,2:end-1) + F2.v(1,2:end-1).*b.R.w(1,2:end-1).^ax).*b.JA.w(1,2:end-1)./b.DXI.rm(1,1:end-1).*b.DETA.zm(1,2:end-1);
        u_se(2:end-2,3:end-1)  = u_se(2:end-2,3:end-1) + u_SSe(1:end-3,3:end-1);
        u_nw(end-1,2:end-2)    = u_nw(end-1,2:end-2) + 0.5*(F1.v(end,2:end-1).*b.ETA_z.w(end,2:end-1) + F2.v(end,2:end-1).*b.R.w(end,2:end-1).^ax).*b.JA.w(end,2:end-1)./b.DXI.rm(end,1:end-1).*b.DETA.zm(end,2:end-1);
        u_nw(2:end-2,2:end-2)  = u_nw(2:end-2,2:end-2) + u_NNw(3:end-1,2:end-2);
        u_ne(end-1,3:end-1)    = u_ne(end-1,3:end-1) + 0.5*(F1.v(end,2:end-1).*b.ETA_z.w(end,2:end-1) + F2.v(end,2:end-1).*b.R.w(end,2:end-1).^ax).*b.JA.w(end,2:end-1)./b.DXI.rm(end,1:end-1).*b.DETA.zm(end,2:end-1);
        u_ne(2:end-2,3:end-1)  = u_ne(2:end-2,3:end-1) + u_NNe(3:end-1,3:end-1);
        
        u_SSw=u_SSw.'; u_SSw=u_SSw(:);
        u_SSe=u_SSe.'; u_SSe=u_SSe(:);
        u_sw=u_sw.'; u_sw=u_sw(:);
        u_se=u_se.'; u_se=u_se(:);
        u_nw=u_nw.'; u_nw=u_nw(:);
        u_ne=u_ne.'; u_ne=u_ne(:);
        u_NNw=u_NNw.'; u_NNw=u_NNw(:);
        u_NNe=u_NNe.'; u_NNe=u_NNe(:);
        
        b.z_m_visc.u=-spdiags([u_SSw u_SSe u_sw u_se u_nw u_ne u_NNw u_NNe],[-2*(b.J+2) -2*(b.J+2)+1 -b.J-2 -b.J-1 0 1 b.J+2 b.J+3],(b.K+2)*(b.J+2),(b.K+1)*(b.J+2));
        
        b.z_m_visc.u(b.J+2:b.J+2:end,:)=[];
        
    % viscous terms for w
        F1.v=mu(b.T.v).*b.JA.v.*(b.XI_r.v.^2+2*b.XI_z.v.^2).*b.R.v.^ax;
        F2.v=-mu(b.T.v).*b.JA.v.*b.XI_r.v.*b.XI_z.v.*b.R.v.^ax;
        F3.v=-2/3*mu(b.T.v).*b.XI_z.v;
        F4.v=2*mu(b.T.v).*b.JA.v.*b.XI_z.v.*b.ETA_z.v.*b.R.v.^ax;
        F5.v=-mu(b.T.v).*b.JA.v.*b.XI_r.v.*b.ETA_z.v.*b.R.v.^ax;
        
        F1.p=2*mu(b.T.p).*b.JA.p.*b.XI_z.p.*b.ETA_z.p.*b.R.p.^ax;
        F2.p=2*mu(b.T.p).*b.JA.p.*b.ETA_z.p.^2.*b.R.p.^ax;
        F3.p=-2/3*mu(b.T.p).*b.ETA_z.p;
        
        w_S=zeros(b.K+2,b.J+1); w_W=w_S; w_C=w_S; w_E=w_S; w_N=w_S; w_SW=w_S; w_SE=w_S; w_NW=w_S; w_NE=w_S;
        SW=0.25*ones(b.K,b.J-1);
        SW(1,:)=0.5;
        SE=-SW;
        NW=-flipud(SW);
        NE=-NW;
        
        w_S(1:end-2,2:end-1) = (F1.v(1:end-1,2:end-1).*b.XI_r.w(1:end-2,2:end-1) + F2.v(1:end-1,2:end-1).*b.XI_z.w(1:end-2,2:end-1)).*b.JA.w(1:end-2,2:end-1)./b.DXI.rm(1:end-1,1:end-1).*b.DETA.zm(:,2:end-1);
        w_W(2:end-1,1:end-2) = (F2.p(:,1:end-1).*b.XI_r.w(2:end-1,1:end-2) + F3.p(:,1:end-1).*b.R.w(2:end-1,1:end-2).^ax).*b.JA.w(2:end-1,1:end-2)./b.DETA.c(:,1:end-1).*b.DXI.zm(:,2:end-1);
        w_E(2:end-1,3:end)   = (F2.p(:,2:end).*b.XI_r.w(2:end-1,3:end) + F3.p(:,2:end).*b.R.w(2:end-1,3:end).^ax).*b.JA.w(2:end-1,3:end)./b.DETA.c(:,2:end).*b.DXI.zm(:,2:end-1);
        w_N(3:end,2:end-1)   = (F1.v(2:end,2:end-1).*b.XI_r.w(3:end,2:end-1) + F2.v(2:end,2:end-1).*b.XI_z.w(3:end,2:end-1)).*b.JA.w(3:end,2:end-1)./b.DXI.rm(2:end,1:end-1).*b.DETA.zm(:,2:end-1);
        w_C(2:end-1,2:end-1) = - ((F1.v(1:end-1,2:end-1)./b.DXI.rm(1:end-1,1:end-1) + F1.v(2:end,2:end-1)./b.DXI.rm(2:end,1:end-1)).*b.XI_r.w(2:end-1,2:end-1) + (F2.v(1:end-1,2:end-1)./b.DXI.rm(1:end-1,1:end-1) + F2.v(2:end,2:end-1)./b.DXI.rm(2:end,1:end-1)).*b.XI_z.w(2:end-1,2:end-1)).*b.JA.w(2:end-1,2:end-1).*b.DETA.zm(:,2:end-1);
        w_C(2:end-1,2:end-1) = w_C(2:end-1,2:end-1) - ((F2.p(:,1:end-1)./b.DETA.c(:,1:end-1) + F2.p(:,2:end)./b.DETA.c(:,2:end)).*b.XI_r.w(2:end-1,2:end-1) + (F3.p(:,1:end-1)./b.DETA.c(:,1:end-1) + F3.p(:,2:end)./b.DETA.c(:,2:end)).*b.R.w(2:end-1,2:end-1).^ax).*b.JA.w(2:end-1,2:end-1).*b.DXI.zm(:,2:end-1);
        
        w_SW(1:end-2,1:end-2) = SW.*(F3.v(1:end-1,2:end-1).*b.R.u(1:end-1,2:end-2).^ax + (F4.v(1:end-1,2:end-1) + F1.p(:,1:end-1)).*b.XI_r.u(1:end-1,2:end-2) + F5.v(1:end-1,2:end-1).*b.XI_z.u(1:end-1,2:end-2)).*b.JA.u(1:end-1,2:end-2);
        w_SE(1:end-2,3:end)   = SE.*(F3.v(1:end-1,2:end-1).*b.R.u(1:end-1,3:end-1).^ax + (F4.v(1:end-1,2:end-1) + F1.p(:,2:end)).*b.XI_r.u(1:end-1,3:end-1) + F5.v(1:end-1,2:end-1).*b.XI_z.u(1:end-1,3:end-1)).*b.JA.u(1:end-1,3:end-1);
        w_NW(3:end,1:end-2)   = NW.*(F3.v(2:end,2:end-1).*b.R.u(2:end,2:end-2).^ax + (F4.v(2:end,2:end-1) + F1.p(:,1:end-1)).*b.XI_r.u(2:end,2:end-2) + F5.v(2:end,2:end-1).*b.XI_z.u(2:end,2:end-2)).*b.JA.u(2:end,2:end-2);
        w_NE(3:end,3:end)     = NE.*(F3.v(2:end,2:end-1).*b.R.u(2:end,3:end-1).^ax + (F4.v(2:end,2:end-1) + F1.p(:,2:end)).*b.XI_r.u(2:end,3:end-1) + F5.v(2:end,2:end-1).*b.XI_z.u(2:end,3:end-1)).*b.JA.u(2:end,3:end-1);
        
        w_S(1:end-2,2:end-1) = w_S(1:end-2,2:end-1) + w_SW(1:end-2,1:end-2);
        w_S(1:end-2,2:end-1) = w_S(1:end-2,2:end-1) + w_SE(1:end-2,3:end);
        w_W(3:end-1,1:end-2) = w_W(3:end-1,1:end-2) + w_SW(2:end-2,1:end-2);
        w_E(3:end-1,3:end)   = w_E(3:end-1,3:end) + w_SE(2:end-2,3:end);
        w_C(3:end-1,2:end-1) = w_C(3:end-1,2:end-1) + w_SW(2:end-2,1:end-2);
        w_C(3:end-1,2:end-1) = w_C(3:end-1,2:end-1) + w_SE(2:end-2,3:end);
        w_N(3:end,2:end-1)   = w_N(3:end,2:end-1) + w_NW(3:end,1:end-2);
        w_N(3:end,2:end-1)   = w_N(3:end,2:end-1) + w_NE(3:end,3:end);
        w_W(2:end-2,1:end-2) = w_W(2:end-2,1:end-2) + w_NW(3:end-1,1:end-2);
        w_E(2:end-2,3:end)   = w_E(2:end-2,3:end) + w_NE(3:end-1,3:end);
        w_C(2:end-2,2:end-1) = w_C(2:end-2,2:end-1) + w_NW(3:end-1,1:end-2);
        w_C(2:end-2,2:end-1) = w_C(2:end-2,2:end-1) + w_NE(3:end-1,3:end);
        
        w_SW=w_SW.'; w_SW=w_SW(:);
        w_S=w_S.'; w_S=w_S(:);
        w_SE=w_SE.'; w_SE=w_SE(:);
        w_W=w_W.'; w_W=w_W(:);
        w_C=w_C.'; w_C=w_C(:);
        w_E=w_E.'; w_E=w_E(:);
        w_NW=w_NW.'; w_NW=w_NW(:);
        w_N=w_N.'; w_N=w_N(:);
        w_NE=w_NE.'; w_NE=w_NE(:);
        
        b.z_m_visc.w=-spdiags([w_SW w_S w_SE w_W w_C w_E w_NW w_N w_NE],[-b.J-2 -b.J-1 -b.J -1 0 1 b.J b.J+1 b.J+2],(b.K+2)*(b.J+1),(b.K+2)*(b.J+1));
                
    end
    
    b.z_m.u=b.z_m_visc.u;
    b.z_m.w=b.z_m_inert.w + b.z_m_visc.w;
    
    % pressure
    % -[J*xi_z*r^ax*p0]_delta.xi*delta.eta -[J*eta_z*r^ax*p0]_delta.eta*delta.xi
        F.p=-b.JA.p.*b.ETA_z.p.*b.R.p.^ax;
        F.v=-b.JA.v.*b.XI_z.v.*b.R.v.^ax;
        
        p_Sse=zeros(b.K,b.J+1); p_Ssw=p_Sse; p_e=p_Sse; p_w=p_Sse; p_Nne=p_Sse; p_Nnw=p_Sse;
    
        p_Ssw(1:end-1,1:end-2) = -b.SW.v(2:end,2:end).*F.v(2:end-1,2:end-1).*b.DETA.zm(2:end,2:end-1);
        p_Sse(1:end-1,2:end-1) = -b.SE.v(2:end,1:end-1).*F.v(2:end-1,2:end-1).*b.DETA.zm(2:end,2:end-1);
        p_Nnw(2:end,1:end-2)   = b.NW.v(1:end-1,2:end).*F.v(2:end-1,2:end-1).*b.DETA.zm(1:end-1,2:end-1);
        p_Nne(2:end,2:end-1)   = b.NE.v(1:end-1,1:end-1).*F.v(2:end-1,2:end-1).*b.DETA.zm(1:end-1,2:end-1);
        p_w(2:end,1:end-2)     = p_w(2:end,1:end-2) - b.NW.v(1:end-1,2:end).*F.v(2:end-1,2:end-1).*b.DETA.zm(2:end,2:end-1);
        p_e(2:end,2:end-1)     = p_e(2:end,2:end-1) - b.NE.v(1:end-1,1:end-1).*F.v(2:end-1,2:end-1).*b.DETA.zm(2:end,2:end-1);
        p_w(1:end-1,1:end-2)   = p_w(1:end-1,1:end-2) + b.SW.v(2:end,2:end).*F.v(2:end-1,2:end-1).*b.DETA.zm(1:end-1,2:end-1);
        p_e(1:end-1,2:end-1)   = p_e(1:end-1,2:end-1) + b.SE.v(2:end,1:end-1).*F.v(2:end-1,2:end-1).*b.DETA.zm(1:end-1,2:end-1);
        p_w(:,1:end-2)         = p_w(:,1:end-2) - F.p(:,1:end-1).*b.DXI.zm(:,2:end-1);
        p_e(:,2:end-1)         = p_e(:,2:end-1) + F.p(:,2:end).*b.DXI.zm(:,2:end-1);

        %{
        p_Ssw(end-1,1:end-2)   = p_Ssw(end-1,1:end-2) + (b.R_cyl.p(end,1:end-1)-b.R_cyl.T(end,2:end-2))./(b.R_cyl.p(end,1:end-1)-b.R_cyl.p(end-1,1:end-1)).*b.W.w(end-2,2:end).*F.v(end,2:end-1).*b.DETA.zm(end,2:end-1);
        p_Sse(end-1,2:end-1)   = p_Sse(end-1,2:end-1) + (b.R_cyl.p(end,2:end)-b.R_cyl.T(end,3:end-1))./(b.R_cyl.p(end,2:end)-b.R_cyl.p(end-1,2:end)).*b.E.w(end-2,1:end-1).*F.v(end,2:end-1).*b.DETA.zm(end,2:end-1);
        p_w(end,1:end-2)       = p_w(end,1:end-2) + (b.R_cyl.T(end,2:end-2)-b.R_cyl.p(end-1,1:end-1))./(b.R_cyl.p(end,1:end-1)-b.R_cyl.p(end-1,1:end-1)).*b.W.w(end-1,2:end).*F.v(end,2:end-1).*b.DETA.zm(end,2:end-1);
        p_e(end,2:end-1)       = p_e(end,2:end-1) + (b.R_cyl.T(end,3:end-1)-b.R_cyl.p(end-1,2:end))./(b.R_cyl.p(end,2:end)-b.R_cyl.p(end-1,2:end)).*b.E.w(end-1,1:end-1).*F.v(end,2:end-1).*b.DETA.zm(end,2:end-1);
        
        p_Nnw(2,1:end-2)       = p_Nnw(2,1:end-2) + (b.R_cyl.p(1,1:end-1)-b.R_cyl.T(1,2:end-2))./(b.R_cyl.p(1,1:end-1)-b.R_cyl.p(2,1:end-1)).*b.W.w(end-2,2:end).*F.v(1,2:end-1).*b.DETA.zm(1,2:end-1);
        p_Nne(2,2:end-1)       = p_Nne(2,2:end-1) + (b.R_cyl.p(1,2:end)-b.R_cyl.T(1,3:end-1))./(b.R_cyl.p(1,2:end)-b.R_cyl.p(2,2:end)).*b.E.w(end-2,1:end-1).*F.v(1,2:end-1).*b.DETA.zm(1,2:end-1);
        p_w(1,1:end-2)         = p_w(1,1:end-2) + (b.R_cyl.T(1,2:end-2)-b.R_cyl.p(2,1:end-1))./(b.R_cyl.p(1,1:end-1)-b.R_cyl.p(2,1:end-1)).*b.W.w(end-1,2:end).*F.v(1,2:end-1).*b.DETA.zm(1,2:end-1);
        p_e(1,2:end-1)         = p_e(1,2:end-1) + (b.R_cyl.T(1,3:end-1)-b.R_cyl.p(2,2:end))./(b.R_cyl.p(1,2:end)-b.R_cyl.p(2,2:end)).*b.E.w(end-1,1:end-1).*F.v(1,2:end-1).*b.DETA.zm(1,2:end-1);
        %}
        %%{
        p_Ssw(end-1,1:end-2)   = p_Ssw(end-1,1:end-2) + (-0.5).*b.W.w(end,2:end).*F.v(end,2:end-1).*b.DETA.zm(end,2:end-1);
        p_Sse(end-1,2:end-1)   = p_Sse(end-1,2:end-1) + (-0.5).*b.E.w(end,1:end-1).*F.v(end,2:end-1).*b.DETA.zm(end,2:end-1);
        p_w(end,1:end-2)       = p_w(end,1:end-2) + 1.5.*b.W.w(end,2:end).*F.v(end,2:end-1).*b.DETA.zm(end,2:end-1);
        p_e(end,2:end-1)       = p_e(end,2:end-1) + 1.5.*b.E.w(end,1:end-1).*F.v(end,2:end-1).*b.DETA.zm(end,2:end-1);
        
        p_Nnw(2,1:end-2)       = p_Nnw(2,1:end-2) - (-0.5).*b.W.w(1,2:end).*F.v(1,2:end-1).*b.DETA.zm(1,2:end-1);
        p_Nne(2,2:end-1)       = p_Nne(2,2:end-1) - (-0.5).*b.E.w(1,1:end-1).*F.v(1,2:end-1).*b.DETA.zm(1,2:end-1);
        p_w(1,1:end-2)         = p_w(1,1:end-2) - 1.5.*b.W.w(1,2:end).*F.v(1,2:end-1).*b.DETA.zm(1,2:end-1);
        p_e(1,2:end-1)         = p_e(1,2:end-1) - 1.5.*b.E.w(1,1:end-1).*F.v(1,2:end-1).*b.DETA.zm(1,2:end-1);
        %}
        
        p_Ssw=p_Ssw.'; p_Ssw=p_Ssw(:);
        p_Sse=p_Sse.'; p_Sse=p_Sse(:);
        p_w=p_w.'; p_w=p_w(:);
        p_e=p_e.'; p_e=p_e(:);
        p_Nnw=p_Nnw.'; p_Nnw=p_Nnw(:);
        p_Nne=p_Nne.'; p_Nne=p_Nne(:);

        b.z_m.p=-spdiags([p_Ssw p_Sse p_w p_e p_Nnw p_Nne],[-2*(b.J+1)-1 -2*(b.J+1) -(b.J+2) -(b.J+1) -1 0],(b.K+2)*(b.J+1),(b.K)*(b.J+1));
        b.z_m.p(:,b.J+1:b.J+1:end)=[];