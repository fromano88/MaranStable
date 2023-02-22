function[b]=r_momentum(b, flowopt)
    ax=flowopt.ax;
    g=flowopt.g;
    inviscid=flowopt.inviscid;
    creeping=flowopt.creeping;
    rho=str2func(b.rho);
    mu=str2func(b.mu);
    
    b.r_m_inert.u=sparse((b.K+1)*(b.J+2),(b.K+1)*(b.J+2));
    b.r_m_inert.w=sparse((b.K+1)*(b.J+2),(b.K+2)*(b.J+1));
    b.r_m_visc.u=sparse((b.K+1)*(b.J+2),(b.K+1)*(b.J+2));
    b.r_m_visc.w=sparse((b.K+1)*(b.J+2),(b.K+2)*(b.J+1));
    
    if creeping~=1
    % inertia terms for u
    % [J^2*r^ax*rho0*(eta_z*u0^2-xi_z*u0*w0)]_delta.xi*delta.eta
        F.p=b.JA.p.^2.*b.R.p.^ax.*rho(b.T.p).*(b.ETA_z.p.*b.u.p-b.XI_z.p.*b.w.p);

        u_S=zeros(b.K+1,b.J+2); u_N=u_S; u_C=u_S;

        u_S(1:end-2,2:end-1) = -F.p(1:end-1,:)/2.*b.DETA.rm(2:end-1,:);
        u_N(3:end,2:end-1)   = F.p(2:end,:)/2.*b.DETA.rm(2:end-1,:);
        u_C(2:end-1,2:end-1) = u_C(2:end-1,2:end-1)+u_S(1:end-2,2:end-1);
        u_C(2:end-1,2:end-1) = u_C(2:end-1,2:end-1)+u_N(3:end,2:end-1);

        u_S=u_S.'; u_S=u_S(:);
        u_N=u_N.'; u_N=u_N(:);
        u_C=u_C.'; u_C=u_C(:);

        b.r_m_inert.u=spdiags([u_S u_C u_N],[-b.J-2 0 b.J+2],(b.K+1)*(b.J+2),(b.K+1)*(b.J+2));
        
    % inertia terms for w
    % [J^2*r^ax*rho0*(eta_z*u0*w0-xi_z*w0^2)]_delta.eta*delta.xi
        F.v=b.JA.v.^2.*b.R.v.^ax.*rho(b.T.v).*(b.ETA_z.v.*b.u.v-b.XI_z.v.*b.w.v);

        w_sw=zeros(b.K+1,b.J+2); w_se=w_sw; w_nw=w_sw; w_ne=w_sw;

        w_sw(2:end-1,2:end-1) = -F.v(2:end-1,1:end-1)/2.*b.DXI.rm(2:end-1,:);
        w_se(2:end-1,2:end-1) = F.v(2:end-1,2:end)/2.*b.DXI.rm(2:end-1,:);
        w_nw(2:end-1,2:end-1) = -F.v(2:end-1,1:end-1)/2.*b.DXI.rm(2:end-1,:);
        w_ne(2:end-1,2:end-1) = F.v(2:end-1,2:end)/2.*b.DXI.rm(2:end-1,:);

        w_sw=w_sw.'; w_sw=w_sw(:);
        w_se=w_se.'; w_se=w_se(:);
        w_nw=w_nw.'; w_nw=w_nw(:);
        w_ne=w_ne.'; w_ne=w_ne(:);

        b.r_m_inert.w=spdiags([w_sw w_se w_nw w_ne],[0 1 b.J+2 b.J+3],(b.K+1)*(b.J+2),(b.K+2)*(b.J+2));
        b.r_m_inert.w(:,1:b.J+2:end)=[];
           
    end
        
    if inviscid~=1
    % viscous terms for u
    
        F1.p=-2/3*mu(b.T.p).*b.XI_r.p;
        F2.p=mu(b.T.p).*b.JA.p.*(2*b.XI_r.p.^2+b.XI_z.p.^2).*b.R.p.^ax;
        F3.p=mu(b.T.p).*b.JA.p.*b.XI_z.p.*b.ETA_z.p.*b.R.p.^ax;

        F1.v=mu(b.T.v).*b.JA.v.*b.XI_z.v.*b.ETA_z.v.*b.R.v.^ax;
        F2.v=mu(b.T.v).*b.JA.v.*b.ETA_z.v.^2.*b.R.v.^ax;

        F1.u=-ax*2*mu(b.T.u).*b.JA.u.^2.*b.ETA_z.u./b.R.u.^ax;
        F2.u=ax*2/3*mu(b.T.u)./b.R.u.^ax;
        
        u_S=zeros(b.K+1,b.J+2); u_N=u_S; u_C=u_S; u_W=u_S; u_E=u_S; u_SW=u_S; u_SE=u_S; u_NW=u_S; u_NE=u_S;
        SW=0.25*ones(b.K-1,b.J);
        SW(:,1)=0.5;
        SE=-fliplr(SW);
        NW=-SW;
        NE=-SE;
        
        u_S(1:end-2,2:end-1) = (F1.p(1:end-1,:).*b.R.u(1:end-2,2:end-1).^ax + F2.p(1:end-1,:).*b.ETA_z.u(1:end-2,2:end-1)).*b.JA.u(1:end-2,2:end-1)./b.DXI.c(1:end-1,:).*b.DETA.rm(2:end-1,:);
        u_S(1:end-2,2:end-1) = u_S(1:end-2,2:end-1) - 0.5*F2.u(2:end-1,2:end-1).*b.JA.p(1:end-1,:).*b.R.p(1:end-1,:).^ax.*b.DETA.rm(2:end-1,:);
        u_W(2:end-1,1:end-2) = F2.v(2:end-1,1:end-1).*b.JA.u(2:end-1,1:end-2).*b.ETA_z.u(2:end-1,1:end-2)./b.DETA.zm(1:end-1,1:end-1).*b.DXI.rm(2:end-1,:);
        u_E(2:end-1,3:end)   = F2.v(2:end-1,2:end).*b.JA.u(2:end-1,3:end).*b.ETA_z.u(2:end-1,3:end)./b.DETA.zm(1:end-1,2:end).*b.DXI.rm(2:end-1,:);
        u_N(3:end,2:end-1)   = (F1.p(2:end,:).*b.R.u(3:end,2:end-1).^ax + F2.p(2:end,:).*b.ETA_z.u(3:end,2:end-1)).*b.JA.u(3:end,2:end-1)./b.DXI.c(2:end,:).*b.DETA.rm(2:end-1,:);
        u_N(3:end,2:end-1)   = u_N(3:end,2:end-1) + 0.5*F2.u(2:end-1,2:end-1).*b.JA.p(2:end,:).*b.R.p(2:end,:).^ax.*b.DETA.rm(2:end-1,:);
        u_C(2:end-1,2:end-1) = u_C(2:end-1,2:end-1) - ((F1.p(1:end-1,:)./b.DXI.c(1:end-1,:) + F1.p(2:end,:)./b.DXI.c(2:end,:)).*b.R.u(2:end-1,2:end-1).^ax + (F2.p(1:end-1,:)./b.DXI.c(1:end-1,:) + F2.p(2:end,:)./b.DXI.c(2:end,:)).*b.ETA_z.u(2:end-1,2:end-1)).*b.JA.u(2:end-1,2:end-1).*b.DETA.rm(2:end-1,:);
        u_C(2:end-1,2:end-1) = u_C(2:end-1,2:end-1) - (F2.v(2:end-1,1:end-1)./b.DETA.zm(1:end-1,1:end-1) + F2.v(2:end-1,2:end)./b.DETA.zm(1:end-1,2:end)).*b.JA.u(2:end-1,2:end-1).*b.ETA_z.u(2:end-1,2:end-1).*b.DXI.rm(2:end-1,:);
        u_C(2:end-1,2:end-1) = u_C(2:end-1,2:end-1) + (F1.u(2:end-1,2:end-1).*b.DXI.rm(2:end-1,:) + 0.5*F2.u(2:end-1,2:end-1).*(b.JA.p(2:end,:).*b.R.p(2:end,:).^ax - b.JA.p(1:end-1,:).*b.R.p(1:end-1,:).^ax)).*b.DETA.rm(2:end-1,:);
        
        u_SW(1:end-2,1:end-2) = SW.*(F3.p(1:end-1,:) + F1.v(2:end-1,1:end-1)).*b.JA.w(2:end-2,1:end-1).*b.ETA_z.w(2:end-2,1:end-1);
        u_SE(1:end-2,3:end)   = SE.*(F3.p(1:end-1,:) + F1.v(2:end-1,2:end)).*b.JA.w(2:end-2,2:end).*b.ETA_z.w(2:end-2,2:end);
        u_NW(3:end,1:end-2)   = NW.*(F3.p(2:end,:) + F1.v(2:end-1,1:end-1)).*b.JA.w(3:end-1,1:end-1).*b.ETA_z.w(3:end-1,1:end-1);
        u_NE(3:end,3:end)     = NE.*(F3.p(2:end,:) + F1.v(2:end-1,2:end)).*b.JA.w(3:end-1,2:end).*b.ETA_z.w(3:end-1,2:end);
        
        u_S(1:end-2,3:end-1) = u_S(1:end-2,3:end-1) + u_SW(1:end-2,2:end-2);
        u_S(1:end-2,2:end-2) = u_S(1:end-2,2:end-2) + u_SE(1:end-2,3:end-1);
        u_W(2:end-1,1:end-2) = u_W(2:end-1,1:end-2) + u_SW(1:end-2,1:end-2) + u_NW(3:end,1:end-2);
        u_C(2:end-1,3:end-1) = u_C(2:end-1,3:end-1) + u_SW(1:end-2,2:end-2) + u_NW(3:end,2:end-2);
        u_C(2:end-1,2:end-2) = u_C(2:end-1,2:end-2) + u_SE(1:end-2,3:end-1) + u_NE(3:end,3:end-1);
        u_E(2:end-1,3:end)   = u_E(2:end-1,3:end) + u_SE(1:end-2,3:end) + u_NE(3:end,3:end);
        u_N(3:end,3:end-1)   = u_N(3:end,3:end-1) + u_NW(3:end,2:end-2);
        u_N(3:end,2:end-2)   = u_N(3:end,2:end-2) + u_NE(3:end,3:end-1);
        
        u_SW=u_SW.'; u_SW=u_SW(:);
        u_S=u_S.'; u_S=u_S(:);
        u_SE=u_SE.'; u_SE=u_SE(:);
        u_W=u_W.'; u_W=u_W(:);
        u_C=u_C.'; u_C=u_C(:);
        u_E=u_E.'; u_E=u_E(:);
        u_NW=u_NW.'; u_NW=u_NW(:);
        u_N=u_N.'; u_N=u_N(:);
        u_NE=u_NE.'; u_NE=u_NE(:);
        
        b.r_m_visc.u=-spdiags([u_SW u_S u_SE u_W u_C u_E u_NW u_N u_NE],[-b.J-3 -b.J-2 -b.J-1 -1 0 1 b.J+1 b.J+2 b.J+3],(b.K+1)*(b.J+2),(b.K+1)*(b.J+2));
        
    % viscous terms for w
    
        F1.p=-2/3*mu(b.T.p).*b.XI_r.p;
        F2.p=-mu(b.T.p).*b.JA.p.*(2*b.XI_r.p.^2+b.XI_z.p.^2).*b.R.p.^ax;
        F3.p=-mu(b.T.p).*b.JA.p.*b.XI_z.p.*b.ETA_z.p.*b.R.p.^ax;
        F4.p=mu(b.T.p).*b.JA.p.*b.XI_r.p.*b.XI_z.p.*b.R.p.^ax;

        F1.v=-mu(b.T.v).*b.JA.v.*b.XI_z.v.*b.ETA_z.v.*b.R.v.^ax;
        F2.v=-mu(b.T.v).*b.JA.v.*b.ETA_z.v.^2.*b.R.v.^ax;
        F3.v=mu(b.T.v).*b.JA.v.*b.XI_r.v.*b.ETA_z.v.*b.R.v.^ax;

        F1.u=0.25*ax*2*mu(b.T.u(2:end-1,2:end-1)).*b.JA.u(2:end-1,2:end-1).^2.*b.XI_z.u(2:end-1,2:end-1)./b.R.u(2:end-1,2:end-1).^ax.*b.DXI.rm(2:end-1,:).*b.DETA.rm(2:end-1,:);
        F2.u=ax*2/3*mu(b.T.u(2:end-1,2:end-1))./b.R.u(2:end-1,2:end-1).^ax.*b.DXI.rm(2:end-1,:);
        
        w_sw=zeros(b.K+1,b.J+2); w_se=w_sw; w_nw=w_sw; w_ne=w_sw; w_Ssw=w_sw; w_Sse=w_sw; w_WsW=w_sw; w_EsE=w_sw; w_WnW=w_sw; w_EnE=w_sw; w_Nnw=w_sw; w_Nne=w_sw;

        C1 = - 0.25*((F2.p(1:end-1,:)./b.DXI.c(1:end-1,:) + F2.p(2:end,:)./b.DXI.c(2:end,:)).*b.XI_z.u(2:end-1,2:end-1) + (F4.p(1:end-1,:)./b.DXI.c(1:end-1,:) + F4.p(2:end,:)./b.DXI.c(2:end,:)).*b.XI_r.u(2:end-1,2:end-1)).*b.JA.u(2:end-1,2:end-1).*b.DETA.rm(2:end-1,:);
        C2 = - 0.25*(F2.v(2:end-1,1:end-1)./b.DETA.zm(1:end-1,1:end-1) + F2.v(2:end-1,2:end)./b.DETA.zm(1:end-1,2:end)).*b.JA.u(2:end-1,2:end-1).*b.XI_z.u(2:end-1,2:end-1).*b.DXI.rm(2:end-1,:);
        
        w_sw(2:end-1,2:end-1) = C1 + C2 + F1.u - 0.5*F2.u.*b.JA.v(2:end-1,1:end-1).*b.R.v(2:end-1,1:end-1).^ax + ((F3.p(1:end-1,:) + F1.v(2:end-1,1:end-1)).*b.XI_z.w(2:end-2,1:end-1) + F1.p(1:end-1,:).*b.R.w(2:end-2,1:end-1).^ax + F3.v(2:end-1,1:end-1).*b.XI_r.w(2:end-2,1:end-1)).*b.JA.w(2:end-2,1:end-1);
        w_se(2:end-1,2:end-1) = C1 + C2 + F1.u + 0.5*F2.u.*b.JA.v(2:end-1,2:end).*b.R.v(2:end-1,2:end).^ax - ((F3.p(1:end-1,:) + F1.v(2:end-1,2:end)).*b.XI_z.w(2:end-2,2:end) + F1.p(1:end-1,:).*b.R.w(2:end-2,2:end).^ax + F3.v(2:end-1,2:end).*b.XI_r.w(2:end-2,2:end)).*b.JA.w(2:end-2,2:end);
        w_nw(2:end-1,2:end-1) = C1 + C2 + F1.u - 0.5*F2.u.*b.JA.v(2:end-1,1:end-1).*b.R.v(2:end-1,1:end-1).^ax - ((F3.p(2:end,:) + F1.v(2:end-1,1:end-1)).*b.XI_z.w(3:end-1,1:end-1) + F1.p(2:end,:).*b.R.w(3:end-1,1:end-1).^ax + F3.v(2:end-1,1:end-1).*b.XI_r.w(3:end-1,1:end-1)).*b.JA.w(3:end-1,1:end-1);
        w_ne(2:end-1,2:end-1) = C1 + C2 + F1.u + 0.5*F2.u.*b.JA.v(2:end-1,2:end).*b.R.v(2:end-1,2:end).^ax + ((F3.p(2:end,:) + F1.v(2:end-1,2:end)).*b.XI_z.w(3:end-1,2:end) + F1.p(2:end,:).*b.R.w(3:end-1,2:end).^ax + F3.v(2:end-1,2:end).*b.XI_r.w(3:end-1,2:end)).*b.JA.w(3:end-1,2:end);
        
        S=0.25*ones(b.K-1,b.J);
        S(1,:)=0.5;
        N=flipud(S);
        
        w_Ssw(2:end-1,2:end-1) = S.*(F2.p(1:end-1,:).*b.XI_z.u(1:end-2,2:end-1) + F4.p(1:end-1,:).*b.XI_r.u(1:end-2,2:end-1)).*b.JA.u(1:end-2,2:end-1)./b.DXI.c(1:end-1,:).*b.DETA.rm(2:end-1,:);
        w_Sse                  = w_Ssw;
        w_Nnw(2:end-1,2:end-1) = N.*(F2.p(2:end,:).*b.XI_z.u(3:end,2:end-1) + F4.p(2:end,:).*b.XI_r.u(3:end,2:end-1)).*b.JA.u(3:end,2:end-1)./b.DXI.c(2:end,:).*b.DETA.rm(2:end-1,:);
        w_Nne                  = w_Nnw;
        w_sw(3:end-1,2:end-1)  = w_sw(3:end-1,2:end-1) + w_Ssw(3:end-1,2:end-1);
        w_se(3:end-1,2:end-1)  = w_se(3:end-1,2:end-1) + w_Sse(3:end-1,2:end-1);
        w_nw(2:end-2,2:end-1)  = w_nw(2:end-2,2:end-1) + w_Nnw(2:end-2,2:end-1);
        w_ne(2:end-2,2:end-1)  = w_ne(2:end-2,2:end-1) + w_Nne(2:end-2,2:end-1);
        
        w_WsW(2:end-1,3:end-1) = 0.25*F2.v(2:end-1,2:end-1).*b.JA.u(2:end-1,2:end-2).*b.XI_z.u(2:end-1,2:end-2)./b.DETA.zm(1:end-1,2:end-1).*b.DXI.rm(2:end-1,2:end);
        w_WnW                  = w_WsW;
        w_EsE(2:end-1,2:end-2) = 0.25*F2.v(2:end-1,2:end-1).*b.JA.u(2:end-1,3:end-1).*b.XI_z.u(2:end-1,3:end-1)./b.DETA.zm(1:end-1,2:end-1).*b.DXI.rm(2:end-1,1:end-1);
        w_EnE                  = w_EsE;
        w_sw(2:end-1,2)        = w_sw(2:end-1,2) + 0.5*F2.v(2:end-1,1).*b.JA.u(2:end-1,1).*b.XI_z.u(2:end-1,1)./b.DETA.zm(1:end-1,1).*b.DXI.rm(2:end-1,1);
        w_sw(2:end-1,3:end-1)  = w_sw(2:end-1,3:end-1) + w_WsW(2:end-1,3:end-1);
        w_se(2:end-1,end-1)    = w_se(2:end-1,end-1) + 0.5*F2.v(2:end-1,end).*b.JA.u(2:end-1,end).*b.XI_z.u(2:end-1,end)./b.DETA.zm(1:end-1,end).*b.DXI.rm(2:end-1,end);
        w_se(2:end-1,2:end-2)  = w_se(2:end-1,2:end-2) + w_EsE(2:end-1,2:end-2);
        w_nw(2:end-1,2)        = w_nw(2:end-1,2) + 0.5*F2.v(2:end-1,1).*b.JA.u(2:end-1,1).*b.XI_z.u(2:end-1,1)./b.DETA.zm(1:end-1,1).*b.DXI.rm(2:end-1,1);
        w_nw(2:end-1,3:end-1)  = w_nw(2:end-1,3:end-1) + w_WnW(2:end-1,3:end-1);
        w_ne(2:end-1,end-1)    = w_ne(2:end-1,end-1) + 0.5*F2.v(2:end-1,end).*b.JA.u(2:end-1,end).*b.XI_z.u(2:end-1,end)./b.DETA.zm(1:end-1,end).*b.DXI.rm(2:end-1,end);
        w_ne(2:end-1,2:end-2)  = w_ne(2:end-1,2:end-2) + w_EnE(2:end-1,2:end-2);
        
        w_Ssw=w_Ssw.'; w_Ssw=w_Ssw(:);
        w_Sse=w_Sse.'; w_Sse=w_Sse(:);
        w_WsW=w_WsW.'; w_WsW=w_WsW(:);
        w_sw=w_sw.'; w_sw=w_sw(:);
        w_se=w_se.'; w_se=w_se(:);
        w_EsE=w_EsE.'; w_EsE=w_EsE(:);
        w_WnW=w_WnW.'; w_WnW=w_WnW(:);
        w_nw=w_nw.'; w_nw=w_nw(:);
        w_ne=w_ne.'; w_ne=w_ne(:);
        w_EnE=w_EnE.'; w_EnE=w_EnE(:);
        w_Nnw=w_Nnw.'; w_Nnw=w_Nnw(:);
        w_Nne=w_Nne.'; w_Nne=w_Nne(:);
        
        b.r_m_visc.w=-spdiags([w_Ssw w_Sse w_WsW w_sw w_se w_EsE w_WnW w_nw w_ne w_EnE w_Nnw w_Nne],[-b.J-2 -b.J-1 -1 0 1 2 b.J+1 b.J+2 b.J+3 b.J+4 2*(b.J+2) 2*(b.J+2)+1],(b.K+1)*(b.J+2),(b.K+2)*(b.J+2));
        b.r_m_visc.w(:,1:b.J+2:end)=[];
        
    end
    
    b.r_m.u=b.r_m_inert.u + b.r_m_visc.u;
    b.r_m.w=b.r_m_inert.w + b.r_m_visc.w;
                
    % pressure
    % -d(J*xi_r*p0)/delta.xi*r^ax*delta.xi*delta.eta
        F.p=-b.JA.p.*b.XI_r.p;
        
        p_s=zeros(b.K,b.J+2); p_n=p_s;
        
        p_s(1:end-1,2:end-1)=-F.p(1:end-1,:).*b.R.u(2:end-1,2:end-1).^ax.*b.DETA.rm(2:end-1,:);
        p_n(2:end,2:end-1)=F.p(2:end,:).*b.R.u(2:end-1,2:end-1).^ax.*b.DETA.rm(2:end-1,:);
        
        p_s=p_s.'; p_s=p_s(:);
        p_n=p_n.'; p_n=p_n(:);
        
        b.r_m.p=-spdiags([p_s p_n],[-(b.J+2) 0],(b.K+1)*(b.J+2),(b.K)*(b.J+2));
        b.r_m.p(:,1:b.J+2:end)=[];
        b.r_m.p(:,b.J+1:b.J+1:end)=[];