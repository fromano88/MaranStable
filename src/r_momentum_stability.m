function[b]=r_momentum_stability(b, flowopt)
    ax=flowopt.ax;
    g=flowopt.g;
    m=flowopt.m;
    energy=flowopt.energy;
    inviscid=flowopt.inviscid;
    creeping=flowopt.creeping;
    iv=flowopt.iv;
    
    rho=str2func(b.rho);
    mu=str2func(b.mu);
    
    drho=str2func(b.drho);
    
    % gamma*J^2*r^ax*eta_z*rho0*u.hat*delta.xi*delta.eta
        F.u=b.JA.u.^2.*b.R.u.^ax.*rho(b.T.u).*b.ETA_z.u;

        u_C=zeros(b.K+1,b.J+2);

        u_C(2:end-1,2:end-1) = F.u(2:end-1,2:end-1).*b.DXI.rm(2:end-1,:).*b.DETA.rm(2:end-1,:);

        u_C=u_C.'; u_C=u_C(:);

        b.st_t.r_m.u=spdiags(u_C, 0,(b.K+1)*(b.J+2),(b.K+1)*(b.J+2));

    % gamma*J^2*r^ax*(-xi_z)*rho0*w.hat*delta.xi*delta.eta
        F.u=-b.JA.u.^2.*b.R.u.^ax.*rho(b.T.u).*b.XI_z.u;

        w_sw=zeros(b.K+1,b.J+2); w_se=w_sw; w_nw=w_sw; w_ne=w_sw;

        w_sw(2:end-1,2:end-1) = F.u(2:end-1,2:end-1)/4.*b.DXI.rm(2:end-1,:).*b.DETA.rm(2:end-1,:);
        w_se(2:end-1,3:end)   = w_sw(2:end-1,2:end-1);
        w_nw(3:end,2:end-1)   = w_sw(2:end-1,2:end-1);
        w_ne(3:end,3:end)     = w_sw(2:end-1,2:end-1);

        w_sw=w_sw.'; w_sw=w_sw(:);
        w_se=w_se.'; w_se=w_se(:);
        w_nw=w_nw.'; w_nw=w_nw(:);
        w_ne=w_ne.'; w_ne=w_ne(:);

        b.st_t.r_m.w=spdiags([w_sw w_se w_nw w_ne],[0 1 b.J+2 b.J+3],(b.K+1)*(b.J+2),(b.K+1)*(b.J+2));
        b.st_t.r_m.w(:,1:b.J+2:end)=[];
        b.st_t.r_m.w=[b.st_t.r_m.w sparse((b.K+1)*(b.J+2),b.J+1)];

    if energy==1
        % gamma*J^2*r^ax*(eta_z*u0-xi_z*w0)*drho0/dT0*T.hat*delta.xi*delta.eta
            F.u=b.JA.u.^2.*b.R.u.^ax.*drho(b.T.u).*(b.ETA_z.u.*b.u.u-b.XI_z.u.*b.w.u);

            T_s=sparse(b.K+1,b.J+2); T_n=T_s;

            T_s(2:end-1,2:end-1)=F.u(2:end-1,2:end-1)/2.*b.DXI.rm(2:end-1,:).*b.DETA.rm(2:end-1,:);
            T_n(2:end-1,2:end-1)=F.u(2:end-1,2:end-1)/2.*b.DXI.rm(2:end-1,:).*b.DETA.rm(2:end-1,:);

            T_s=T_s.'; T_s=T_s(:);
            T_n=T_n.'; T_n=T_n(:);

            b.st_t.r_m.T=spdiags([T_s T_n],[0 b.J+2],(b.K+1)*(b.J+2),(b.K+2)*(b.J+2));

    end
    
    if inviscid~=1
    % viscous terms for u
        F.u=-m^2*b.JA.u.^2.*b.ETA_z.u./b.R.u.^ax.*mu(b.T.u);
        
        u_C=zeros(b.K+1,b.J+2);
        
        u_C(2:end-1,2:end-1)=F.u(2:end-1,2:end-1).*b.DXI.rm(2:end-1,:).*b.DETA.rm(2:end-1,:);
        
        u_C=u_C.'; u_C=u_C(:);

        b.st.r_m.u=b.Ja.r_m.u - spdiags(u_C,0,(b.K+1)*(b.J+2),(b.K+1)*(b.J+2));
        
    % viscous terms for w
        F.u=m^2*b.JA.u.^2.*b.XI_z.u./b.R.u.^ax.*mu(b.T.u);
        
        w_sw=zeros(b.K+1,b.J+2); w_se=w_sw; w_nw=w_sw; w_ne=w_sw;

        w_sw(2:end-1,2:end-1) = F.u(2:end-1,2:end-1)/4.*b.DXI.rm(2:end-1,:).*b.DETA.rm(2:end-1,:);
        w_se(2:end-1,3:end)   = w_sw(2:end-1,2:end-1);
        w_nw(3:end,2:end-1)   = w_sw(2:end-1,2:end-1);
        w_ne(3:end,3:end)     = w_sw(2:end-1,2:end-1);

        w_sw=w_sw.'; w_sw=w_sw(:);
        w_se=w_se.'; w_se=w_se(:);
        w_nw=w_nw.'; w_nw=w_nw(:);
        w_ne=w_ne.'; w_ne=w_ne(:);

        b.st.r_m.w=-spdiags([w_sw w_se w_nw w_ne],[0 1 b.J+2 b.J+3],(b.K+1)*(b.J+2),(b.K+1)*(b.J+2));
        b.st.r_m.w(:,1:b.J+2:end)=[];
        b.st.r_m.w=[b.st.r_m.w sparse((b.K+1)*(b.J+2),b.J+1)];
        b.st.r_m.w=b.Ja.r_m.w + b.st.r_m.w;
            
    else
        b.st.r_m.u=b.Ja.r_m.u;
        b.st.r_m.w=b.Ja.r_m.w;
        
    end
    
    F.u=0; F1.u=0; F2.u=0; F1.T=0;
    
    if creeping~=1
    % m*J^2*rho0*(eta_z*u0-xi_z*w0)*delta.xi*delta.eta
        F.u=m*b.JA.u.^2.*rho(b.T.u).*(b.ETA_z.u.*b.u.u-b.XI_z.u.*b.w.u);
        
    else
        F.u=sparse(b.K+1,b.J+2);
        
    end
    
    if inviscid~=1
    % viscous terms for v
    % sign changed because of viscous terms!!!
        F1.u=-m*b.JA.u.*b.XI_r.u.*b.R.u.^ax.*mu(b.T.u);
        F2.u=ax*4/3*m*b.JA.u./b.R.u.^ax.*mu(b.T.u);

        F1.T=2/3*m*b.JA.T.*b.XI_r.T.*mu(b.T.T);

    else
        F1.u=sparse(b.K+1,b.J+2);
        F2.u=sparse(b.K+1,b.J+2);
        F1.T=sparse(b.K+2,b.J+2);
        
    end
    
        v_s=sparse(b.K+1,b.J+2); v_n=v_s;

        v_s(2:end-1,2:end-1)=((F.u(2:end-1,2:end-1)+F2.u(2:end-1,2:end-1))./2.*b.DXI.rm(2:end-1,:) - F1.T(2:end-2,2:end-1) - F1.u(2:end-1,2:end-1)./b.R.T(2:end-2,2:end-1).^ax).*b.DETA.rm(2:end-1,:);
        v_n(2:end-1,2:end-1)=((F.u(2:end-1,2:end-1)+F2.u(2:end-1,2:end-1))./2.*b.DXI.rm(2:end-1,:) + F1.T(3:end-1,2:end-1) + F1.u(2:end-1,2:end-1)./b.R.T(3:end-1,2:end-1).^ax).*b.DETA.rm(2:end-1,:);

        v_s=v_s.'; v_s=v_s(:);
        v_n=v_n.'; v_n=v_n(:);

        b.st.r_m.v=spdiags(iv*[v_s v_n],[0 b.J+2],(b.K+1)*(b.J+2),(b.K+2)*(b.J+2));
        