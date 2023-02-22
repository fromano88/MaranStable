function[b]=phi_momentum_stability(b, flowopt)
    ax=flowopt.ax;
    m=flowopt.m;
    energy=flowopt.energy;
    inviscid=flowopt.inviscid;
    creeping=flowopt.creeping;
    iv=flowopt.iv;
    
    rho=str2func(b.rho);
    mu=str2func(b.mu);
    
    dmu=str2func(b.dmu);
    
    % gamma*J*r^ax*rho0*v.hat*delta.xi*delta.eta
        F.T=b.JA.T.*b.R.T.^ax.*rho(b.T.T);

        v_C=sparse(b.K+2,b.J+2);

        v_C(2:end-1,2:end-1)=F.T(2:end-1,2:end-1).*b.DXI.c.*b.DETA.c;

        v_C=v_C.'; v_C=v_C(:);

        b.st_t.phi_m.v=spdiags(iv*v_C,0,(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
        
    b.st.phi_m.v=sparse((b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
        
    if creeping~=1
    % inertia terms for v
    % [J*r^ax*rho0*u0*v.hat]_delta.xi*delta.eta + [J*r^ax*rho0*w0*v.hat]_delta.eta*delta.xi + J^2*rho0*(eta_z*u0-xi_z*w0)*v.hat*delta.xi*delta.eta
        F.u=b.JA.u.*b.R.u.^ax.*rho(b.T.u).*b.u.u;
        F.w=b.JA.w.*b.R.w.^ax.*rho(b.T.w).*b.w.w;
        F.T=ax*b.JA.T.^2.*(b.ETA_z.T.*b.u.T-b.XI_z.T.*b.w.T).*rho(b.T.T);

        v_S=zeros(b.K+2,b.J+2); v_W=v_S; v_E=v_S; v_N=v_S; v_C=v_S;

        v_S(1:end-2,2:end-1) = -b.S.u(:,2:end-1).*F.u(1:end-1,2:end-1).*b.DETA.c;
        v_W(2:end-1,1:end-2) = -b.W.w(2:end-1,:).*F.w(2:end-1,1:end-1).*b.DXI.c;
        v_C(2:end-1,2:end-1) = F.T(2:end-1,2:end-1).*b.DXI.c.*b.DETA.c;
        v_E(2:end-1,3:end)   = b.E.w(2:end-1,:).*F.w(2:end-1,2:end).*b.DXI.c;
        v_N(3:end,2:end-1)   = b.N.u(:,2:end-1).*F.u(2:end,2:end-1).*b.DETA.c;
        v_C(3:end-1,2:end-1) = v_C(3:end-1,2:end-1)-b.N.u(1:end-1,2:end-1).*F.u(2:end-1,2:end-1).*b.DETA.c(2:end,:);
        v_C(2:end-1,3:end-1) = v_C(2:end-1,3:end-1)-b.E.w(2:end-1,1:end-1).*F.w(2:end-1,2:end-1).*b.DXI.c(:,2:end);
        v_C(2:end-1,2:end-2) = v_C(2:end-1,2:end-2)+b.W.w(2:end-1,2:end).*F.w(2:end-1,2:end-1).*b.DXI.c(:,1:end-1);
        v_C(2:end-2,2:end-1) = v_C(2:end-2,2:end-1)+b.S.u(2:end,2:end-1).*F.u(2:end-1,2:end-1).*b.DETA.c(1:end-1,:);
        
        v_S=v_S.'; v_S=v_S(:);
        v_W=v_W.'; v_W=v_W(:);
        v_E=v_E.'; v_E=v_E(:);
        v_N=v_N.'; v_N=v_N(:);
        v_C=v_C.'; v_C=v_C(:);

        b.st.phi_m.v=spdiags(iv*[v_S v_W v_C v_E v_N],[-b.J-2 -1 0 1 b.J+2],(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
        
    end
    
    if inviscid~=1
    % viscous terms for u
        F.u=-m*b.JA.u.^2.*b.XI_r.u.*b.ETA_z.u.*b.R.u.^ax.*mu(b.T.u);
        
        F1.T=-ax*2*m*b.JA.T.^2.*b.ETA_z.T./b.R.T.^ax.*mu(b.T.T);
        F2.T=2/3*m./b.R.T.^ax.*mu(b.T.T);
        
        u_s=zeros(b.K+2,b.J+2); u_n=u_s;

        u_s(1:end-2,2:end-1) = (-1./b.R.T(2:end-1,2:end-1).^ax.*F.u(1:end-1,2:end-1) + F1.T(2:end-1,2:end-1)./2.*b.DXI.c - F2.T(2:end-1,2:end-1).*b.JA.u(1:end-1,2:end-1).*b.R.u(1:end-1,2:end-1).^ax).*b.DETA.c;
        u_n(2:end-1,2:end-1) = (1./b.R.T(2:end-1,2:end-1).^ax.*F.u(2:end,2:end-1) + F1.T(2:end-1,2:end-1)./2.*b.DXI.c + F2.T(2:end-1,2:end-1).*b.JA.u(2:end,2:end-1).*b.R.u(2:end,2:end-1).^ax).*b.DETA.c;

        u_s=u_s.'; u_s=u_s(:);
        u_n=u_n.'; u_n=u_n(:);

        b.st.phi_m.u=-spdiags([u_s u_n],[-b.J-2 0],(b.K+2)*(b.J+2),(b.K+1)*(b.J+2));
    
    % viscous terms for v
        %F.u=b.JA.u.*b.XI_r.u.^2.*(b.R.u.^ax).^2.*mu(b.T.u);
        F.u=b.JA.u.*b.XI_r.u.^2.*(b.R.u.^ax).^3.*mu(b.T.u);
        F1.u=b.JA.u.*b.XI_z.u.^2.*b.R.u.^ax.*mu(b.T.u);
        F2.u=b.JA.u.*b.XI_z.u.*b.ETA_z.u.*b.R.u.^ax.*mu(b.T.u);
        %F3.u=-b.JA.u.*b.XI_r.u.*b.R.u.^ax.*mu(b.T.u);
        
        F1.T=-4/3*m^2*b.JA.T./b.R.T.^ax.*mu(b.T.T);
        
        F1.w=b.JA.w.*b.XI_z.w.*b.ETA_z.w.*b.R.w.^ax.*mu(b.T.w);
        F2.w=b.JA.w.*b.ETA_z.w.^2.*b.R.w.^ax.*mu(b.T.w);
        
        v_S=zeros(b.K+2,b.J+2); v_W=v_S; v_E=v_S; v_N=v_S; v_C=v_S; v_SW=v_S; v_SE=v_S; v_NW=v_S; v_NE=v_S;
        invradius=1./b.R.T(:,2:end-1).^ax;
        invradius(isinf(invradius)==1)=0;
        v_SW(1:end-2,1:end-2)= b.SW.v.*(F2.u(1:end-1,2:end-1)+F1.w(2:end-1,1:end-1));
        %v_S(1:end-2,2:end-1) = ((1./b.R.T(2:end-1,2:end-1).^ax.*F.u(1:end-1,2:end-1) + F1.u(1:end-1,2:end-1))./b.DXI.rm(1:end-1,:) - 1./b.R.T(2:end-1,2:end-1).^ax.*b.S.u(:,2:end-1).*F3.u(1:end-1,2:end-1)).*b.DETA.c;
        v_S(1:end-2,2:end-1) = (1./b.R.T(2:end-1,2:end-1).^ax.*F.u(1:end-1,2:end-1).*invradius(1:end-2,:) + F1.u(1:end-1,2:end-1))./b.DXI.rm(1:end-1,:).*b.DETA.c;
        v_SE(1:end-2,3:end)  = -b.SE.v.*(F2.u(1:end-1,2:end-1)+F1.w(2:end-1,2:end));
        v_W(2:end-1,1:end-2) = F2.w(2:end-1,1:end-1)./b.DETA.zm(:,1:end-1).*b.DXI.c;
        v_E(2:end-1,3:end)   = F2.w(2:end-1,2:end)./b.DETA.zm(:,2:end).*b.DXI.c;
        v_NW(3:end,1:end-2)  = -b.NW.v.*(F2.u(2:end,2:end-1)+F1.w(2:end-1,1:end-1));
        %v_N(3:end,2:end-1)   = ((1./b.R.T(2:end-1,2:end-1).^ax.*F.u(2:end,2:end-1) + F1.u(2:end,2:end-1))./b.DXI.rm(2:end,:) + 1./b.R.T(2:end-1,2:end-1).^ax.*b.N.u(:,2:end-1).*F3.u(2:end,2:end-1)).*b.DETA.c;
        v_N(3:end,2:end-1)   = (1./b.R.T(2:end-1,2:end-1).^ax.*F.u(2:end,2:end-1).*invradius(3:end,:) + F1.u(2:end,2:end-1))./b.DXI.rm(2:end,:).*b.DETA.c;
        v_NE(3:end,3:end)    = b.NE.v.*(F2.u(2:end,2:end-1)+F1.w(2:end-1,2:end));
        v_C(2:end-1,2:end-1) = F1.T(2:end-1,2:end-1).*b.DXI.c.*b.DETA.c - (1./b.R.T(2:end-1,2:end-1).^ax.*F.u(1:end-1,2:end-1)./b.R.T(2:end-1,2:end-1).^ax + F1.u(1:end-1,2:end-1))./b.DXI.rm(1:end-1,:).*b.DETA.c - (1./b.R.T(2:end-1,2:end-1).^ax.*F.u(2:end,2:end-1)./b.R.T(2:end-1,2:end-1).^ax + F1.u(2:end,2:end-1))./b.DXI.rm(2:end,:).*b.DETA.c - v_W(2:end-1,1:end-2) - v_E(2:end-1,3:end);
        %v_C(2:end-1,2:end-1) = (F1.T(2:end-1,2:end-1).*b.DXI.c - (1./b.R.T(2:end-1,2:end-1).^ax.*F.u(1:end-1,2:end-1) + F1.u(1:end-1,2:end-1))./b.DXI.rm(1:end-1,:) - (1./b.R.T(2:end-1,2:end-1).^ax.*F.u(2:end,2:end-1)./b.R.T(2:end-1,2:end-1).^ax + F1.u(2:end,2:end-1))./b.DXI.rm(2:end,:)).*b.DETA.c - v_W(2:end-1,1:end-2) - v_E(2:end-1,3:end);
        %v_C(3:end-1,2:end-1) = v_C(3:end-1,2:end-1) - 1./b.R.T(3:end-1,2:end-1).^ax.*b.N.u(1:end-1,2:end-1).*F3.u(2:end-1,2:end-1).*b.DETA.c(2:end,:);
        %v_C(2:end-2,2:end-1) = v_C(2:end-2,2:end-1) + 1./b.R.T(2:end-2,2:end-1).^ax.*b.S.u(2:end,2:end-1).*F3.u(2:end-1,2:end-1).*b.DETA.c(1:end-1,:);
        
        v_S(1:end-2,3:end-1) = v_S(1:end-2,3:end-1) + b.SE.v(:,1:end-1).*(F2.u(1:end-1,3:end-1)+F1.w(2:end-1,2:end-1));
        v_S(1:end-2,2:end-2) = v_S(1:end-2,2:end-2) - b.SW.v(:,2:end).*(F2.u(1:end-1,2:end-2)+F1.w(2:end-1,2:end-1));
        v_W(3:end-1,1:end-2) = v_W(3:end-1,1:end-2) + b.NW.v(1:end-1,:).*(F2.u(2:end-1,2:end-1)+F1.w(3:end-1,1:end-1));
        v_E(3:end-1,3:end)   = v_E(3:end-1,3:end) - b.NE.v(1:end-1,:).*(F2.u(2:end-1,2:end-1)+F1.w(3:end-1,2:end));
        v_C(3:end-1,3:end-1) = v_C(3:end-1,3:end-1) + b.NE.v(1:end-1,1:end-1).*(F2.u(2:end-1,3:end-1)+F1.w(3:end-1,2:end-1));
        v_C(3:end-1,2:end-2) = v_C(3:end-1,2:end-2) - b.NW.v(1:end-1,2:end).*(F2.u(2:end-1,2:end-2)+F1.w(3:end-1,2:end-1));
        v_N(3:end,3:end-1)   = v_N(3:end,3:end-1) - b.NE.v(:,1:end-1).*(F2.u(2:end,3:end-1)+F1.w(2:end-1,2:end-1));
        v_N(3:end,2:end-2)   = v_N(3:end,2:end-2) + b.NW.v(:,2:end).*(F2.u(2:end,2:end-2)+F1.w(2:end-1,2:end-1));
        v_W(2:end-2,1:end-2) = v_W(2:end-2,1:end-2) - b.SW.v(2:end,:).*(F2.u(2:end-1,2:end-1)+F1.w(2:end-2,1:end-1));
        v_E(2:end-2,3:end)   = v_E(2:end-2,3:end) + b.SE.v(2:end,:).*(F2.u(2:end-1,2:end-1)+F1.w(2:end-2,2:end));
        v_C(2:end-2,3:end-1) = v_C(2:end-2,3:end-1) - b.SE.v(2:end,1:end-1).*(F2.u(2:end-1,3:end-1)+F1.w(2:end-2,2:end-1));
        v_C(2:end-2,2:end-2) = v_C(2:end-2,2:end-2) + b.SW.v(2:end,2:end).*(F2.u(2:end-1,2:end-2)+F1.w(2:end-2,2:end-1));
        
        v_SW=v_SW.'; v_SW=v_SW(:);
        v_S=v_S.'; v_S=v_S(:);
        v_SE=v_SE.'; v_SE=v_SE(:);
        v_W=v_W.'; v_W=v_W(:);
        v_E=v_E.'; v_E=v_E(:);
        v_NW=v_NW.'; v_NW=v_NW(:);
        v_N=v_N.'; v_N=v_N(:);
        v_NE=v_NE.'; v_NE=v_NE(:);
        v_C=v_C.'; v_C=v_C(:);

        b.st.phi_m.v=b.st.phi_m.v - spdiags(iv*[v_SW v_S v_SE v_W v_C v_E v_NW v_N v_NE],[-b.J-3 -b.J-2 -b.J-1 -1 0 1 b.J+1 b.J+2 b.J+3],(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
    
    % viscous terms for w
        F.u=m*b.JA.u.^2.*b.XI_r.u.*b.XI_z.u.*b.R.u.^ax.*mu(b.T.u);
        F1.u=-m*b.JA.u.^2.*b.XI_r.u.*b.XI_z.u.*mu(b.T.u);
        
        F1.T=ax*2*m*b.JA.T.^2.*b.XI_z.T./b.R.T.^ax.*mu(b.T.T);
        F2.T=2/3*m./b.R.T.^ax.*mu(b.T.T);
        
        F.w=-m*b.JA.w.^2.*b.XI_r.w.*b.ETA_z.w.*mu(b.T.w);
        
        w_w=zeros(b.K+2,b.J+2); w_e=w_w; w_Ssw=w_w; w_Sse=w_w; w_Nnw=w_w; w_Nne=w_w;
        
        w_Ssw(1:end-2,2:end-1) = -b.S.u(:,2:end-1)./2.*(1./b.R.T(2:end-1,2:end-1).^ax.*F.u(1:end-1,2:end-1) + F1.u(1:end-1,2:end-1)).*b.DETA.c;
        w_Sse(1:end-2,3:end)   = w_Ssw(1:end-2,2:end-1);
        w_Nnw(3:end,2:end-1)   = b.N.u(:,2:end-1)./2.*(1./b.R.T(2:end-1,2:end-1).^ax.*F.u(2:end,2:end-1) + F1.u(2:end,2:end-1)).*b.DETA.c;
        w_Nne(3:end,3:end)     = w_Nnw(3:end,2:end-1);
        w_w(3:end-1,2:end-1)   = w_Ssw(2:end-2,2:end-1);
        w_e(3:end-1,3:end)     = w_Ssw(2:end-2,2:end-1);
        w_w(2:end-2,2:end-1)   = w_w(2:end-2,2:end-1) + w_Nnw(3:end-1,2:end-1);
        w_e(2:end-2,3:end)     = w_e(2:end-2,3:end) + w_Nnw(3:end-1,2:end-1);
        
        w_w(2:end-1,2:end-1)   = w_w(2:end-1,2:end-1) + (F1.T(2:end-1,2:end-1)./2.*b.DETA.c - F2.T(2:end-1,2:end-1).*b.JA.w(2:end-1,1:end-1).*b.R.w(2:end-1,1:end-1).^ax - F.w(2:end-1,1:end-1)).*b.DXI.c;
        w_e(2:end-1,3:end)     = w_e(2:end-1,3:end) + (F1.T(2:end-1,2:end-1)./2.*b.DETA.c + F2.T(2:end-1,2:end-1).*b.JA.w(2:end-1,2:end).*b.R.w(2:end-1,2:end).^ax + F.w(2:end-1,2:end)).*b.DXI.c;
        
        w_Ssw=w_Ssw.'; w_Ssw=w_Ssw(:);
        w_Sse=w_Sse.'; w_Sse=w_Sse(:);
        w_w=w_w.'; w_w=w_w(:);
        w_e=w_e.'; w_e=w_e(:);
        w_Nnw=w_Nnw.'; w_Nnw=w_Nnw(:);
        w_Nne=w_Nne.'; w_Nne=w_Nne(:);
        
        b.st.phi_m.w=-spdiags([w_Ssw  w_Sse w_w w_e w_Nnw w_Nne],[-b.J-2 -b.J-1 0 1 b.J+2 b.J+3],(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
        b.st.phi_m.w(:,1:b.J+2:end)=[];
        
    % viscous terms for T
        if energy==1
            F1.T=-ax*2*m*b.JA.T.^2.*b.ETA_z.T./b.R.T.^ax.*b.u.T.*dmu(b.T.T);
            F2.T=2/3*m./b.R.T(2:end-1,2:end-1).^ax.*(b.JA.u(2:end,2:end-1).*b.R.u(2:end,2:end-1).^ax.*b.u.u(2:end,2:end-1) - b.JA.u(1:end-1,2:end-1).*b.R.u(1:end-1,2:end-1).^ax.*b.u.u(1:end-1,2:end-1))./b.DXI.c.*dmu(b.T.T(2:end-1,2:end-1));
            F3.T=ax*2*m*b.JA.T.^2.*b.XI_z.T./b.R.T.^ax.*b.w.T.*dmu(b.T.T);
            F4.T=2/3*m./b.R.T(2:end-1,2:end-1).^ax.*(b.JA.w(2:end-1,2:end).*b.R.w(2:end-1,2:end).^ax.*b.w.w(2:end-1,2:end) - b.JA.w(2:end-1,1:end-1).*b.R.w(2:end-1,1:end-1).^ax.*b.w.w(2:end-1,1:end-1))./b.DETA.c.*dmu(b.T.T(2:end-1,2:end-1));

            T_C=sparse(b.K+2,b.J+2);

            T_C(2:end-1,2:end-1)=(F1.T(2:end-1,2:end-1) + F2.T + F3.T(2:end-1,2:end-1) + F4.T).*b.DXI.c.*b.DETA.c;

            T_C=T_C.'; T_C=T_C(:);

            b.st.phi_m.T=-spdiags(T_C,0,(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));

        end
        
    else
        b.st.phi_m.u=sparse((b.K+2)*(b.J+2),(b.K+1)*(b.J+2));
        b.st.phi_m.w=sparse((b.K+2)*(b.J+2),(b.K+2)*(b.J+1));
        b.st.phi_m.T=sparse((b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
        
    end

    % m*J*p0*delta.xi*delta.eta
        p_C=sparse(b.K,b.J+2);

        p_C(:,2:end-1)=m*b.JA.p.*b.DXI.c.*b.DETA.c;

        p_C=p_C.'; p_C=p_C(:);

        b.st.phi_m.p=-spdiags(p_C,-(b.J+2),(b.K+2)*(b.J+2),(b.K)*(b.J+2));
        b.st.phi_m.p(:,1:b.J+2:end)=[];
        b.st.phi_m.p(:,b.J+1:b.J+1:end)=[];