function[b]=thermal_energy(b, flowopt)
    ax=flowopt.ax;
    rho=str2func(b.rho);
    lambda=str2func(b.lambda);
    cp=str2func(b.cp);
    
    % convective terms
    % [J*r^ax*rho0*cp0*u0*T0]_delta.xi*delta.eta + [J*r^ax*rho0*cp0*w0*dT0]_delta.eta*delta.xi
        F.u = b.JA.u.*b.R.u.^ax.*rho(b.T.u).*cp(b.T.u).*b.u.u;
        F.w = b.JA.w.*b.R.w.^ax.*rho(b.T.w).*cp(b.T.w).*b.w.w;
        
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
        
        b.th_e.T=spdiags([T_S T_W T_C T_E T_N],[-b.J-2 -1 0 1 b.J+2],(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
        
    % conductive terms
    % [J*r^ax*(xi_r^2+xi_z^2)*lambda0*dT0/dxi]_delta.xi*delta.eta + [J*r^ax*eta_z*xi_z*lambda0*dT0/deta]_delta.xi*delta.eta + [J*r^ax*eta_z*xi_z*lambda0*dT0/dxi]_delta.eta*delta.xi + [J*r^ax*eta_z^2*lambda0*dT0/deta]_delta.eta*delta.xi
        F1.u = b.JA.u.*b.R.u.^ax.*(b.XI_r.u.^2+b.XI_z.u.^2).*lambda(b.T.u);
        F2.u = b.JA.u.*b.R.u.^ax.*b.ETA_z.u.*b.XI_z.u.*lambda(b.T.u);
        F1.w = b.JA.w.*b.R.w.^ax.*b.ETA_z.w.*b.XI_z.w.*lambda(b.T.w);
        F2.w = b.JA.w.*b.R.w.^ax.*b.ETA_z.w.^2.*lambda(b.T.w);
        
        T_S=zeros(b.K+2,b.J+2); T_W=T_S; T_E=T_S; T_N=T_S; T_C=T_S; T_SW=T_S; T_SE=T_S; T_NW=T_S; T_NE=T_S;
        
        T_SW(1:end-2,1:end-2)= b.SW.v.*(F2.u(1:end-1,2:end-1)+F1.w(2:end-1,1:end-1));
        T_S(1:end-2,2:end-1) = F1.u(1:end-1,2:end-1)./b.DXI.rm(1:end-1,:).*b.DETA.c;
        T_SE(1:end-2,3:end)  = -b.SE.v.*(F2.u(1:end-1,2:end-1)+F1.w(2:end-1,2:end));
        T_W(2:end-1,1:end-2) = F2.w(2:end-1,1:end-1)./b.DETA.zm(:,1:end-1).*b.DXI.c;
        T_E(2:end-1,3:end)   = F2.w(2:end-1,2:end)./b.DETA.zm(:,2:end).*b.DXI.c;
        T_NW(3:end,1:end-2)  = -b.NW.v.*(F2.u(2:end,2:end-1)+F1.w(2:end-1,1:end-1));
        T_N(3:end,2:end-1)   =  F1.u(2:end,2:end-1)./b.DXI.rm(2:end,:).*b.DETA.c;
        T_NE(3:end,3:end)    = b.NE.v.*(F2.u(2:end,2:end-1)+F1.w(2:end-1,2:end));
        T_C(2:end-1,2:end-1) = -T_S(1:end-2,2:end-1) - T_N(3:end,2:end-1) - T_W(2:end-1,1:end-2) - T_E(2:end-1,3:end);
        
        T_S(1:end-2,3:end-1) = T_S(1:end-2,3:end-1) + b.SE.v(:,1:end-1)      .*(F2.u(1:end-1,3:end-1)+F1.w(2:end-1,2:end-1));
        T_S(1:end-2,2:end-2) = T_S(1:end-2,2:end-2) - b.SW.v(:,2:end)        .*(F2.u(1:end-1,2:end-2)+F1.w(2:end-1,2:end-1));
        T_W(3:end-1,1:end-2) = T_W(3:end-1,1:end-2) + b.NW.v(1:end-1,:)      .*(F2.u(2:end-1,2:end-1)+F1.w(3:end-1,1:end-1));
        T_E(3:end-1,3:end)   = T_E(3:end-1,3:end)   - b.NE.v(1:end-1,:)      .*(F2.u(2:end-1,2:end-1)+F1.w(3:end-1,2:end));
        T_C(3:end-1,3:end-1) = T_C(3:end-1,3:end-1) + b.NE.v(1:end-1,1:end-1).*(F2.u(2:end-1,3:end-1)+F1.w(3:end-1,2:end-1));
        T_C(3:end-1,2:end-2) = T_C(3:end-1,2:end-2) - b.NW.v(1:end-1,2:end)  .*(F2.u(2:end-1,2:end-2)+F1.w(3:end-1,2:end-1));
        T_N(3:end,3:end-1)   = T_N(3:end,3:end-1)   - b.NE.v(:,1:end-1)      .*(F2.u(2:end,3:end-1)  +F1.w(2:end-1,2:end-1));
        T_N(3:end,2:end-2)   = T_N(3:end,2:end-2)   + b.NW.v(:,2:end)        .*(F2.u(2:end,2:end-2)  +F1.w(2:end-1,2:end-1));
        T_W(2:end-2,1:end-2) = T_W(2:end-2,1:end-2) - b.SW.v(2:end,:)        .*(F2.u(2:end-1,2:end-1)+F1.w(2:end-2,1:end-1));
        T_E(2:end-2,3:end)   = T_E(2:end-2,3:end)   + b.SE.v(2:end,:)        .*(F2.u(2:end-1,2:end-1)+F1.w(2:end-2,2:end));
        T_C(2:end-2,3:end-1) = T_C(2:end-2,3:end-1) - b.SE.v(2:end,1:end-1)  .*(F2.u(2:end-1,3:end-1)+F1.w(2:end-2,2:end-1));
        T_C(2:end-2,2:end-2) = T_C(2:end-2,2:end-2) + b.SW.v(2:end,2:end)    .*(F2.u(2:end-1,2:end-2)+F1.w(2:end-2,2:end-1));
                
        T_SW=T_SW.'; T_SW=T_SW(:);
        T_S=T_S.'; T_S=T_S(:);
        T_SE=T_SE.'; T_SE=T_SE(:);
        T_W=T_W.'; T_W=T_W(:);
        T_C=T_C.'; T_C=T_C(:);
        T_E=T_E.'; T_E=T_E(:);
        T_NW=T_NW.'; T_NW=T_NW(:);
        T_N=T_N.'; T_N=T_N(:);
        T_NE=T_NE.'; T_NE=T_NE(:);
        
        b.th_e.T=b.th_e.T-spdiags([T_SW T_S T_SE T_W T_C T_E T_NW T_N T_NE],[-b.J-3 -b.J-2 -b.J-1 -1 0 1 b.J+1 b.J+2 b.J+3],(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));