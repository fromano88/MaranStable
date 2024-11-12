function[b]=thermal_energy_stability_V3d2(b, flowopt)
    ax=flowopt.ax;
    m=flowopt.m;
    iv=flowopt.iv;
    
    rho=str2func(b.rho);
    lambda=str2func(b.lambda);
    cp=str2func(b.cp);
    
    drho=str2func(b.drho);
    dcp=str2func(b.dcp);
    
    % gamma*J*r^ax*d(rho0*cp0*T0)/dT0*T.hat*delta.xi*delta.eta
        F.T=b.JA.T.*b.R.T.^ax.*(rho(b.T.T).*cp(b.T.T)+b.T.T.*drho(b.T.T).*cp(b.T.T));

        T_C=sparse(b.K+2,b.J+2);

        T_C(2:end-1,2:end-1)=F.T(2:end-1,2:end-1).*b.DXI.c.*b.DETA.c;

        T_C=T_C.'; T_C=T_C(:);

        b.st_t.th_e.T=spdiags(T_C,0,(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));

    % m*J*rho0*cp0*T0*v.hat*delta.xi*delta.eta
        F.T=m*b.JA.T.*rho(b.T.T).*cp(b.T.T).*b.T.T;

        v_C=sparse(b.K+2,b.J+2);

        v_C(2:end-1,2:end-1)=F.T(2:end-1,2:end-1).*b.DXI.c.*b.DETA.c;

        v_C=v_C.'; v_C=v_C(:);

        b.st.th_e.v=spdiags(iv*v_C,0,(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
        
    % conductive terms
        F.T=-m^2*b.JA.T./b.R.T.^ax.*lambda(b.T.T);

        T_C=sparse(b.K+2,b.J+2);

        T_C(2:end-1,2:end-1)=F.T(2:end-1,2:end-1).*b.DXI.c.*b.DETA.c;

        T_C=T_C.'; T_C=T_C(:);

        b.st.th_e.T=b.Ja.th_e.T-spdiags(T_C,0,(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
