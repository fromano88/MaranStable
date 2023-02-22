function [thermal_production, phi] = thermal_energy_production(b,gamma,t,m,Nphi)
    % this function renders the thermal energy production as 3D data

    if nargin == 4
        Nphi = 100;
    end

    phi  = 0:2*pi/Nphi:2*pi;

    rho = str2func(b.rho);
    cp  = str2func(b.cp);

    cpT_xi.u  = (cp(b.T.T(2:end,:)).*b.T.T(2:end,:)-cp(b.T.T(1:end-1,:)).*b.T.T(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm b.DXI.rm(:,end)];
    cpT_eta.w = (cp(b.T.T(:,2:end)).*b.T.T(:,2:end)-cp(b.T.T(:,1:end-1)).*b.T.T(:,1:end-1))./[b.DETA.zm(1,:); b.DETA.zm; b.DETA.zm(end,:)];
    cpT_xi.v  = interp2(b.ETA.u,b.XI.u,cpT_xi.u,b.ETA.v,b.XI.v,'spline');
    cpT_eta.v = interp2(b.ETA.w,b.XI.w,cpT_eta.w,b.ETA.v,b.XI.v,'spline');

    rho_cpT_xi_3D  = repmat(rho(b.T.v(:)).*cpT_xi.v(:),1,length(phi));
    rho_cpT_eta_3D = repmat(rho(b.T.v(:)).*cpT_eta.v(:),1,length(phi));
    thermal_production = rho_cpT_xi_3D.*plus_cc(b.T_hat.v(:),m,phi,gamma,t).*plus_cc(b.u_hat.v(:),m,phi,gamma,t);
    thermal_production = -thermal_production-rho_cpT_eta_3D.*plus_cc(b.T_hat.v(:),m,phi,gamma,t).*plus_cc(b.w_hat.v(:),m,phi,gamma,t);