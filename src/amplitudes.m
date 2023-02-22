function b = amplitudes(eigenvector,vector_index,flowopt,b)

ax = flowopt.ax;

b.u_hat.u = reshape(eigenvector(vector_index+1:vector_index+(b.K+1)*(b.J+2),1),[b.J+2,b.K+1]); b.u_hat.u = b.u_hat.u.';
b.u_hat.v = interp2(b.ETA.u,b.XI.u,b.u_hat.u,b.ETA.v,b.XI.v,'spline');
b.u_hat.T = interp2(b.ETA.u,b.XI.u,b.u_hat.u,b.ETA.T,b.XI.T,'spline');
b.u_hat.p = b.u_hat.T(2:end-1,2:end-1);
b.v_hat.T = reshape(eigenvector(vector_index+(b.K+1)*(b.J+2)+1:vector_index+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2),1),[b.J+2,b.K+2]); b.v_hat.T = b.v_hat.T.'*flowopt.iv/1i;
b.v_hat.p = b.v_hat.T(2:end-1,2:end-1);
b.v_hat.v = interp2(b.ETA.T,b.XI.T,b.v_hat.T,b.ETA.v,b.XI.v,'spline');
b.w_hat.w = reshape(eigenvector(vector_index+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+1:vector_index+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1),1),[b.J+1,b.K+2]); b.w_hat.w = b.w_hat.w.';
b.w_hat.v = interp2(b.ETA.w,b.XI.w,b.w_hat.w,b.ETA.v,b.XI.v,'spline');
b.w_hat.T = interp2(b.ETA.w,b.XI.w,b.w_hat.w,b.ETA.T,b.XI.T,'spline');
b.w_hat.p = b.w_hat.T(2:end-1,2:end-1);
b.p_hat.p = reshape(eigenvector(vector_index+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)+1:vector_index+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J,1),[b.J,b.K]); b.p_hat.p = b.p_hat.p.';
b.p_hat.v = interp1(b.xi.w(2:end-1).',b.p_hat.p,b.xi.u.','spline','extrap');
b.p_hat.v = interp1(b.eta.u(2:end-1).',b.p_hat.v.',b.eta.w.','spline','extrap').';

if flowopt.energy==1
    b.T_hat.T = reshape(eigenvector(vector_index+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J+1:vector_index+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J+(b.K+2)*(b.J+2),1),[b.J+2,b.K+2]); b.T_hat.T = b.T_hat.T.';
    b.T_hat.p = b.T_hat.T(2:end-1,2:end-1);
    b.T_hat.v = interp2(b.ETA.T,b.XI.T,b.T_hat.T,b.ETA.v,b.XI.v,'spline');
    b.T_hat.u = interp2(b.ETA.T,b.XI.T,b.T_hat.T,b.ETA.u,b.XI.u,'spline');
    b.T_hat.w = interp2(b.ETA.T,b.XI.T,b.T_hat.T,b.ETA.w,b.XI.w,'spline');
    
    rho = str2func(b.rho);
    cp  = str2func(b.cp);
    b.conv_energy_u_hat = rho(b.T.u).*cp(b.T.u).*b.R.u.^ax.*b.JA.u.*b.T.u.*b.u_hat.u;
    b.conv_energy_w_hat = rho(b.T.w).*cp(b.T.w).*b.R.w.^ax.*b.JA.w.*b.T.w.*b.w_hat.w;
    b.conv_energy_hat = 1./b.JA.p./b.R.p.^ax.*((b.conv_energy_u_hat(2:end,2:end-1)-b.conv_energy_u_hat(1:end-1,2:end-1))./b.DXI.c+(b.conv_energy_w_hat(2:end-1,2:end)-b.conv_energy_w_hat(2:end-1,1:end-1))./b.DETA.c);

end