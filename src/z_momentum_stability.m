function[b]=z_momentum_stability(b, flowopt)
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
    
    % gamma*J^2*r^ax*xi_r*rho0*w.hat*delta.xi*delta.eta
        F.w=b.JA.w.^2.*b.R.w.^ax.*rho(b.T.w).*b.XI_r.w;

        w_C=zeros(b.K+2,b.J+1);

        w_C(2:end-1,2:end-1) = F.w(2:end-1,2:end-1).*b.DXI.zm(:,2:end-1).*b.DETA.zm(:,2:end-1);

        w_C=w_C.'; w_C=w_C(:);

        b.st_t.z_m.w=spdiags(w_C, 0,(b.K+2)*(b.J+1),(b.K+2)*(b.J+1));

    if energy==1
        % gamma*J^2*r^ax*xi_r*w0)*drho0/dT0*T.hat*delta.xi*delta.eta
            F.w=b.JA.w.^2.*b.R.w.^ax.*drho(b.T.w).*b.XI_r.w.*b.w.w;

            T_w=sparse(b.K+2,b.J+2); T_e=T_w;

            T_w(2:end-1,2:end-2)=F.w(2:end-1,2:end-1)/2.*b.DXI.zm(:,2:end-1).*b.DETA.zm(:,2:end-1);
            T_e(2:end-1,3:end-1)=F.w(2:end-1,2:end-1)/2.*b.DXI.zm(:,2:end-1).*b.DETA.zm(:,2:end-1);

            T_w=T_w.'; T_w=T_w(:);
            T_e=T_e.'; T_e=T_e(:);

            b.st_t.z_m.T=spdiags([T_w T_e],[0 1],(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
            b.st_t.z_m.T(b.J+2:b.J+2:end,:)=[];

    end
    
    F.w=0; F1.w=0; F.v=0;
    
    if creeping~=1
    % m*J^2*rho0*xi_r*w0*delta.xi*delta.eta
        F.w=m*b.JA.w.^2.*rho(b.T.w).*b.XI_r.w.*b.w.w;
        
    else
        F.w=sparse(b.K+2,b.J+1);
        
    end
    
    if inviscid~=1
    % viscous terms for v
    % sign changed because of viscous terms!!!
        F1.w=-m*b.JA.w.*b.ETA_z.w.*mu(b.T.w);
        F2.w=-m*b.JA.w.*b.XI_z.w.*mu(b.T.w);

        F.v=2/3*m*b.JA.v.*b.XI_z.v.*mu(b.T.v);
        F.p=2/3*m*b.JA.p.*b.ETA_z.p.*mu(b.T.p);
        
    else
        F1.w=sparse(b.K+2,b.J+1);
        F2.w=sparse(b.K+2,b.J+1);
        F.v=sparse(b.K+1,b.J+1);
        F.p=sparse(b.K,b.J);
                
    end
        
        v_w=sparse(b.K+2,b.J+2); v_e=v_w; v_Ssw=v_w; v_Sse=v_w; v_Nnw=v_w; v_Nne=v_w;

        v_Ssw(2:end-2,2:end-2)=-0.25*(F2.w(3:end-1,2:end-1) + F.v(2:end-1,2:end-1)).*b.DETA.zm(2:end,2:end-1);
        v_Sse(2:end-2,3:end-1)=-0.25*(F2.w(3:end-1,2:end-1) + F.v(2:end-1,2:end-1)).*b.DETA.zm(2:end,2:end-1);
        v_Nnw(3:end-1,2:end-2)=0.25*(F2.w(2:end-2,2:end-1) + F.v(2:end-1,2:end-1)).*b.DETA.zm(1:end-1,2:end-1);
        v_Nne(3:end-1,3:end-1)=0.25*(F2.w(2:end-2,2:end-1) + F.v(2:end-1,2:end-1)).*b.DETA.zm(1:end-1,2:end-1);
        v_w(2:end-1,2:end-2)=(F.w(2:end-1,2:end-1)./2.*b.DETA.zm(:,2:end-1) - (F1.w(2:end-1,2:end-1) + F.p(:,1:end-1))).*b.DXI.zm(:,2:end-1) + v_Ssw(1:end-2,2:end-2) + v_Nnw(3:end,2:end-2);
        v_e(2:end-1,3:end-1)=(F.w(2:end-1,2:end-1)./2.*b.DETA.zm(:,2:end-1) + (F1.w(2:end-1,2:end-1) + F.p(:,2:end))).*b.DXI.zm(:,2:end-1) + v_Sse(1:end-2,3:end-1) + v_Nne(3:end,3:end-1);
        
        v_Ssw(1,2:end-2)=-0.5*(F2.w(2,2:end-1) + F.v(1,2:end-1)).*b.DETA.zm(1,2:end-1);
        v_Sse(1,3:end-1)=-0.5*(F2.w(2,2:end-1) + F.v(1,2:end-1)).*b.DETA.zm(1,2:end-1);
        v_Nnw(end,2:end-2)=0.5*(F2.w(end-1,2:end-1) + F.v(end,2:end-1)).*b.DETA.zm(end,2:end-1);
        v_Nne(end,3:end-1)=0.5*(F2.w(end-1,2:end-1) + F.v(end,2:end-1)).*b.DETA.zm(end,2:end-1);
        
        v_Ssw=v_Ssw.'; v_Ssw=v_Ssw(:);
        v_Sse=v_Sse.'; v_Sse=v_Sse(:);
        v_w=v_w.'; v_w=v_w(:);
        v_e=v_e.'; v_e=v_e(:);
        v_Nnw=v_Nnw.'; v_Nnw=v_Nnw(:);
        v_Nne=v_Nne.'; v_Nne=v_Nne(:);

        b.st.z_m.v=spdiags(iv*[v_Ssw v_Sse v_w v_e v_Nnw v_Nne],[-b.J-2 -b.J-1 0 1 b.J+2 b.J+3],(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
        b.st.z_m.v(b.J+2:b.J+2:end,:)=[];
        
    if inviscid~=1
    % viscous terms for w
        F.w=-m^2*b.JA.w.^2.*b.XI_r.w./b.R.w.^ax.*mu(b.T.w);
        
        w_C=zeros(b.K+2,b.J+1);
        
        w_C(2:end-1,2:end-1)=F.w(2:end-1,2:end-1).*b.DXI.zm(:,2:end-1).*b.DETA.zm(:,2:end-1);
        
        w_C=w_C.'; w_C=w_C(:);

        b.st.z_m.w=b.Ja.z_m.w - spdiags(w_C,0,(b.K+2)*(b.J+1),(b.K+2)*(b.J+1));
    
    end