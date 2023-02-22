function[b]=continuity_stability(b, flowopt)
    ax=flowopt.ax;
    m=flowopt.m;
    energy=flowopt.energy;
    iv=flowopt.iv;
    
    rho=str2func(b.rho);
    
    drho=str2func(b.drho);
    
    if energy==1
        % gamma*J*r^ax*drho0/dT0*T.hat*delta.xi*delta.eta
            T_C=sparse(b.K+2,b.J+2);

            T_C(2:end-1,2:end-1)=b.JA.p.*drho(b.T.p).*b.R.p.^ax.*b.DXI.c.*b.DETA.c;

            T_C=T_C.'; T_C=T_C(:);

            b.st_t.C.T=spdiags(T_C,0,(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
            b.st_t.C.T(1:b.J+2,:)=[];
            b.st_t.C.T(end-b.J-2:end,:)=[];
            b.st_t.C.T(1:b.J+2:end,:)=[];
            b.st_t.C.T(b.J+1:b.J+1:end,:)=[];
            
    end

    % m*J*rho0*v.hat*delta.xi*delta.eta
        v_C=sparse(b.K+2,b.J+2);

        v_C(2:end-1,2:end-1)=m*b.JA.p.*rho(b.T.p).*b.DXI.c.*b.DETA.c;

        v_C=v_C.'; v_C=v_C(:);

        b.st.C.v=spdiags(iv*v_C,0,(b.K+2)*(b.J+2),(b.K+2)*(b.J+2));
        b.st.C.v(1:b.J+2,:)=[];
        b.st.C.v(end-b.J-2:end,:)=[];
        b.st.C.v(1:b.J+2:end,:)=[];
        b.st.C.v(b.J+1:b.J+1:end,:)=[];
        
    if ((max(strncmp('o', b.bc.z(:),1)) || max(strncmp('o', b.bc.r(:),1)))~=1) && m==0
         ri=find(abs(b.r.w(2:end-1)-b.rp)==min(min(abs(b.r.w(2:end-1)-b.rp))));
         zi=find(abs(b.z.u(2:end-1)-b.zp)==min(min(abs(b.z.u(2:end-1)-b.zp))));
         b.st.C.v((ri(1)-1)*b.J+zi(1),:)=0;
         b.st_t.C.T((ri(1)-1)*b.J+zi(1),:)=0;
    end
    