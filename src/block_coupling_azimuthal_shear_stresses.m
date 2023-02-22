if o==1
    % lower boundary
        index_bc_z=0;
        index_bc_r=1;
        K=1;
        R_z.u=b.R1_z.u;

elseif o==2
    % upper boundary
        index_bc_z=length(b.bc.z)/2;
        index_bc_r=length(b.bc.r)/2;
        K=b.K+2;
        R_z.u=b.Rend_z.u;

end

if max(strncmp('s', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},1))
    N.u=sqrt(1+R_z.u(1,:).^2);

    % azimuthal shear stresses
        u_c=[0 -mu(b.T.T(K,2:end-1))./b.R.T(K,2:end-1).^ax.*m./b.XI_r.T(K,2:end-1).*b.dz.c 0];

        F.T=mu(b.T.T(K,:)).*b.XI_r.T(K,:)./[b.DXI.c(1,1) b.DXI.c(1,:) b.DXI.c(1,end)].*2;
        F.T=F.T(1,2:end-1).*b.dz.c;
        F1.T=-mu(b.T.T(K,:)).*R_z.u(1,:).*b.XI_z.T(K,:)./[b.DXI.c(1,1) b.DXI.c(1,:) b.DXI.c(1,end)].*2;
        F1.T=F1.T(1,2:end-1).*b.dz.c;
        F2.T=-mu(b.T.T(K,:)).*R_z.u(1,:).*b.ETA_z.T(K,:)./[b.DETA.zm(1,1) b.DETA.c(1,:) b.DETA.zm(1,end)];
        F2.T=F2.T(1,2:end-1).*b.dz.c;

        v_Wsw = [0 -F2.T/2 0];
        v_Ese = [0 F2.T/2 0];
        v_s   = (-1)^o*[0 F1.T+F.T 0];
        v_C   = -(-1)^o*[0 F1.T+F.T.*b.R.T(K,2:end-1).^ax./b.R.T(K-(-1)^o,2:end-1).^ax 0];
        v_Wsw(1,2) = -F2.T(1,1);
        v_s(1,2)   =  v_s(1,2)+F2.T(1,1)/2;
        v_Ese(1,end-1) = F2.T(1,end);
        v_s(1,end-1)   = v_s(1,end-1)-F2.T(1,end)/2;

    % thermocapillary azimuthal stresses
        if energy==1 && thermcapcon==1 && isfield(b, 'sigma')==1
            T_c=[0 (-1)^o*m./b.R.T(K,2:end-1).^ax.*dsigma(b.T.T(K,2:end-1)).*b.dz.c.*N.u(1,2:end-1) 0];

        else
            T_c=sparse(1,b.J+2);

        end

end
if o==1 && max(strncmp('s', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},1))
    st.ass1.u = spdiags(u_c.',0,(b.J+2),(b.K+1)*(b.J+2));
    st.ass1.v = spdiags(iv*[v_Wsw.' v_s.' v_Ese.' v_C.'],[-1 0 1 b.J+2],(b.J+2),(b.K+2)*(b.J+2));
    if energy==1
        st.ass1.T=spdiags(T_c.',0,(b.J+2),(b.K+2)*(b.J+2));

    end
    % equal velocity
        b.e1.v=speye(b.J+2,(b.K+2)*(b.J+2));

elseif o==2 && max(strncmp('s', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},1))
    st.assend.u = spdiags(u_c.',b.K*(b.J+2),(b.J+2),(b.K+1)*(b.J+2));
    st.assend.v = spdiags(iv*[v_C.' v_Wsw.' v_s.' v_Ese.'],(b.K+1)*(b.J+2)+[-(b.J+2) -1 0 1],(b.J+2),(b.K+2)*(b.J+2));
    if energy==1
        st.assend.T=spdiags(T_c.',(b.K+1)*(b.J+2),(b.J+2),(b.K+2)*(b.J+2));

    end
    % equal velocity
        b.eend.v=[sparse(b.J+2,(b.K+1)*(b.J+2)) speye(b.J+2)];

end