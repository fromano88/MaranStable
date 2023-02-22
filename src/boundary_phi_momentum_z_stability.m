v_WW=zeros(1,J+2); v_Wsw=v_WW; v_s=v_WW; v_Ese=v_WW; v_EE=v_WW; v_C=v_WW; v_N=v_WW;

if o==1
    % lower boundary
        index_bc_z=0;
        index_bc_r=1;
        index_geom_r=1;
        K=1;
        N.u=N1.u;
        R_z.T=b.R1_z.T;

elseif o==2
    % upper boundary
        index_bc_z=length(b.bc.z)/2;
        index_bc_r=length(b.bc.r)/2;
        index_geom_r=length(b.geom.r);
        K=b.K+2;
        N.u=Nend.u;
        R_z.T=b.Rend_z.T;

end

% extrapolated
    if ((strcmp(b.bc.z{index_bc_z+n}(1),'w')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n}(2),'e')==1)
        v_s=ones(1,J+2);
        v_C=-3/2*ones(1,J+2);
        v_N=1/2*ones(1,J+2);

        if n==1
            if (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r}(2),'n')==1 || strcmp(b.bc.r{index_bc_r}(2),'m')==1)
                v_C(1,1)   = 0;
                v_N(1,1)   = 0;

            end
        end
        if n==length(b.bc.z)/2
            if (strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(2),'n')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(2),'m')==1)
                v_C(1,end)   = 0;
                v_N(1,end)   = 0;

            end
        end

% slip or axis
    elseif ((strcmp(b.bc.z{index_bc_z+n}(1),'w')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n}(2),'s')==1) || (strcmp(b.bc.z{index_bc_z+n}(1),'a')==1 && stability==1 && (m==1 || ax==0))
        F.T=1./N.u.*(b.XI_r.T(K,:)-b.R_z.T(K,:).*b.XI_z.T(K,:))./[b.DXI.rm(K-(o-1),1) b.DXI.rm(K-(o-1),:) b.DXI.rm(K-(o-1),end)];
        F.T=(-1)^o*F.T(1,b.Z.T(1,:)>b.geom.z(n) & b.Z.T(1,:)<b.geom.z(n+1));
        F1.T=-1./N.u.*b.R_z.T(K,:).*b.ETA_z.T(K,:)./[b.DETA.zm(1,1) b.DETA.c(1,:) b.DETA.zm(1,end)];
        F1.T=F1.T(1,b.Z.T(1,:)>b.geom.z(n) & b.Z.T(1,:)<b.geom.z(n+1));


        v_Wsw = [0 -F1.T/2 0];
        v_Ese = [0 F1.T/2 0];
        v_s   = [0 F.T 0];
        v_C   = [0 -F.T 0];

        if n==1
            v_Wsw(1,2) = -F1.T(1,1);
            v_s(1,2)   =  v_s(1,2)+F1.T(1,1)/2;
            if (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'s')==1
                v_s(1,1)   = 1;
                v_C(1,1)   = -1;

            elseif (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'e')==1
                v_s(1,1)   = 1;
                v_Ese(1,1) = -3/2;
                v_EE(1,1)  = 1/2;

            elseif (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r}(2),'n')==1 || strcmp(b.bc.r{index_bc_r}(2),'m')==1)
                v_s(1,1)   = 1;

            end
        end

        if n==length(b.bc.z)/2
            v_Ese(1,end-1) = F1.T(1,end);
            v_s(1,end-1)   = v_s(1,end-1)-F1.T(1,end)/2;
            if (strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(2),'s')==1
                v_s(1,end)   = 1;
                v_C(1,end)   = -1;

            elseif (strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(2),'e')==1
                v_s(1,end)   = 1;
                v_Wsw(1,end) = -3/2;
                v_WW(1,end)  = 1/2;

            elseif (strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(2),'n')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(2),'m')==1)
                v_s(1,end)   = 1;

            end
        end                        

% no slip, moving wall or axis
    elseif ((strcmp(b.bc.z{index_bc_z+n}(1),'w')==1 ||strcmp(b.bc.z{index_bc_z+n}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'o')==1) && (strcmp(b.bc.z{index_bc_z+n}(2),'n')==1 || strcmp(b.bc.z{index_bc_z+n}(2),'m')==1)) || (strcmp(b.bc.z{index_bc_z+n}(1),'a')==1 && (stability==1 && (m~=1 || ax==0)))
        v_s=ones(1,J+2);

% placeholder surface, interior, periodic
    else
        if n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r}(2),'n')==1 || strcmp(b.bc.r{index_bc_r}(2),'m')==1)
            v_s(1,1)=1;

        elseif n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'s')==1
            v_s(1,1)=-b.ETA_z.T(K,1)./b.DETA.zm(1,1)+(-1)^o*b.XI_z.T(K,1)./b.DXI.rm(K-(o-1),1);
            v_Ese(1,1)=b.ETA_z.T(K,1)./b.DETA.zm(1,1);
            v_C(1,1)=-(-1)^o*b.XI_z.T(K,1)./b.DXI.rm(K-(o-1),1);

        elseif n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'e')==1
            v_s(1,1)=1;
            v_Ese(1,1)=-3/2;
            v_EE(1,1)=1/2;

        end
        if n==length(b.bc.z)/2 && (strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(2),'n')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(2),'m')==1)
            v_s(1,end)=1;

        elseif n==length(b.bc.z)/2 && (strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(2),'s')==1
            v_s(1,end)=b.ETA_z.T(K,end)./b.DETA.zm(1,end)+(-1)^o*b.XI_z.T(K,end)./b.DXI.rm(K-(o-1),end);
            v_Wsw(1,end)=-b.ETA_z.T(K,end)./b.DETA.zm(1,end);
            v_C(1,end)=-(-1)^o*b.XI_z.T(K,end)./b.DXI.rm(K-(o-1),end);

        elseif n==length(b.bc.z)/2 && (strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+length(b.bc.r)/2}(2),'e')==1
            v_s(1,end)   = 1;
            v_Wsw(1,end) = -3/2;
            v_WW(1,end)  = 1/2;

        end

    end

if o==1
    st.phim.v = spdiags(iv*[v_WW.' v_Wsw.' v_s.' v_Ese.' v_EE.' v_C.' v_N.'],length(b.Z.T(b.Z.T(1,:)<=b.geom.z(n)))-1+[-2 -1 0 1 2 b.J+2 2*(b.J+2)],(J+2),(b.K+2)*(b.J+2));

else
    st.phim.v = spdiags(iv*[v_N.' v_C.' v_WW.' v_Wsw.' v_s.' v_Ese.' v_EE.'],(b.K+2)*(b.J+2)-length(b.Z.T(b.Z.T(1,:)>b.geom.z(n)))-1+[-2*(b.J+2) -(b.J+2) -2 -1 0 1 2],(J+2),(b.K+2)*(b.J+2));

end

if n~=1
   st.phim.v(1,:)=[];

end
if n~=length(b.bc.z)/2
   st.phim.v(end,:)=[];

end

if o==1
    bc.st.phim.v1 = [bc.st.phim.v1; st.phim.v];

else
    bc.st.phim.vend = [bc.st.phim.vend; st.phim.v];

end