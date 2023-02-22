v_Ssw=zeros(K,b.J+2); v_Nnw=v_Ssw; v_w=v_Ssw; v_Cl=v_Ssw; v_E=v_Ssw; v_Sse=v_Ssw; v_Nne=v_Ssw; v_e=v_Ssw; v_Cr=v_Ssw; v_W=v_Ssw;

if o==1
    % left boundary
        index_bc_r=0;
        J=1;

elseif o==2
    % right boundary
        index_bc_r=length(b.bc.r)/2;
        J=b.J+2;

end

% extrapolated
    if ((strcmp(b.bc.r{index_bc_r+n}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n}(2),'e')==1)
        v_w(:,J)   = ones(K,1);
        v_Cl(:,J)  = -3/2*ones(K,1);
        v_E(:,J)   = 1/2*ones(K,1);


% slip
    elseif ((strcmp(b.bc.r{index_bc_r+n}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n}(2),'s')==1)
        F.T=b.ETA_z.T(:,J)./[b.DETA.zm(1,J-1*(o-1)); b.DETA.zm(:,J-1*(o-1)); b.DETA.zm(end,J-1*(o-1))];
        F.T=(-1)^o*F.T(b.R_cyl.T(:,1)>b.geom.r(n) & b.R_cyl.T(:,1)<b.geom.r(n+1),1);
        F1.T=b.XI_z.T(:,J)./[b.DXI.rm(1,1); b.DXI.c(:,1); b.DXI.rm(end,1)];
        F1.T=F1.T(b.R_cyl.T(:,1)>b.geom.r(n) & b.R_cyl.T(:,1)<b.geom.r(n+1),1);

        v_Ssw(:,J) = -F1.T/2;
        v_Nnw(:,J) = F1.T/2;
        v_w(:,J)   = F.T;
        v_Cl(:,J)  = -F.T;

        if n==1
            v_Ssw(1,J) = -F1.T(1,1);
            v_w(1,J)   = v_w(1,1)+F1.T(1,1)/2;

        end

        if n==length(b.bc.r)/2
            v_Nnw(end,J) = F1.T(end,1);
            v_w(end,J)   = v_w(end,1)-F1.T(end,1)/2;

        end

% no slip
    elseif (strcmp(b.bc.r{index_bc_r+n}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r+n}(2),'n')==1 || strcmp(b.bc.r{index_bc_r+n}(2),'m')==1)
        v_w(:,J)=ones(K,1);

% placeholder surface, interior, periodic
    else

    end

if o==1
    phim_v_Ssw = [phim_v_Ssw; v_Ssw];
    phim_v_w   = [phim_v_w; v_w];
    phim_v_Cl  = [phim_v_Cl; v_Cl];
    phim_v_Nnw = [phim_v_Nnw; v_Nnw];
    phim_v_E   = [phim_v_E; v_E];

else
    phim_v_Sse = [phim_v_Sse; v_Ssw];
    phim_v_e   = [phim_v_e; v_w];
    phim_v_Cr  = [phim_v_Cr; v_Cl];
    phim_v_Nne = [phim_v_Nne; v_Nnw];
    phim_v_W   = [phim_v_W; v_E];

end