w_Ssw=zeros(K,b.J+1); w_Nnw=w_Ssw; w_w=w_Ssw; w_Cl=w_Ssw; w_E=w_Ssw; p_Cl=w_Ssw;
rhs1=zeros(K*(b.J+1),1);

if o==1
    % left boundary
        index_bc_r=0;
        J=1;

elseif o==2
    % right boundary
        index_bc_r=length(b.bc.r)/2;
        J=b.J+1;

end
% wall
    if strcmp(b.bc.r{index_bc_r+n}(1),'w')==1
        w_w(:,J)=1;

% inlet
    elseif strcmp(b.bc.r{index_bc_r+n}(1),'i')==1
        w_w(:,J)=1;

        rhs=str2func(b.bc.rhs.w.r{index_bc_r+n});
        rhs=rhs(b.R_cyl.w(b.R_cyl.w(:,1)>b.geom.r(n) & b.R_cyl.w(:,1)<b.geom.r(n+1),1)).*ones(K,1).*b.ETA_z.w(b.R_cyl.w(:,1)>b.geom.r(n) & b.R_cyl.w(:,1)<b.geom.r(n+1),J);
        rhs1(J:b.J+1:end,1)=rhs;

% outlet
    elseif strcmp(b.bc.r{index_bc_r+n}(1),'o')==1
        p_Cl(:,J)  = 1;

        rhs1(J:b.J+1:K*(b.J+1),1)=b.p_fluid;

% placeholder surface, interior, periodic
    else

    end

if o==1
    zm_w_Ssw = [zm_w_Ssw; w_Ssw];
    zm_w_w   = [zm_w_w; w_w];
    zm_w_Cl  = [zm_w_Cl; w_Cl];
    zm_w_Nnw = [zm_w_Nnw; w_Nnw];
    zm_w_E   = [zm_w_E; w_E];
    zm_p_Cl  = [zm_p_Cl; p_Cl];

    bc.rhs.wl  = [bc.rhs.wl; rhs1];

else
    zm_w_Sse = [zm_w_Sse; w_Ssw];
    zm_w_e   = [zm_w_e; w_w];
    zm_w_Cr  = [zm_w_Cr; w_Cl];
    zm_w_Nne = [zm_w_Nne; w_Nnw];
    zm_w_W   = [zm_w_W; w_E];
    zm_p_Cr  = [zm_p_Cr; p_Cl];

    bc.rhs.wr  = [bc.rhs.wr; rhs1];

end