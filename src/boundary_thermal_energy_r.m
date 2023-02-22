T_Ssw=zeros(K,b.J+2); T_Nnw=T_Ssw; T_w=T_Ssw; T_Cl=T_Ssw; T_E=T_Ssw; T_Sse=T_Ssw; T_Nne=T_Ssw; T_e=T_Ssw; T_Cr=T_Ssw; T_W=T_Ssw;
rhs1=zeros(K*(b.J+2),1);

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
    if ((strcmp(b.bc.r{index_bc_r+n}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n}(3),'e')==1)
        T_w(:,J)   = ones(K,1);
        T_Cl(:,J)  = -1.5;%-(b.Z.T(1:K,J)-b.Z.T(1:K,J-2*(-1)^o))./(b.Z.T(1:K,J-1*(-1)^o)-b.Z.T(1:K,J-2*(-1)^o));
        T_E(:,J)   = 0.5;%-(b.Z.T(1:K,J-1*(-1)^o)-b.Z.T(1:K,J))./(b.Z.T(1:K,J-1*(-1)^o)-b.Z.T(1:K,J-2*(-1)^o));

% adiabatic
    elseif ((strcmp(b.bc.r{index_bc_r+n}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n}(3),'a')==1)
        F.T=b.ETA_z.T(:,J)./[b.DETA.zm(1,J-1*(o-1)); b.DETA.zm(:,J-1*(o-1)); b.DETA.zm(end,J-1*(o-1))].*[b.dr.rm(1) b.dr.c b.dr.rm(end)].';
        F.T=(-1)^o*F.T(b.R_cyl.T(:,1)>b.geom.r(n) & b.R_cyl.T(:,1)<b.geom.r(n+1),1);

        F1.T=b.XI_z.T(:,J)./[b.DXI.rm(1,1); b.DXI.c(:,1); b.DXI.rm(end,1)];
        F1.T=F1.T(2:end-1,:).*b.dr.c.';
        F2.T=[0; -b.S.u(:,1).*F1.T; 0];
        F3.T=[0; b.N.u(:,1).*F1.T; 0];
        F4.T=[0; -[0; b.N.u(1:end-1,1)].*F1.T+[b.S.u(2:end,1); 0].*F1.T; 0];

        F2.T=F2.T(b.R_cyl.T(:,1)>b.geom.r(n) & b.R_cyl.T(:,1)<b.geom.r(n+1),1);
        F3.T=F3.T(b.R_cyl.T(:,1)>b.geom.r(n) & b.R_cyl.T(:,1)<b.geom.r(n+1),1);
        F4.T=F4.T(b.R_cyl.T(:,1)>b.geom.r(n) & b.R_cyl.T(:,1)<b.geom.r(n+1),1);

        T_Ssw(:,J) = F2.T;
        T_Nnw(:,J) = F3.T;
        T_w(:,J)   = F.T + F4.T;
        T_Cl(:,J)  = -F.T;

% conductive
    elseif (strcmp(b.bc.r{index_bc_r+n}(1),'w')==1 ||strcmp(b.bc.r{index_bc_r+n}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n}(3),'c')==1
        T_w(:,J)=ones(K,1);

        rhs=str2func(b.bc.rhs.T.r{index_bc_r+n});
        rhs=rhs(b.R_cyl.T(b.R_cyl.T(:,1)>b.geom.r(n) & b.R_cyl.T(:,1)<b.geom.r(n+1),1)).*ones(K,1);
        rhs1(J:b.J+2:end,1)=rhs;

% placeholder surface, interior, periodic
    else

    end

if o==1
    th_e_T_Ssw = [th_e_T_Ssw; T_Ssw];
    th_e_T_w   = [th_e_T_w; T_w];
    th_e_T_Cl  = [th_e_T_Cl; T_Cl];
    th_e_T_Nnw = [th_e_T_Nnw; T_Nnw];
    th_e_T_E   = [th_e_T_E; T_E];

    bc.rhs.Tl  = [bc.rhs.Tl; rhs1];

else
    th_e_T_Sse = [th_e_T_Sse; T_Sse];
    th_e_T_e   = [th_e_T_e; T_w];
    th_e_T_Cr  = [th_e_T_Cr; T_Cl];
    th_e_T_Nne = [th_e_T_Nne; T_Nne];
    th_e_T_W   = [th_e_T_W; T_E];

    bc.rhs.Tr  = [bc.rhs.Tr; rhs1];       
end