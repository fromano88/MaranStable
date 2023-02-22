u_Ssw=zeros(K+1,b.J+2); u_Nnw=u_Ssw; u_w=u_Ssw; u_Cl=u_Ssw; u_E=u_Ssw;
w_sw = zeros(K+1,b.J+3); w_se=w_sw; w_nw=w_sw; w_ne=w_sw;
rhs1=0; rhs=zeros((K+1)*(b.J+2),1);

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
    if (strcmp(b.bc.r{index_bc_r+n}(1),'w')==1 ||strcmp(b.bc.r{index_bc_r+n}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n}(2),'e')==1

        u_w(:,J)   = ones(K+1,1);
        u_Cl(:,J)  = -3/2*ones(K+1,1);
        u_E(:,J)   = 1/2*ones(K+1,1);

        if n>1 && ((strcmp(b.bc.r{index_bc_r+n-1}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n-1}(2),'e')==1)
            u_Ssw(1,J)  = -1/2;
            u_w(1,J)    = 1;
            u_Nnw(1,J)  = -1/2;
            u_Cl(1,J)   = 0;
            u_E(1,J)    = 0;

        end
        if n<length(b.bc.r)/2 && (strcmp(b.bc.r{index_bc_r+n+1}(1),'s')==1 || strcmp(b.bc.r{index_bc_r+n+1}(1),'c')==1 || strcmp(b.bc.r{index_bc_r+n+1}(1),'p')==1)
            u_Ssw(end,J)  = -1/2;
            u_w(end,J)    = 1;
            u_Nnw(end,J)  = -1/2;
            u_Cl(end,J)   = 0;
            u_E(end,J)    = 0;

        end

% slip
    elseif (strcmp(b.bc.r{index_bc_r+n}(1),'w')==1 ||strcmp(b.bc.r{index_bc_r+n}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n}(2),'s')==1
        F.u=b.ETA_z.u(:,J)./b.DETA.rm(:,J-2*(o-1))*2.*b.dr.rm.';
        F1.u=b.XI_z.u(:,J)./b.DXI.rm(:,J-2*(o-1)).*b.dr.rm.';

        F.w=-b.JA.w(:,J-(o-1)).*b.XI_z.w(:,J-(o-1));
        F.v=-b.JA.v(:,J-(o-1):(-1)^(o-1):J-(o-1)-(-1)^o).*b.XI_z.v(:,J-(o-1):(-1)^(o-1):J-(o-1)-(-1)^o);
        F.v(2:end-1,:)=(-1)^o*F.v(2:end-1,:)/2;

        w_sw        = [F.u(1:end-1,:)/2.*F.v(1:end-1,1); 0]-F1.u.*F.w(1:end-1,1);
        w_sw        = [w_sw(b.R_cyl.u(:,1)>=b.geom.r(n) & b.R_cyl.u(:,1)<=b.geom.r(n+1),1) sparse(K+1,b.J+2)];
        w_se        = [-F.u(1:end-1,:)/2.*F.v(1:end-1,2); 0];
        w_se        = [w_se(b.R_cyl.u(:,1)>=b.geom.r(n) & b.R_cyl.u(:,1)<=b.geom.r(n+1),1) sparse(K+1,b.J+2)];
        w_nw        = [0; F.u(2:end,:)/2.*F.v(2:end,1)]+F1.u.*F.w(2:end,1);
        w_nw        = [w_nw(b.R_cyl.u(:,1)>=b.geom.r(n) & b.R_cyl.u(:,1)<=b.geom.r(n+1),1) sparse(K+1,b.J+2)];
        w_ne        = [0; -F.u(2:end,:)/2.*F.v(2:end,2)];
        w_ne        = [w_ne(b.R_cyl.u(:,1)>=b.geom.r(n) & b.R_cyl.u(:,1)<=b.geom.r(n+1),1) sparse(K+1,b.J+2)];

        F.u=F.u(b.R_cyl.u(:,1)>=b.geom.r(n) & b.R_cyl.u(:,1)<=b.geom.r(n+1),1);
        F1.u=F1.u(b.R_cyl.u(:,1)>=b.geom.r(n) & b.R_cyl.u(:,1)<=b.geom.r(n+1),1);

        F.w=b.JA.w(:,J-(o-1)).*b.ETA_z.w(:,J-(o-1));
        F1.w=F.w(b.R_cyl.w(:,1)>b.geom.r(n) & b.R_cyl.w(:,1)<b.geom.r(n+1),1);
        F2.w=F.w(b.R_cyl.u(:,1)==b.geom.r(n),1);
        F2.u=b.JA.u(:,J:(-1)^(o-1):J-(-1)^o).*b.ETA_z.u(:,J:(-1)^(o-1):J-(-1)^o);
        F2.u=(-1)^o*F2.u(b.R_cyl.u(:,1)>=b.geom.r(n) & b.R_cyl.u(:,1)<=b.geom.r(n+1),1:2);

        u_Ssw(:,J) = [0; -F1.u(2:end,1).*F1.w/2];
        u_Nnw(:,J) = [F1.u(1:end-1,1).*F1.w/2; 0];
        u_w(:,J)   = F.u.*F2.u(:,1);
        u_w        = u_w + u_Ssw + u_Nnw;
        u_Cl(:,J)  = -F.u.*F2.u(:,2);

       if n>1 && ((strcmp(b.bc.r{index_bc_r+n-1}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n-1}(2),'s')==1)
            u_Ssw(1,1) = -F1.u(2:end,1).*F2.w/2;
            u_w(1,1)   = u_Ssw(1,1);

        elseif n>1
            u_Ssw(1,J)  = -1/2;
            u_w(1,J)    = 1;
            u_Nnw(1,J)  = -1/2;
            u_Cl(1,J)   = 0;

            w_sw(1,J)   = 0;
            w_se(1,J)   = 0;
            w_nw(1,J)   = 0;
            w_ne(1,J)   = 0;

        end

        if n<length(b.bc.r)/2 && (strcmp(b.bc.r{index_bc_r+n+1}(1),'s')==1 || strcmp(b.bc.r{index_bc_r+n+1}(1),'c')==1 || strcmp(b.bc.r{index_bc_r+n+1}(1),'p')==1)
            u_Ssw(end,J)  = -1/2;
            u_w(end,J)    = 1;
            u_Nnw(end,J)  = -1/2;
            u_Cl(end,J)   = 0;

            w_sw(end,J)   = 0;
            w_se(end,J)   = 0;
            w_nw(end,J)   = 0;
            w_ne(end,J)   = 0;

        end

% no slip
    elseif (strcmp(b.bc.r{index_bc_r+n}(1),'w')==1 ||strcmp(b.bc.r{index_bc_r+n}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n}(2),'n')==1
        u_w(:,J)        = ones(K+1,1);

        w_sw(:,J)   = -b.XI_z.u(:,J)./b.ETA_z.u(:,J)/2;
        w_nw        = w_sw;
        w_sw(1,J)   = w_sw(1,J)*2;
        w_sw(end,J) = 0;
        w_nw(1,J)   = 0;
        w_nw(end,J) = w_nw(end,J)*2;

        if n>1 && ((strcmp(b.bc.r{index_bc_r+n-1}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r+n-1}(2),'n')==1 || strcmp(b.bc.r{index_bc_r+n-1}(2),'m')==1))
            % do nothing

        elseif n>1
            u_Ssw(1,J)  = -1/2;
            u_w(1,J)    = 1;
            u_Nnw(1,J)  = -1/2;

            w_sw(1,J)   = 0;
            w_nw (1,J)  = 0;

        end

        if n<length(b.bc.r)/2 && (strcmp(b.bc.r{index_bc_r+n+1}(1),'s')==1 || strcmp(b.bc.r{index_bc_r+n+1}(1),'c')==1 || strcmp(b.bc.r{index_bc_r+n+1}(1),'p')==1)
            u_Ssw(end,J)  = -1/2;
            u_w(end,J)    = 1;
            u_Nnw(end,J)  = -1/2;

            w_sw(end,J)   = 0;
            w_nw (end,J)  = 0;

        end

        if n>1 && ((strcmp(b.bc.r{index_bc_r+n-1}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n-1}(2),'m')==1)
            rhs1=str2func(b.bc.rhs.u.r{index_bc_r+n-1});
            rhs1=rhs1(b.geom.r(n)).*b.XI_r.u(b.R_cyl.u(:,1)==b.geom.r(n),J)/2;

        end

% moving wall
    elseif (strcmp(b.bc.r{index_bc_r+n}(1),'w')==1 ||strcmp(b.bc.r{index_bc_r+n}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n}(2),'m')==1
        u_w(:,J)        = ones(K+1,1);

        w_sw(:,J)   = -b.XI_z.u(:,J)./b.ETA_z.u(:,J)/2;
        w_nw        = w_sw;
        w_sw(1,J)   = w_sw(1,J)*2;
        w_sw(end,J) = 0;
        w_nw(1,J)   = 0;
        w_nw(end,J) = w_nw(end,J)*2;

        if n>1 && ((strcmp(b.bc.r{index_bc_r+n-1}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r+n-1}(2),'n')==1 || strcmp(b.bc.r{index_bc_r+n-1}(2),'m')==1))
            % do nothing

        elseif n>1
            u_Ssw(1,J)  = -1/2;
            u_w(1,J)    = 1;
            u_Nnw(1,J)  = -1/2;

            w_sw(1,J)   = 0;
            w_nw (1,J)  = 0;

        end
        if n<length(b.bc.r)/2 && (strcmp(b.bc.r{index_bc_r+n+1}(1),'s')==1 || strcmp(b.bc.r{index_bc_r+n+1}(1),'c')==1 || strcmp(b.bc.r{index_bc_r+n+1}(1),'p')==1)
            u_Ssw(end,J)  = -1/2;
            u_w(end,J)    = 1;
            u_Nnw(end,J)  = -1/2;

            w_sw(end,J)   = 0;
            w_nw (end,J)  = 0;

        end

        rhs=str2func(b.bc.rhs.u.r{index_bc_r+n});
        rhs=rhs(b.R_cyl.u(b.R_cyl.u(:,1)>=b.geom.r(n) & b.R_cyl.u(:,1)<=b.geom.r(n+1),1)).*b.XI_r.u(b.R_cyl.u(:,1)>=b.geom.r(n) & b.R_cyl.u(:,1)<=b.geom.r(n+1),J);
        if n>1 && ((strcmp(b.bc.r{index_bc_r+n-1}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n-1}(2),'m')==1)
            rhs1=str2func(b.bc.rhs.u.r{index_bc_r+n-1});
            rhs1=(rhs1(b.geom.r(n)).*b.XI_r.u(b.R_cyl.u(:,1)==b.geom.r(n),J)+rhs(1))/2;

        elseif n>1 && ((strcmp(b.bc.r{index_bc_r+n-1}(1),'w')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'i')==1 || strcmp(b.bc.r{index_bc_r+n-1}(1),'o')==1) && strcmp(b.bc.r{index_bc_r+n-1}(2),'n')==1)
            rhs1=rhs(1)/2;

        end

        rhs2=zeros((K+1)*(b.J+2),1);
        rhs2(J+b.J+2:b.J+2:end-(b.J+2),1)=rhs(2:end-1);
        rhs=rhs2;

% placeholder surface, interior, periodic
    else

    end

rhs(J)=rhs1;

if n==1 || (n>1 && (strcmp(b.bc.r{index_bc_r+n}(1),'s')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'c')==1 || strcmp(b.bc.r{index_bc_r+n}(1),'p')==1))
    u_Ssw(1,:) = [];
    u_w(1,:)   = [];
    u_Cl(1,:)  = [];
    u_Nnw(1,:) = [];
    u_E(1,:)   = [];

    w_sw(1,:)  = [];
    w_se(1,:)  = [];
    w_nw(1,:)  = [];
    w_ne(1,:)  = [];

    rhs(1:b.J+2)   = [];

end
if n==length(b.bc.r)/2 || (n<length(b.bc.r)/2 && (strcmp(b.bc.r{index_bc_r+n+1}(1),'s')~=1 || strcmp(b.bc.r{index_bc_r+n+1}(1),'c')~=1 || strcmp(b.bc.r{index_bc_r+n+1}(1),'p')~=1))
    u_Ssw(end,:) = [];
    u_w(end,:)   = [];
    u_Cl(end,:)  = [];
    u_Nnw(end,:) = [];
    u_E(end,:)   = [];

    w_sw(end,:)  = [];
    w_se(end,:)  = [];
    w_nw(end,:)  = [];
    w_ne(end,:)  = [];

    rhs(end+1-(b.J+2):end)     = [];

end

if o==1
    rm_u_Ssw=[rm_u_Ssw; u_Ssw];
    rm_u_w=[rm_u_w; u_w];
    rm_u_Cl=[rm_u_Cl; u_Cl];
    rm_u_Nnw=[rm_u_Nnw; u_Nnw];
    rm_u_E=[rm_u_E; u_E];

    rm_w_swl=[rm_w_swl; w_sw];
    rm_w_sel=[rm_w_sel; w_se];
    rm_w_nwl=[rm_w_nwl; w_nw];
    rm_w_nel=[rm_w_nel; w_ne];

    bc.rhs.ul=[bc.rhs.ul; rhs];

else
    rm_u_Sse=[rm_u_Sse; u_Ssw];
    rm_u_e=[rm_u_e; u_w];
    rm_u_Cr=[rm_u_Cr; u_Cl];
    rm_u_Nne=[rm_u_Nne; u_Nnw];
    rm_u_W=[rm_u_W; u_E];

    rm_w_swr=[rm_w_swr; w_se];
    rm_w_ser=[rm_w_ser; w_sw];
    rm_w_nwr=[rm_w_nwr; w_ne];
    rm_w_ner=[rm_w_ner; w_nw];

    bc.rhs.ur=[bc.rhs.ur; rhs];

end