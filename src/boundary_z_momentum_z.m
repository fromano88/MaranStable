w_s=zeros(1,J+1); w_C=w_s; w_N=w_s; w_Wsw=w_s; w_Ese=w_s; w_EE=w_s; w_WW=w_s;
u_sw=zeros(1,J+1); u_se=u_sw; u_nw=u_sw; u_ne=u_sw;
rhs1=0; rhs2=0; rhs=sparse(J+1,1);

if o==1
    % lower boundary
        index_bc_z=0;
        index_bc_r=1;
        index_geom_r=1;
        K=1;
        N.u=N1.u;
        N.w=N1.w;
        R_z.u=b.R1_z.u;
        R_z.w=b.R1_z.w;
        R_z.v=b.R1_z.v;

elseif o==2
    % upper boundary
        index_bc_z=length(b.bc.z)/2;
        index_bc_r=length(b.bc.r)/2;
        index_geom_r=length(b.geom.r);
        K=b.K+2;
        N.u=Nend.u;
        N.w=Nend.w;
        R_z.u=b.Rend_z.u;
        R_z.w=b.Rend_z.w;
        R_z.v=b.Rend_z.v;

end

% extrapolated
    if ((strcmp(b.bc.z{index_bc_z+n}(1),'w')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n}(2),'e')==1)
        w_s=ones(1,J+1);
        w_C=-3/2*ones(1,J+1);
        w_N=1/2*ones(1,J+1);

        if n>1 && ((strcmp(b.bc.z{index_bc_z+n-1}(1),'w')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n-1}(2),'e')~=1)
            w_Wsw(1,1)  = -1/2;
            w_s(1,1)    = 1;
            w_Ese(1,1)  = -1/2;
            w_C(1,1)    = 0;
            w_N(1,1)    = 0;

        elseif n==1 && strcmp(b.bc.r{index_bc_r}(1),'o')==1
            w_s(1,1)   = 1;
            w_Ese(1,1) = -1;
            w_EE(1,1)  = 0.5;
            w_C(1,1)    = -0.75;
            w_N(1,1)    = 0.25;

        elseif n==1
            w_s(1,1)    = 1;
            w_C(1,1)    = 0;
            w_N(1,1)    = 0;

        end

        if n<length(b.bc.z)/2 && (strcmp(b.bc.z{index_bc_z+n+1}(1),'s')==1 || strcmp(b.bc.z{index_bc_z+n+1}(1),'c')==1 || strcmp(b.bc.z{index_bc_z+n+1}(1),'p')==1)
            w_s(1,end)    = 1;
            w_Wsw(1,end)  = -1/2;
            w_Ese(1,end)  = -1/2;
            w_C(1,end)    = 0;
            w_N(1,end)    = 0;

        elseif n==length(b.bc.z)/2 && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1
            w_s(1,end)   = 1;
            w_Wsw(1,end) = -1;
            w_WW(1,end)  = 0.5;
            w_C(1,end)    = -0.75;
            w_N(1,end)    = 0.25;

        elseif n==length(b.bc.z)/2
            w_s(1,end)    = 1;
            w_C(1,end)    = 0;
            w_N(1,end)    = 0;

        end

        if n==1 && strcmp(b.bc.r{index_bc_r}(1),'i')==1
            rhs1=str2func(b.bc.rhs.w.r{index_bc_r});
            rhs1=rhs1(b.geom.r(index_geom_r)).*b.ETA_z.w(K,1);

        end

        if n==length(b.bc.z)/2 && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1
            rhs2=str2func(b.bc.rhs.w.r{length(b.bc.r)/2+index_bc_r});
            rhs2=rhs2(b.geom.r(index_geom_r)).*b.ETA_z.w(K,end);

        end

% slip or axis
    elseif ((strcmp(b.bc.z{index_bc_z+n}(1),'w')==1 ||strcmp(b.bc.z{index_bc_z+n}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n}(2),'s')==1) || (strcmp(b.bc.z{index_bc_z+n}(1),'a')==1 && (stability~=1 || (stability==1 && (m==0 || ax==0))))
        F.u=N.u.*b.JA.u(K-(o-1),:).*b.XI_r.u(K-(o-1),:);
        F.w=N.w.*b.XI_r.w(K,:)./b.DXI.zm(K-2*(o-1),:)*2;
        F.v=([1;1]*(1./N.w)).*b.JA.v(K-(o-1):(-1)^(o-1):K-(o-1)-(-1)^o,:).*R_z.v(K-(o-1):(-1)^(o-1):K-(o-1)-(-1)^o,:).*b.ETA_z.v(K-(o-1):(-1)^(o-1):K-(o-1)-(-1)^o,:);
        F.v(:,2:end-1)=(-1)^o*F.v(:,2:end-1)/2;
        
        F1.w=-1./N.w.*R_z.w(K,:).*b.ETA_z.w(K,:)./b.DETA.zm(K-2*(o-1),:);
        
        F2.w=(-1)^o*([1;1]*(1./N.w)).*b.JA.w(K:(-1)^(o-1):K-(-1)^o,:).*(b.XI_r.w(K:(-1)^(o-1):K-(-1)^o,:)-R_z.w(K:(-1)^(o-1):K-(-1)^o,:).*b.XI_z.w(K:(-1)^(o-1):K-(-1)^o,:));
        F2.u=F.u(:,b.Z.w(1,:)==b.geom.z(n));
        
        F3.u=1./N.u.*b.JA.u(K-(o-1),:).*R_z.u(K-(o-1),:).*b.ETA_z.u(K-(o-1),:);

        u_sw  = ([F.w(1,1:end-1)/2.*F.v(1,1:end-1) 0] -F1.w.*F3.u(1,1:end-1)).*b.dz.zm.*N.w;
        u_se  = ([0 F.w(1,2:end)/2.*F.v(1,2:end)] +F1.w.*F3.u(1,2:end)).*b.dz.zm.*N.w;
        u_nw  = ([-F.w(1,1:end-1)/2.*F.v(2,1:end-1) 0]).*b.dz.zm.*N.w;
        u_ne  = ([0 -F.w(1,2:end)/2.*F.v(2,2:end)]).*b.dz.zm.*N.w;

        w_Wsw = -F1.w.*F.u(2:end)/2.*b.dz.zm.*N.w;
        w_Ese = F1.w.*F.u(1:end-1)/2.*b.dz.zm.*N.w;
        w_s   = F.w.*F2.w(1,:).*b.dz.zm.*N.w+w_Wsw+w_Ese;
        w_C   = -F.w.*F2.w(2,:).*b.dz.zm.*N.w;
        
        u_sw  = u_sw(:,b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1));
        u_se  = u_se(:,b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1));
        u_nw  = u_nw(:,b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1));
        u_ne  = u_ne(:,b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1));
        
        w_Wsw = w_Wsw(:,b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1));
        w_Ese = w_Ese(:,b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1));
        w_s   = w_s(:,b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1));
        w_C   = w_C(:,b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1));

        if n>1 && (((strcmp(b.bc.z{index_bc_z+n-1}(1),'w')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n-1}(2),'s')==1) || strcmp(b.bc.z{index_bc_z+n-1}(1),'a')==1)
            w_Wsw(1,1)  = -F1.w(1,1).*F2.u/2.*b.dz.zm(1,1).*N.w(1,1);
            w_s(1,1)    = w_s(1,1)+w_Wsw(1,1);

        elseif n>1
            u_sw(1,1)   = 0;
            u_se(1,1)   = 0;
            u_nw(1,1)   = 0;
            u_ne(1,1)   = 0;

            w_Wsw(1,1)  = -1/2;
            w_s(1,1)    = 1;
            w_Ese(1,1)  = -1/2;
            w_C(1,1)    = 0;

        elseif n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1)
            u_sw(1,1)   = 0;
            u_se(1,1)   = 0;
            u_nw(1,1)   = 0;
            u_ne(1,1)   = 0;

            w_Wsw(1,1)  = 0;
            w_s(1,1)    = 1;
            w_Ese(1,1)  = 0;
            w_C(1,1)    = 0;

            rhs1=str2func(b.bc.rhs.w.r{index_bc_r});
            rhs1=rhs1(b.geom.r(index_geom_r)).*b.ETA_z.w(K,1);

        end

        if n<length(b.bc.z)/2 && (strcmp(b.bc.z{index_bc_z+n+1}(1),'s')==1 || strcmp(b.bc.z{index_bc_z+n+1}(1),'c')==1 || strcmp(b.bc.z{index_bc_z+n+1}(1),'p')==1)
            u_sw(1,end)   = 0;
            u_se(1,end)   = 0;
            u_nw(1,end)   = 0;
            u_ne(1,end)   = 0;

            w_s(1,end)    = 1;
            w_Wsw(1,end)  = -1/2;
            w_Ese(1,end)  = -1/2;
            w_C(1,end)    = 0;

        elseif n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1)
            u_sw(1,end)   = 0;
            u_se(1,end)   = 0;
            u_nw(1,end)   = 0;
            u_ne(1,end)   = 0;

            w_Wsw(1,end)  = 0;
            w_s(1,end)    = 1;
            w_C(1,end)    = 0;
            w_Ese(1,end)  = 0;

            rhs2=str2func(b.bc.rhs.w.r{length(b.bc.r)/2+index_bc_r});
            rhs2=rhs2(b.geom.r(index_geom_r)).*b.ETA_z.w(K,end);

        end

% no slip or axis
    elseif ((strcmp(b.bc.z{index_bc_z+n}(1),'w')==1 ||strcmp(b.bc.z{index_bc_z+n}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n}(2),'n')==1) || (strcmp(b.bc.z{index_bc_z+n}(1),'a')==1 && stability==1 && (m~=0 || ax==0))
        u_sw        = R_z.w(K,:).*b.ETA_z.w(K,:)./N.w.^2./b.XI_r.w(K,:)/2;
        u_se        = u_sw;
        u_sw(1,1)   = u_sw(1,1)*2;
        u_sw(1,end) = 0;
        u_se(1,end) = u_se(1,end)*2;
        u_se(1,1)   = 0;
        u_sw        = u_sw(1,b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1));
        u_se        = u_se(1,b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1));

        w_s        = ones(1,J+1);

        if n>1 && ((strcmp(b.bc.z{index_bc_z+n-1}(1),'w')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'o')==1) && (strcmp(b.bc.z{index_bc_z+n-1}(2),'n')==1 || strcmp(b.bc.z{index_bc_z+n-1}(2),'m')==1))
            % do nothing

        elseif n>1
            u_sw(1,1)   = 0;
            u_se(1,1)   = 0;

            w_s(1,1)    = 1;
            w_Wsw(1,1)  = -1/2;
            w_Ese(1,1)  = -1/2;

        end

        if n<length(b.bc.z)/2 && (strcmp(b.bc.z{index_bc_z+n+1}(1),'s')==1 || strcmp(b.bc.z{index_bc_z+n+1}(1),'c')==1 || strcmp(b.bc.z{index_bc_z+n+1}(1),'p')==1)
            u_sw(1,end)   = 0;
            u_se(1,end)   = 0;

            w_s(1,end)    = 1;
            w_Wsw(1,end)  = -1/2;
            w_Ese(1,end)  = -1/2;

        end

        if n>1 && ((strcmp(b.bc.z{index_bc_z+n-1}(1),'w')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n-1}(2),'m')==1)
            rhs1=str2func(b.bc.rhs.w.z{index_bc_z+n-1});
            rhs1=rhs1(b.geom.z(n))./N.w(1,b.Z.w(1,:)==b.geom.z(n)).*b.ETA_z.w(K,b.Z.w(1,:)==b.geom.z(n))/2;

        elseif n==1 && strcmp(b.bc.r{index_bc_r}(1),'i')==1
            w_s(1,1)    = 2*w_s(1,1);
            rhs1=str2func(b.bc.rhs.w.r{index_bc_r});
            rhs1=rhs1(b.geom.r(index_geom_r)).*b.ETA_z.w(K,1);

        end

        if n==length(b.bc.z)/2 && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1
            w_s(1,end)    = 2*w_s(1,end);
            rhs2=str2func(b.bc.rhs.w.r{length(b.bc.r)/2+index_bc_r});
            rhs2=rhs2(b.geom.r(index_geom_r)).*b.ETA_z.w(K,end);

        end  

% moving wall
    elseif (strcmp(b.bc.z{index_bc_z+n}(1),'w')==1 ||strcmp(b.bc.z{index_bc_z+n}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n}(2),'m')==1
        u_sw        = R_z.w(K,:).*b.ETA_z.w(K,:)./N.w.^2./b.XI_r.w(K,:)/2;
        u_se        = u_sw;
        u_sw(1,1)   = u_sw(1,1)*2;
        u_sw(1,end) = 0;
        u_se(1,end) = u_se(1,end)*2;
        u_se(1,1)   = 0;
        u_sw        = u_sw(1,b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1));
        u_se        = u_se(1,b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1));

        w_s        = ones(1,J+1);

        if n>1 && ((strcmp(b.bc.z{index_bc_z+n-1}(1),'w')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'o')==1) && (strcmp(b.bc.z{index_bc_z+n-1}(2),'n')==1 || strcmp(b.bc.z{index_bc_z+n-1}(2),'m')==1))
            % do nothing

        elseif n>1
            u_sw(1,1)   = 0;
            u_se(1,1)   = 0;

            w_Wsw(1,1)  = -1/2;
            w_Ese(1,1)  = -1/2;

        end

        if n<length(b.bc.z)/2 && (strcmp(b.bc.z{index_bc_z+n+1}(1),'s')==1 || strcmp(b.bc.z{index_bc_z+n+1}(1),'c')==1 || strcmp(b.bc.z{index_bc_z+n+1}(1),'p')==1)
            u_sw(1,end)   = 0;
            u_se(1,end)   = 0;

            w_Wsw(1,end)  = -1/2;
            w_Ese(1,end)  = -1/2;

        end

        rhs_w=str2func(b.bc.rhs.w.z{index_bc_z+n});
        rhs=(rhs_w(b.Z.w(1,:))./N.w.*b.ETA_z.w(K,:)).';
        rhs=rhs(b.Z.w(1,:)>=b.geom.z(n) & b.Z.w(1,:)<=b.geom.z(n+1),1);

        if n>1 && ((strcmp(b.bc.z{index_bc_z+n-1}(1),'w')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n-1}(2),'m')==1)
            rhs1=str2func(b.bc.rhs.w.z{index_bc_z+n-1});
            rhs1=(rhs1(b.geom.z(n))./N.w(1,b.Z.w(1,:)==b.geom.z(n)).*b.ETA_z.w(K,b.Z.w(1,:)==b.geom.z(n))+rhs(1))/2;

        elseif n>1 && ((strcmp(b.bc.z{index_bc_z+n-1}(1),'w')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n-1}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n-1}(2),'n')==1)
            rhs1=rhs(1)/2;

        elseif n==1 && strcmp(b.bc.r{index_bc_r}(1),'i')==1
            w_s(1,1)    = 2*w_s(1,1);
            rhs1=str2func(b.bc.rhs.w.r{index_bc_r});
            rhs1=rhs1(b.geom.r(index_geom_r)).*b.ETA_z.w(K,1)+rhs(1);

        elseif n==1 && strcmp(b.bc.r{index_bc_r}(1),'w')==1
            w_s(1,1)    = 2*w_s(1,1);
            rhs1=rhs(1);

        elseif n==1
            rhs1=rhs(1);

        end
        if n==length(b.bc.z)/2 && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1
            w_s(1,end)    = 2*w_s(1,end);
            rhs2=str2func(b.bc.rhs.w.r{length(b.bc.r)/2+index_bc_r});
            rhs2=rhs2(b.geom.r(index_geom_r)).*b.ETA_z.w(K,end)+rhs(end);

        elseif n==length(b.bc.z)/2 && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1
            w_s(1,end)    = 2*w_s(1,end);
            rhs2=rhs(end);

        elseif n==length(b.bc.z)/2
            rhs2=rhs(end);

        end

% placeholder surface, interior, periodic
    else
        if n==1 && strcmp(b.bc.r{index_bc_r}(1),'w')==1
            w_s(1,1)=1;

        elseif n==1 && strcmp(b.bc.r{index_bc_r}(1),'i')==1
            w_s(1,1)=1;
            rhs1=str2func(b.bc.rhs.w.r{index_bc_r});
            rhs1=rhs1(b.geom.r(index_geom_r)).*b.ETA_z.w(K,1);

        elseif n==1 && strcmp(b.bc.r{index_bc_r}(1),'o')==1
            w_s(1,1)=1;
            w_Ese(1,1)=-2;
            w_EE(1,1)=1;

        end
        if n==length(b.bc.z)/2 && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1
            w_s(1,end)=1;

        elseif n==length(b.bc.z)/2 && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1
            w_s(1,end)=1;
            rhs2=str2func(b.bc.rhs.w.r{length(b.bc.r)/2+index_bc_r});
            rhs2=rhs2(b.geom.r(index_geom_r)).*b.ETA_z.w(K,end);

        elseif n==length(b.bc.z)/2 && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1
            w_s(1,end)=1;
            w_Wsw(1,end)=-2;
            w_WW(1,end)=1;

        end

    end

if o==1
    bc.zm.u  = spdiags([u_sw.' u_se.' u_nw.' u_ne.'],length(b.Z.w(b.Z.w(1,:)<b.geom.z(n)))+[0 1 b.J+2 b.J+3],(J+1),(b.K+1)*(b.J+2));
    bc.zm.w  = spdiags([w_WW.' w_Wsw.' w_s.' w_Ese.' w_EE.' w_C.' w_N.'],length(b.Z.w(b.Z.w(1,:)<b.geom.z(n)))+[-2 -1 0 1 2 b.J+1 2*(b.J+1)],(J+1),(b.K+2)*(b.J+1));
    bc.rhs.w = [rhs1; rhs(2:end-1); rhs2];

else
    bc.zm.u  = spdiags([u_nw.' u_ne.' u_sw.' u_se.'],b.K*(b.J+2)-1-length(b.Z.u(b.Z.u(1,:)>b.geom.z(n)))+[0 1 b.J+2 b.J+3],(J+1),(b.K+1)*(b.J+2));
    bc.zm.w  = spdiags([w_N.' w_C.' w_WW.' w_Wsw.' w_s.' w_Ese.' w_EE.'],(b.K+2)*(b.J+1)-length(b.Z.w(b.Z.w(1,:)>=b.geom.z(n)))+[-2*(b.J+1) -(b.J+1) -2 -1 0 1 2],(J+1),(b.K+2)*(b.J+1));
    bc.rhs.w = [rhs1; rhs(2:end-1); rhs2];

end
if n>1 && (strcmp(b.bc.z{index_bc_z+n}(1),'s')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'c')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'p')==1)
    bc.zm.u(1,:)  = [];
    bc.zm.w(1,:)  = [];
    bc.rhs.w(1,:) = [];

end
if n<length(b.bc.z)/2 && (strcmp(b.bc.z{index_bc_z+n+1}(1),'s')~=1 && strcmp(b.bc.z{index_bc_z+n+1}(1),'c')~=1 && strcmp(b.bc.z{index_bc_z+n+1}(1),'p')~=1)
    bc.zm.u(end,:)  = [];
    bc.zm.w(end,:)  = [];
    bc.rhs.w(end,:) = [];

end
if o==1
    bc.zm.u1 = [bc.zm.u1; bc.zm.u];
    bc.zm.w1  = [bc.zm.w1; bc.zm.w];
    bc.rhs.w1 = [bc.rhs.w1; bc.rhs.w];

else
    bc.zm.uend = [bc.zm.uend; bc.zm.u];
    bc.zm.wend  = [bc.zm.wend; bc.zm.w];
    bc.rhs.wend = [bc.rhs.wend; bc.rhs.w];

end