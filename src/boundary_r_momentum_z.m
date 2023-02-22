p_s=zeros(J,1);
rm.p= sparse(J+2,b.K*b.J);
u_Wsw=zeros(1,J+2); u_s=u_Wsw; u_Ese=u_Wsw; u_WW=u_Wsw; u_EE=u_Wsw; u_C=u_Wsw; u_N=u_Wsw;
rm.w = sparse(J+2,(b.K+2)*(b.J+1));
rhs1=0; rhs2=0; rhs=zeros(J,1);
C.R=[];
Ja.C.R=[];
counter_C.R=[];
Ja.counter_C.R=[];

if o==1
    % lower boundary
        index_bc_z=0;
        index_bc_r=1;
        index_geom_r=1;
        K=1;
        N.u=N1.u;
        R_z.u=b.R1_z.u;

elseif o==2
    % upper boundary
        index_bc_z=length(b.bc.z)/2;
        index_bc_r=length(b.bc.r)/2;
        index_geom_r=length(b.geom.r);
        K=b.K+1;
        N.u=Nend.u;
        R_z.u=b.Rend_z.u;

end

% wall or axis
    if strcmp(b.bc.z{index_bc_z+n}(1),'w')==1 || (strcmp(b.bc.z{index_bc_z+n}(1),'a')==1 && (stability~=1 || (stability==1 && (m~=1 || ax==0))))
        u_s=ones(1,J+2);

        if n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'m')==1
            rm.w(1,(b.K+1)*(b.J+1)*(o-1)+1)=-b.XI_z.u(K,1)./b.ETA_z.u(K,1)/2;
            rhs1=str2func(b.bc.rhs.u.r{index_bc_r});
            rhs1=rhs1(b.geom.r(index_geom_r)).*b.XI_r.u(K,1)/2;

        elseif n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'n')==1
            rm.w(1,(b.K+1)*(b.J+1)*(o-1)+1)=-b.XI_z.u(K,1)./b.ETA_z.u(K,1)/2;

        end
        if n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'m')==1
            rm.w(J+2,(b.K+1)*(b.J+1)*(o-1)+b.J+1)=-b.XI_z.u(K,end)./b.ETA_z.u(K,end)/2;
            rhs2=str2func(b.bc.rhs.u.r{length(b.bc.r)/2+index_bc_r});
            rhs2=rhs2(b.geom.r(index_geom_r)).*b.XI_r.u(K,end)/2;

        elseif n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'n')==1
            rm.w(J+2,(b.K+1)*(b.J+1)*(o-1)+b.J+1)=-b.XI_z.u(K,end)./b.ETA_z.u(K,end)/2;

        end

% inlet
    elseif strcmp(b.bc.z{index_bc_z+n}(1),'i')==1
        u_s=ones(1,J+2);

        rhs=str2func(b.bc.rhs.u.z{index_bc_z+n});
        if n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'m')==1
            rm.w(1,(b.K+1)*(b.J+1)*(o-1)+1)=-b.XI_z.u(K,1)./b.ETA_z.u(K,1)/2;
            rhs1=str2func(b.bc.rhs.u.r{index_bc_r});
            rhs1=(rhs1(b.geom.r(index_geom_r))+rhs(b.geom.z(1)).*N.u(1,1)).*b.XI_r.u(K,1)/2;

        elseif n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'n')==1
            rm.w(1,(b.K+1)*(b.J+1)*(o-1)+1)=-b.XI_z.u(K,1)./b.ETA_z.u(K,1)/2;
            rhs1=rhs(b.geom.z(1)).*N.u(1,1).*b.XI_r.u(K,1)/2;

        elseif n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r}(2),'s')==1 || strcmp(b.bc.r{index_bc_r}(2),'e')==1)
            rhs1=rhs(b.geom.z(1)).*N.u(1,1).*b.XI_r.u(K,1);

        end
        if n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'m')==1
            rm.w(J+2,(b.K+1)*(b.J+1)*(o-1)+b.J+1)=-b.XI_z.u(K,end)./b.ETA_z.u(K,end)/2;
            rhs2=str2func(b.bc.rhs.u.r{length(b.bc.r)/2+index_bc_r});
            rhs2=(rhs2(b.geom.r(index_geom_r))+rhs(b.geom.z(end))).*N.u(1,end).*b.XI_r.u(K,end)/2;

        elseif n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'n')==1
            rm.w(J+2,(b.K+1)*(b.J+1)*(o-1)+b.J+1)=-b.XI_z.u(K,end)./b.ETA_z.u(K,end)/2;
            rhs2=rhs(b.geom.z(end)).*N.u(1,end).*b.XI_r.u(K,end)/2;

        elseif n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'s')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'e')==1)
            rhs2=rhs(b.geom.z(end)).*N.u(1,end).*b.XI_r.u(K,end);

        end
        rhs=rhs(b.Z.u(1,b.Z.u(1,:)>b.geom.z(n) & b.Z.u(1,:)<b.geom.z(n+1))).'.*ones(J,1).*N.u(1,b.Z.u(1,:)>b.geom.z(n) & b.Z.u(1,:)<b.geom.z(n+1)).'.*b.XI_r.u(K,b.Z.u(1,:)>b.geom.z(n) & b.Z.u(1,:)<b.geom.z(n+1)).';

% outlet
    elseif strcmp(b.bc.z{index_bc_z+n}(1),'o')==1
        rm.u1 = sparse(J+2,(b.K+1)*(b.J+2));
        rm.u1(1,1)=1;
        rm.u1(J+2,b.J+2)=1;
        p_s=ones(J,1);

        if n==1
            if (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'m')==1
                rm.w(1,(b.K+1)*(b.J+1)*(o-1)+1)=-b.XI_z.u(K,1)./b.ETA_z.u(K,1);
                rhs1=str2func(b.bc.rhs.u.r{index_bc_r});
                rhs1=rhs1(b.geom.r(index_geom_r)).*b.XI_r.u(K,1);

            elseif (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'n')==1
                rm.w(1,(b.K+1)*(b.J+1)*(o-1)+1)=-b.XI_z.u(K,1)./b.ETA_z.u(K,1);

            elseif (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r}(2),'s')==1 || strcmp(b.bc.r{index_bc_r}(2),'e')==1)
                u_s(1,1)=1;
                u_Ese(1,1)=-0.75;
                u_EE(1,1)=0.25;
                u_C(1,1)=-1;
                u_N(1,1)=0.5;

            end
        end
        if n==length(b.bc.z)/2
            if (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'m')==1
                rm.w(J+2,(b.K+1)*(b.J+1)*(o-1)+b.J+1)=-b.XI_z.u(K,end)./b.ETA_z.u(K,end);
                rhs2=str2func(b.bc.rhs.u.r{length(b.bc.r)/2+index_bc_r});
                rhs2=rhs2(b.geom.r(index_geom_r)).*b.XI_r.u(K,end);

            elseif (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'n')==1
                rm.w(J+2,(b.K+1)*(b.J+1)*(o-1)+b.J+1)=-b.XI_z.u(K,end)./b.ETA_z.u(K,end);

            elseif (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'s')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'e')==1)
                u_s(1,end)=1;
                u_Wsw(1,end)=-0.75;
                u_WW(1,end)=0.25;
                u_C(1,end)=-1;
                u_N(1,end)=0.5;

            end
        end                        
        if n==1 && strcmp(b.bc.r{index_bc_r}(1),'o')==1
            rm.p(2,(b.K-1)*b.J*(o-1)+1) = -1;
            rm.p(2,(b.K-3)*b.J*(o-1)+b.J+2) = 1;

        end
        if n==length(b.bc.z)/2 && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1
            rm.p(J+1,(b.K-1)*b.J*(o-1)+b.J) = -1;
            rm.p(J+1,(b.K-3)*b.J*(o-1)+2*b.J-1) = 1;

        end
        rhs  = ones(J,1)*b.p_fluid;

% axis
    elseif (strcmp(b.bc.z{index_bc_z+n}(1),'a')==1 && stability==1 && m==1 && ax==1)
        F.u=1./N.u.*(b.XI_r.u(K,:)-R_z.u(1,:).*b.XI_z.u(K,:))./[b.DXI.c(K-(o-1),1) b.DXI.c(K-(o-1),:) b.DXI.c(K-(o-1),end)];
        F.u=(-1)^o*F.u(1,b.Z.u(1,:)>b.geom.z(n) & b.Z.u(1,:)<b.geom.z(n+1));
        F1.u=-1./N.u.*R_z.u(1,:).*b.ETA_z.u(K,:)./[b.DETA.zm(K-(o-1),1) b.DETA.c(K-(o-1),:) b.DETA.zm(K-(o-1),end)];
        F1.u=F1.u(1,b.Z.u(1,:)>b.geom.z(n) & b.Z.u(1,:)<b.geom.z(n+1));

        u_Wsw = [0 -F1.u/2 0];
        u_Ese = [0 F1.u/2 0];
        u_s   = [0 -F.u 0];
        u_C   = [0 F.u 0];

        if n==1
            u_Wsw(1,2) = -F1.u(1,1);
            u_s(1,2)   =  u_s(1,2)+F1.u(1,1)/2;
            if (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r}(2),'s')==1 || strcmp(b.bc.r{index_bc_r}(2),'e')==1)
                u_s(1,1)   = 1;
                u_C(1,1)   = -1;

            elseif (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'m')==1
                u_s(1,1)   = 1;
                rm.w(1,(b.K+1)*(b.J+1)*(o-1)+1)=-b.XI_z.u(K,1)./b.ETA_z.u(K,1);
                rhs1=str2func(b.bc.rhs.u.r{index_bc_r});
                rhs1=rhs1(b.geom.r(index_geom_r)).*b.XI_r.u(K,1);

            elseif (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'n')==1
                u_s(1,1)   = 1;
                rm.w(1,(b.K+1)*(b.J+1)*(o-1)+1)=-b.XI_z.u(K,1)./b.ETA_z.u(K,1);

            end
        end

        if n==length(b.bc.z)/2
            u_Ese(1,end-1) = F1.u(1,end);
            u_s(1,end-1)   = u_s(1,end-1)-F1.u(1,end)/2;
            if (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'s')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'e')==1)
                u_s(1,end)   = 1;
                u_C(1,end)   = -1;

            elseif (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'m')==1
                u_s(1,end)   = 1;
                rm.w(J+2,(b.K+1)*(b.J+1)*(o-1)+b.J+1)=-b.XI_z.u(K,end)./b.ETA_z.u(K,end);
                rhs2=str2func(b.bc.rhs.u.r{length(b.bc.r)/2+index_bc_r});
                rhs2=rhs2(b.geom.r(index_geom_r)).*b.XI_r.u(K,end);

            elseif (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'n')==1
                u_s(1,end)   = 1;
                rm.w(J+2,(b.K+1)*(b.J+1)*(o-1)+b.J+1)=-b.XI_z.u(K,end)./b.ETA_z.u(K,end);

            end
        end                        

% placeholder surface, interior, periodic
    else
        if n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'m')==1
            u_s(1,1)=1;
            rm.w(1,(b.K+1)*(b.J+1)*(o-1)+1)=-b.XI_z.u(K,1)./b.ETA_z.u(K,1);
            rhs1=str2func(b.bc.rhs.u.r{index_bc_r});
            rhs1=rhs1(b.geom.r(index_geom_r)).*b.XI_r.u(K,1);

        elseif n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'n')==1
            u_s(1,1)=1;
            rm.w(1,(b.K+1)*(b.J+1)*(o-1)+1)=-b.XI_z.u(K,1)./b.ETA_z.u(K,1);

        elseif n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'s')==1
            u_s(1,1)=-b.ETA_z.u(K,1)./b.DETA.zm(K-(o-1),1)+(-1)^o*b.XI_z.u(K,1)./b.DXI.c(K-(o-1),1);
            u_Ese(1,1)=b.ETA_z.u(K,1)./b.DETA.zm(K-(o-1),1);
            u_C(1,1)=-(-1)^o*b.XI_z.u(K,1)./b.DXI.c(K-(o-1),1);

        elseif n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'e')==1
            u_s(1,1)=1;
            u_Ese(1,1)=-3/2;
            u_EE(1,1)=1/2;

        end
        if n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'m')==1
            u_s(1,end)=1;
            rm.w(J+2,(b.K+1)*(b.J+1)*(o-1)+b.J+1)=-b.XI_z.u(K,end)./b.ETA_z.u(K,end);
            rhs2=str2func(b.bc.rhs.u.r{length(b.bc.r)/2+index_bc_r});
            rhs2=rhs2(b.geom.r(index_geom_r)).*b.XI_r.u(K,end);

        elseif n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'n')==1
            u_s(1,end)=1;
            rm.w(J+2,(b.K+1)*(b.J+1)*(o-1)+b.J+1)=-b.XI_z.u(K,end)./b.ETA_z.u(K,end);

        elseif n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'s')==1
            u_s(1,end)=b.ETA_z.u(K,end)./b.DETA.zm(K-(o-1),end)+(-1)^o*b.XI_z.u(K,end)./b.DXI.c(K-(o-1),end);
            u_Wsw(1,end)=-b.ETA_z.u(K,end)./b.DETA.zm(K-(o-1),end);
            u_C(1,end)=-(-1)^o*b.XI_z.u(K,end)./b.DXI.c(K-(o-1),end);

        elseif n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'e')==1
            u_s(1,end)=1;
            u_Wsw(1,end)=-3/2;
            u_WW(1,end)=1/2;

        end

        if o==1 && strcmp(b.bc.z{index_bc_z+n}(1),'s')==1 && strcmp(b.bc.z{index_bc_z+n}(2),'s')==1 && eval(['b.' b.bc.z{index_bc_z+n} '.o'])==1
            sizeR1=length(b.z.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)));
            if energy==1
                eval(['Rs.' b.bc.z{index_bc_z+n} '= [b.vector_index+b.sizeR1+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J+(b.K+2)*(b.J+2)+1 b.vector_index+b.sizeR1+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J+(b.K+2)*(b.J+2)+sizeR1];'])

            else
                eval(['Rs.' b.bc.z{index_bc_z+n} '= [b.vector_index+b.sizeR1+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J+1 b.vector_index+b.sizeR1+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J+sizeR1];'])

            end
            b.sizeR1=b.sizeR1+sizeR1;
            b.dbc.R1=[b.dbc.R1; [sparse(sizeR1, length(b.z.u(b.z.u<=b.geom.z(n)))) speye(sizeR1, (b.K+1)*(b.J+2)-length(b.z.u(b.z.u<=b.geom.z(n))))]];

        elseif strcmp(b.bc.z{index_bc_z+n}(1),'s')==1 && strcmp(b.bc.z{index_bc_z+n}(2),'s')==1 && eval(['b.' b.bc.z{index_bc_z+n} '.o'])==1
            sizeRend=length(b.z.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)));
            if energy==1
                eval(['Rs.' b.bc.z{length(b.bc.z)/2+n} '= [b.vector_index+b.sizeRend+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J+(b.K+2)*(b.J+2)+1 b.vector_index+b.sizeRend+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J+(b.K+2)*(b.J+2)+sizeRend];'])

            else
                eval(['Rs.' b.bc.z{length(b.bc.z)/2+n} '= [b.vector_index+b.sizeRend+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J+1 b.vector_index+b.sizeRend+(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J+sizeRend];'])

            end
            b.sizeRend=b.sizeRend+sizeRend;
            b.dbc.Rend=[b.dbc.Rend; [sparse(sizeRend, b.K*(b.J+2)+length(b.z.u(b.z.u<=b.geom.z(n)))) speye(sizeRend, (b.J+2)-length(b.z.u(b.z.u<=b.geom.z(n))))]];

        end

        if (max(strncmp('o', b.bc.z(:),1)) || max(strncmp('o', b.bc.r(:),1)))~=1 && isfield(b, 'V_l')==1 && eval(['b.' b.bc.z{index_bc_z+n} '.o==1'])
            C.R=(-1)^o*(pi*b.R.u(K,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1))).^ax.*b.dz.c(b.Z.p(1,:)>b.geom.z(n) & b.Z.p(1,:)<b.geom.z(n+1));
            Ja.C.R=(-1)^o*(2*pi*b.R.u(K,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1))).^ax.*b.dz.c(b.Z.p(1,:)>b.geom.z(n) & b.Z.p(1,:)<b.geom.z(n+1));

        elseif (max(strncmp('o', b.bc.z(:),1)) || max(strncmp('o', b.bc.r(:),1)))~=1 && isfield(b, 'V_l')==1 && eval(['b.' b.bc.z{index_bc_z+n} '.o==0'])
            counter_C.R=(-1)^o*(pi*b.R.u(K,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1))).^ax.*b.dz.c(b.Z.p(1,:)>b.geom.z(n) & b.Z.p(1,:)<b.geom.z(n+1));
            Ja.counter_C.R=(-1)^o*(2*pi*b.R.u(K,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1))).^ax.*b.dz.c(b.Z.p(1,:)>b.geom.z(n) & b.Z.p(1,:)<b.geom.z(n+1));

        end

    end

if max(strncmp('ss', b.bc.z(:),2))==1 && strncmp('ss', b.bc.z(index_bc_z+n),2)~=1 && (max(strncmp('o', b.bc.z(:),1)) || max(strncmp('o', b.bc.r(:),1)))~=1 && isfield(b, 'V_l')==1
    bc.rhs.C((ri(1)-1)*b.J+zi(1))=bc.rhs.C((ri(1)-1)*b.J+zi(1))-(-1)^o*(pi*b.R.u(K,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1))).^ax.*b.R.u(K,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1))*b.dz.c(b.Z.p(1,:)>b.geom.z(n) & b.Z.p(1,:)<b.geom.z(n+1)).';
end
if o==1
    rm.u = spdiags([u_WW.' u_Wsw.' u_s.' u_Ese.' u_EE.' u_C.' u_N.'],length(b.Z.u(b.Z.u(1,:)<=b.geom.z(n)))-1+[-2 -1 0 1 2 b.J+2 2*(b.J+2)],(J+2),(b.K+1)*(b.J+2));
    rm.p = rm.p + spdiags([0; p_s; 0],length(b.Z.u(b.Z.u(1,:)<=b.geom.z(n)))-2,(J+2),b.K*b.J);

else
    rm.u = spdiags([u_N.' u_C.' u_WW.' u_Wsw.' u_s.' u_Ese.' u_EE.'],(b.K+1)*(b.J+2)-length(b.Z.u(b.Z.u(1,:)>b.geom.z(n)))-1+[-2*(b.J+2) -(b.J+2) -2 -1 0 1 2],(J+2),(b.K+1)*(b.J+2));
    rm.p = rm.p + spdiags([0; p_s; 0],b.K*b.J-length(b.Z.u(b.Z.u(1,:)>b.geom.z(n))),(J+2),b.K*b.J);

end

rhs=[rhs1; rhs; rhs2];

if n~=1
    rm.u(1,:)=[];
    rm.w(1,:)=[];
    rm.p(1,:)=[];
    rhs(1)=[];

end
if n~=length(b.bc.z)/2
    rm.u(end,:)=[];
    rm.w(end,:)=[];
    rm.p(end,:)=[];
    rhs(end)=[];

end
if o==1
    bc.rm.u1          = [bc.rm.u1; rm.u];
    bc.rm.w1          = [bc.rm.w1; rm.w];
    bc.rm.p1          = [bc.rm.p1; rm.p];
    bc.rhs.u1         = [bc.rhs.u1; rhs];
    b.C.R1            = [b.C.R1 C.R];
    b.counter.C.R1    = [b.counter.C.R1 counter_C.R];
    b.Ja.C.R1         = [b.Ja.C.R1 Ja.C.R];
    b.Ja.counter.C.R1 = [b.Ja.counter.C.R1 Ja.counter_C.R];

else
    bc.rm.uend          = [bc.rm.uend; rm.u];
    bc.rm.wend          = [bc.rm.wend; rm.w];
    bc.rm.pend          = [bc.rm.pend; rm.p];
    bc.rhs.uend         = [bc.rhs.uend; rhs];
    b.C.Rend            = [b.C.Rend C.R];
    b.counter.C.Rend    = [b.counter.C.Rend counter_C.R];
    b.Ja.C.Rend         = [b.Ja.C.Rend Ja.C.R];
    b.Ja.counter.C.Rend = [b.Ja.counter.C.Rend Ja.counter_C.R];

end