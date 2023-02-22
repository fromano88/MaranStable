T_WW=zeros(1,J+2); T_Wsw=T_WW; T_s=T_WW; T_Ese=T_WW; T_EE=T_WW; T_C=T_WW; T_N=T_WW;
rhs1=0; rhs2=0; rhs=sparse(J,1);

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
    if ((strcmp(b.bc.z{index_bc_z+n}(1),'w')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n}(3),'e')==1)
        T_s=ones(1,J+2);
        T_C=-1.5*ones(1,J+2);%-(b.R_cyl.T(K,1:J+2)-b.R_cyl.T(K-2*(-1)^o,1:J+2))./(b.R_cyl.T(K-(-1)^o,1:J+2)-b.R_cyl.T(K-2*(-1)^o,1:J+2));
        T_N=0.5*ones(1,J+2);%-(b.R_cyl.T(K-(-1)^o,1:J+2)-b.R_cyl.T(K,1:J+2))./(b.R_cyl.T(K-(-1)^o,1:J+2)-b.R_cyl.T(K-2*(-1)^o,1:J+2));

        if n==1 && ((strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(3),'c')==1)
            T_C(1,1)   = 0;
            T_N(1,1)   = 0;
            rhs1=str2func(b.bc.rhs.T.r{index_bc_r});
            rhs1=rhs1(b.geom.r(index_geom_r));

        end
        if n==length(b.bc.z)/2 && ((strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(3),'c')==1)
            T_C(1,end)   = 0;
            T_N(1,end)   = 0;
            rhs2=str2func(b.bc.rhs.T.r{length(b.bc.r)/2+index_bc_r});
            rhs2=rhs2(b.geom.r(index_geom_r));

        end                        

% adiabatic or axis
    elseif ((strcmp(b.bc.z{index_bc_z+n}(1),'w')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n}(3),'a')==1) || (strcmp(b.bc.z{index_bc_z+n}(1),'a')==1 && (stability~=1 || (stability==1 && (m==0 || ax==0))))
        F.T=N.u.^2.*b.XI_r.T(K,:)./[b.DXI.rm(K-(o-1),1) b.DXI.rm(K-(o-1),:) b.DXI.rm(K-(o-1),end)];
        F.T=(-1)^o*F.T(2:end-1).*b.dz.c;
        F1.T=-R_z.T(K,:).*b.ETA_z.T(K,:)./[b.DETA.zm(1,1) b.DETA.c(1,:) b.DETA.zm(1,end)];
        F1.T=F1.T(2:end-1).*b.dz.c;

        T_Wsw = [0 -b.W.w(1,:).*F1.T 0];
        T_Wsw = [0 T_Wsw(1,b.Z.T(1,:)>b.geom.z(n) & b.Z.T(1,:)<b.geom.z(n+1)) 0];
        T_Ese = [0 b.E.w(1,:).*F1.T 0];
        T_Ese = [0 T_Ese(1,b.Z.T(1,:)>b.geom.z(n) & b.Z.T(1,:)<b.geom.z(n+1)) 0];
        T_s   = [0 (F.T - [0 b.E.w(1,1:end-1)].*F1.T + [b.W.w(1,2:end) 0].*F1.T) 0];
        T_s   = [0 T_s(1,b.Z.T(1,:)>b.geom.z(n) & b.Z.T(1,:)<b.geom.z(n+1)) 0];
        T_C   = [0 -F.T 0];
        T_C   = [0 T_C(1,b.Z.T(1,:)>b.geom.z(n) & b.Z.T(1,:)<b.geom.z(n+1)) 0];
        T_EE  = zeros(1,J+2);
        T_WW  = zeros(1,J+2);
        if n==1 && ((strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(3),'a')==1)
            T_s(1,1)   = 1;
            T_C(1,1)   = -1;

        elseif n==1 && ((strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(3),'e')==1)
            T_s(1,1)   = 1;
            T_Ese(1,1) = -1.5;%-(b.Z.T(1,1)-b.Z.T(1,3))./(b.Z.T(1,2)-b.Z.T(1,3));
            T_EE(1,1)  = 0.5;%-(b.Z.T(1,2)-b.Z.T(1,1))./(b.Z.T(1,2)-b.Z.T(1,3));

        elseif n==1 && ((strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(3),'c')==1)
            T_s(1,1)   = 1;
            rhs1=str2func(b.bc.rhs.T.r{index_bc_r});
            rhs1=rhs1(b.geom.r(index_geom_r));

        end

        if n==length(b.bc.z)/2 && ((strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(3),'a')==1)
            T_s(1,end)   = 1;
            T_C(1,end)   = -1;

        elseif n==length(b.bc.z)/2 && ((strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(3),'e')==1)
            T_s(1,end)   = 1;
            T_Wsw(1,end) = -1.5;%-(b.Z.T(1,end)-b.Z.T(1,end-2))./(b.Z.T(1,end-1)-b.Z.T(1,end-2));
            T_WW(1,end)  = 0.5;%-(b.Z.T(1,end-1)-b.Z.T(1,end))./(b.Z.T(1,end-1)-b.Z.T(1,end-2));

        elseif n==length(b.bc.z)/2 && ((strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(3),'c')==1)
            T_s(1,end)   = 1;
            rhs2=str2func(b.bc.rhs.T.r{length(b.bc.r)/2+index_bc_r});
            rhs2=rhs2(b.geom.r(index_geom_r));

        end

% conductive
    elseif ((strcmp(b.bc.z{index_bc_z+n}(1),'w')==1 ||strcmp(b.bc.z{index_bc_z+n}(1),'i')==1 || strcmp(b.bc.z{index_bc_z+n}(1),'o')==1) && strcmp(b.bc.z{index_bc_z+n}(3),'c')==1) || (strcmp(b.bc.z{index_bc_z+n}(1),'a')==1 && stability==1 && m~=0)
        T_s=ones(1,J+2);

        rhs=str2func(b.bc.rhs.T.z{index_bc_z+n});
        if n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(3),'c')==1
            rhs1=str2func(b.bc.rhs.T.r{index_bc_r});
            rhs1=(rhs1(b.geom.r(index_geom_r))+rhs(b.geom.z(1)))/2;

        elseif n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && (strcmp(b.bc.r{index_bc_r}(3),'a')==1 || strcmp(b.bc.r{index_bc_r}(3),'e')==1)
            rhs1=rhs(b.geom.z(1));

        end
        if n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(3),'c')==1
            rhs2=str2func(b.bc.rhs.T.r{length(b.bc.r)/2+index_bc_r});
            rhs2=(rhs2(b.geom.r(index_geom_r))+rhs(b.geom.z(end)))/2;

        elseif n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(3),'a')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(3),'e')==1)
            rhs2=rhs(b.geom.z(end));

        end
        rhs=rhs(b.Z.u(1,b.Z.u(1,:)>b.geom.z(n) & b.Z.u(1,:)<b.geom.z(n+1))).'.*ones(J,1);

% placeholder surface, interior, periodic
    else
        if o == 2 && strcmp(b.bc.z{index_bc_z+n}(3),'c') == 1 && str2double(b.bc.z{o}(4)) == 0  % SFM
            l_lb  = diff(b.geom.z);
            r_i   = b.r.u(end);
            r_c   = b.r.u(1);
            T_d1l = str2func(b.bc.rhs.T.r{1}); T_d1l = T_d1l(r_i);
            T_d2l = str2func(b.bc.rhs.T.r{2}); T_d2l = T_d2l(r_i);
            T_ref = min(T_d1l,T_d2l);
            
            lambda = str2func(b.lambda);
            if strcmp(b.bc.rhs.T.z{index_bc_z+n}(1),'@')
                rhs = str2func(b.bc.rhs.T.z{index_bc_z+n});
                rhs = rhs(b.Z.T(end,2:end-1))'.*lambda(b.T.T(end,2:end-1))'/l_lb.*T_ref.*b.dz.c'; % heat transfer coefficient
            else
                subfolder = b.bc.rhs.T.z{index_bc_z+n};

                Gamma  = l_lb/r_i;
                T_0    = (T_d1l+T_d2l)/2;
                rho    = str2func(b.rho);
                mu     = str2func(b.mu);
                cp     = str2func(b.cp);
                dsigma = str2func(b.dsigma);

                Re  = abs(dsigma(T_0)*(T_d1l-T_d2l))*l_lb*rho(T_0)/mu(T_0)^2;
                Pr  = mu(T_0)*cp(T_0)/lambda(T_0);
                Bi  = Biot_func(Pr,flowopt.g,Gamma,Re,subfolder);% % % 
                rhs = Bi(0.5-b.Z.T(end,2:end-1)/l_lb).*lambda(b.T.T(end,2:end-1))/l_lb; % heat transfer coefficient

                % smoothing
                sm.bound  = 0.49;
                sm.z_orig = (l_lb/2-(b.Z.T(end,2:end-1)))/l_lb;
                sm.z      = sm.z_orig(sm.z_orig<=sm.bound & sm.z_orig>=-sm.bound);
                sm.rhs    = smooth(sm.z,rhs((sm.z_orig<=sm.bound & sm.z_orig>-sm.bound)),5,'moving');
                rhs = [rhs(sm.z_orig>sm.bound) sm.rhs' rhs(sm.z_orig<-sm.bound)];
                
                % constant Bi function
                %{
                %Bi = load(['BiotFunction_coeffs/' subfolder '/Bi_avg_Pr28.84_Gamma' num2str(Gamma,'%3.2f') '_g0.dat']);
                %Bi = load(['BiotFunction_coeffs/' subfolder '/Bi_avg_lin_Pr28.84_Gamma' num2str(Gamma,'%3.2f') '_g0.dat']);
                %Bi_avg = interp1(Bi(:,1),Bi(:,2),Re,'linear');
                Bi_avg = BiotAvgLin(Pr,flowopt.g,Gamma,Re,subfolder);
                rhs = Bi_avg*ones(size(b.Z.T(end,2:end-1))).*lambda(b.T.T(end,2:end-1))/l_lb;
                %}
            
                % linear Bi function
                %{
                Bi = load(['BiotFunction_coeffs/' subfolder '/Bi_lin_Pr28.84_Gamma' num2str(Gamma,'%3.2f') '_g0.dat']);            
                k = interp1(Bi(:,1),Bi(:,3),Re,'linear');
                d = interp1(Bi(:,1),Bi(:,2),Re,'linear');
                Bi_lin = @(z)k*z+d;
                rhs = Bi_lin(0.5-b.Z.T(end,2:end-1)/l_lb).*lambda(b.T.T(end,2:end-1))/l_lb;
                %}
                
                rhs = rhs'*T_ref.*b.dz.c';
            end
        end
        if n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(3),'c')==1
            T_s(1,1)=1;
            rhs1=str2func(b.bc.rhs.T.r{index_bc_r});
            rhs1=rhs1(b.geom.r(index_geom_r));

        elseif n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(3),'a')==1
            T_s(1,1)   = -b.ETA_z.T(K,1)./b.DETA.zm(1,1)+(-1)^o*b.XI_z.T(K,1)./b.DXI.rm(K-(o-1),1);
            T_Ese(1,1) = b.ETA_z.T(K,1)./b.DETA.zm(1,1);
            T_C(1,1)   = -(-1)^o*b.XI_z.T(K,1)./b.DXI.rm(K-(o-1),1);

        elseif n==1 && (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(3),'e')==1
            T_s(1,1)   = 1;
            T_Ese(1,1) = -1.5;%-(b.Z.T(1,1)-b.Z.T(1,3))./(b.Z.T(1,2)-b.Z.T(1,3));
            T_EE(1,1)  = 0.5;%-(b.Z.T(1,2)-b.Z.T(1,1))./(b.Z.T(1,2)-b.Z.T(1,3));

        end
        if n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(3),'c')==1
            T_s(1,end)=1;
            rhs2=str2func(b.bc.rhs.T.r{length(b.bc.r)/2+index_bc_r});
            rhs2=rhs2(b.geom.r(index_geom_r));

        elseif n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(3),'a')==1
            T_s(1,end)   = b.ETA_z.T(K,end)./b.DETA.zm(1,end)+(-1)^o*b.XI_z.T(K,end)./b.DXI.rm(K-(o-1),end);
            T_Wsw(1,end) = -b.ETA_z.T(K,end)./b.DETA.zm(1,end);
            T_C(1,end)   = -(-1)^o*b.XI_z.T(K,end)./b.DXI.rm(K-(o-1),end);

        elseif n==length(b.bc.z)/2 && (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(3),'e')==1
            T_s(1,end)   = 1;
            T_Wsw(1,end) = -1.5;%-(b.Z.T(1,end)-b.Z.T(1,end-2))./(b.Z.T(1,end-1)-b.Z.T(1,end-2));
            T_WW(1,end)  = 0.5;%-(b.Z.T(1,end-1)-b.Z.T(1,end))./(b.Z.T(1,end-1)-b.Z.T(1,end-2));

        end

    end

if o==1
    th_e.T = spdiags([T_WW.' T_Wsw.' T_s.' T_Ese.' T_EE.' T_C.' T_N.'],length(b.Z.T(b.Z.T(1,:)<=b.geom.z(n)))-1+[-2 -1 0 1 2 b.J+2 2*(b.J+2)],(J+2),(b.K+2)*(b.J+2));

else
    th_e.T = spdiags([T_N.' T_C.' T_WW.' T_Wsw.' T_s.' T_Ese.' T_EE.'],(b.K+1)*(b.J+2)+length(b.Z.T(b.Z.T(1,:)<=b.geom.z(n)))-1+[-2*(b.J+2) -(b.J+2) -2 -1 0 1 2],(J+2),(b.K+2)*(b.J+2));

end

rhs_T  = [rhs1; rhs; rhs2];
if n~=1
   th_e.T(1,:)=[];
   rhs_T(1,:)=[];

end
if n~=length(b.bc.z)/2
   th_e.T(end,:)=[];
   rhs_T(end,:)=[];

end

if o==1
    bc.th_e.T1 = [bc.th_e.T1; th_e.T];
    bc.rhs.T1  = [bc.rhs.T1; rhs_T];

else
    bc.th_e.Tend = [bc.th_e.Tend; th_e.T];
    bc.rhs.Tend  = [bc.rhs.Tend; rhs_T];

end