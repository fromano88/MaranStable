if o==1
    % lower boundary
        index_bc_z=0;
        index_bc_r=1;
        K=1;
        R_z.u=b.R1_z.u;
        R_z.T=b.R1_z.T;

elseif o==2
    % upper boundary
        index_bc_z=length(b.bc.z)/2;
        index_bc_r=length(b.bc.r)/2;
        K=b.K+2;
        R_z.u=b.Rend_z.u;
        R_z.T=b.Rend_z.T;

end

if max(strncmp('s', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},1))
    N.u=sqrt(1+R_z.u(1,:).^2);

    % normal heat transfer
        F.T=lambda(b.T.T(K,:)).*N.u.^2.*b.XI_r.T(K,:)./[b.DXI.c(K-2*(o-1),1) b.DXI.c(K-2*(o-1),:) b.DXI.c(K-2*(o-1),end)].*2;
        F.T=-(-1)^o*F.T(1,2:end-1).*b.dz.c;
        F1.T=-lambda(b.T.T(K,:)).*R_z.T(1,:).*b.ETA_z.T(K,:)./[b.DETA.zm(K-2*(o-1),1) b.DETA.c(K-2*(o-1),:) b.DETA.zm(K-2*(o-1),end)];
        F1.T=F1.T(1,2:end-1).*b.dz.c;

        T_Wsw = [0 -b.W.w(1,:).*F1.T 0];
        T_Ese = [0 b.E.w(1,:).*F1.T 0];
        T_s   = [0 (-F.T - [0 b.E.w(1,1:end-1)].*F1.T + [b.W.w(1,2:end) 0].*F1.T)  0];
        T_C   = [0 F.T 0];

    % Jacobi of normal heat transfer with respect to T
        Ja_T_s=dlambda(b.T.T(K,:)).*N.u.^2.*b.XI_r.T(K,:).*(b.T.T(K+(2-o),:)-b.T.T(K+(1-o),:))./[b.DXI.c(K-2*(o-1),1) b.DXI.c(K-2*(o-1),:) b.DXI.c(K-2*(o-1),end)].*2-dlambda(b.T.T(K,:)).*R_z.T(1,:).*b.ETA_z.T(K,:).*[(b.T.T(K,2)-b.T.T(K,1)) (b.T.w(K,2:end)-b.T.w(K,1:end-1)) (b.T.T(K,end)-b.T.T(K,end-1))]./[b.DETA.zm(K-2*(o-1),1) b.DETA.c(K-2*(o-1),:) b.DETA.zm(K-2*(o-1),end)];
        Ja_T_s=Ja_T_s(1,2:end-1).*b.dz.c;

    if stability~=1 && max(strncmp('ss', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},2))
        for pp=1:2
            if max(strncmp('ss', {b.bc.z{(pp-1)*length(b.bc.z)/2+1:(pp-1)*length(b.bc.z)/2+length(b.bc.z)/2}},2))
                if pp==1
                    % derivative with respect to lower boundary
                        R_z_rlb.T=b.R1_z_rlb1.T;
                        XI_r_rlb.T=b.XI_r_rlb1.T;
                        XI_z_rlb_C.T=b.XI_z_rlb1_C.T;
                        XI_z_rlb_E.T=b.XI_z_rlb1_E.T;
                        XI_z_rlb_E.u=b.XI_z_rlb1_E.u;

                elseif pp==2
                    % derivative with respect to upper boundary
                        R_z_rlb.T=b.Rend_z_rlbend.T;
                        XI_r_rlb.T=b.XI_r_rlbend.T;
                        XI_z_rlb_C.T=b.XI_z_rlbend_C.T;
                        XI_z_rlb_E.T=b.XI_z_rlbend_E.T;
                        XI_z_rlb_E.u=b.XI_z_rlbend_E.u;

                end

            % derivative of normal heat transfer with respect to Rs
                if o==pp
                    N_inv_rlb.u=-0.5./N.u.^3.*2.*R_z.u(1,:).*R_z_rlb.T(1,:);
                    N_rlb.u=0.5./N.u.*2.*R_z.u(1,:).*R_z_rlb.T(1,:);

                    % lambda*N^2*dxi_r/dRs*dT/delta.xi
                        F.u=(b.T.T(K+(2-o),:)-b.T.T(K+(1-o),:))./[b.DXI.c(K-2*(o-1),1) b.DXI.c(K-2*(o-1),:) b.DXI.c(K-2*(o-1),end)].*2;
                        FC.u=lambda(b.T.T(K,:)).*N.u.^2.*XI_r_rlb.T(K,:).*F.u;

                    % lambda*(-dr_z/dRs*eta_z*dT/delta.eta)
                        F1.u=[(b.T.T(K,2)-b.T.T(K,1)) (b.T.w(K,2:end)-b.T.w(K,1:end-1)) (b.T.T(K,end)-b.T.T(K,end-1))]./[b.DETA.zm(K-2*(o-1),1) b.DETA.c(K-2*(o-1),:) b.DETA.zm(K-2*(o-1),end)];
                        FE.u=lambda(b.T.T(K,:)).*(-R_z_rlb.T(1,:).*b.ETA_z.T(K,:).*F1.u);
                        
                    % lambda*(d(N^2)/dRs*xi_r*dT/delta.xi)
                        FE.u=FE.u + lambda(b.T.T(K,:)).*2.*N.u.*N_rlb.u.*b.XI_r.T(K,:).*F.u;

                        Rs_W=-[0 b.W.w(1,2:end)].*FE.u(2:end-1).*b.dz.c;
                        Rs_E=[b.E.w(1,1:end-1) 0].*FE.u(2:end-1).*b.dz.c;
                        Rs_ad_W=- [0 b.E.w(1,1:end-1)].*FE.u(2:end-1).*b.dz.c;
                        Rs_ad_E=[b.W.w(1,2:end) 0].*FE.u(2:end-1).*b.dz.c;
                        Rs_C=FC.u(2:end-1).*b.dz.c + Rs_ad_W + Rs_ad_E;

                    % derivative of adiabatic boundary with respect to Rs
                        % (-dr_z/dRs*eta_z*dT/delta.eta)
                            Ja.ad.Rs_W=-0.5.*(-R_z_rlb.T(1,:).*b.ETA_z.T(K,:).*F1.u).*[b.dz.zm(1) b.dz.c b.dz.zm(end)];

                        % 1*(d(N^2)/dRs*xi_r*dT/delta.xi)
                            Ja.ad.Rs_W=Ja.ad.Rs_W-0.5.*2.*N.u.*N_rlb.u.*b.XI_r.T(K,:).*F.u.*[b.dz.zm(1) b.dz.c b.dz.zm(end)];
                            Ja.ad.Rs_E=-Ja.ad.Rs_W;

                    % Jacobi of normal heat transfer with respect to Rs
                        Ja.nht.Rs=sparse(b.J+2,b.J+2);
                        Ja.nht.Rs(2:end-1,:)=spdiags([Rs_W.' Rs_C.' Rs_E.'],[0 1 2],b.J,b.J+2);

                        if (strcmp(b.bc.r{index_bc_r}(1),'w')==1 || strcmp(b.bc.r{index_bc_r}(1),'i')==1 || strcmp(b.bc.r{index_bc_r}(1),'o')==1) && strcmp(b.bc.r{index_bc_r}(2),'a')==1
                            Ja.nht.Rs(1,1)=XI_z_rlb_E.u(K+(1-o),1).*(b.T.T(K+(2-o),1)-b.T.T(K+(1-o),1))./b.DXI.rm(K+(1-o),1);

                        end
                        if (strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'w')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'i')==1 || strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(1),'o')==1) && strcmp(b.bc.r{length(b.bc.r)/2+index_bc_r}(2),'a')==1
                            Ja.nht.Rs(end,end)=XI_z_rlb_E.u(K+(1-o),end).*(b.T.T(K+(2-o),end)-b.T.T(K+(1-o),end))./b.DXI.rm(K+(1-o),end);

                        end

                    if o==1 && o==pp
                        % Jacobi of normal heat transfer with respect to Rs
                            Ja.nht.Rs1=Ja.nht.Rs;
                            Ja.nht.Rs1_ad_w=Rs_ad_W;
                            Ja.nht.Rs1_ad_e=Rs_ad_E;
                            Ja.ad.Rs1_E=Ja.ad.Rs_E;
                            Ja.ad.Rs1_W=Ja.ad.Rs_W;

                    elseif o==2 && o==pp
                        % Jacobi of normal heat transfer with respect to Rs
                            Ja.nht.Rsend=Ja.nht.Rs;
                            Ja.nht.Rsend_ad_w=Rs_ad_W;
                            Ja.nht.Rsend_ad_e=Rs_ad_E;
                            Ja.ad.Rsend_E=Ja.ad.Rs_E;
                            Ja.ad.Rsend_W=Ja.ad.Rs_W;

                    end

                end

            end

        end

    end

end

if o==1 && max(strncmp('s', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},1))
    b.nht1.T = spdiags([T_Wsw.' T_s.' T_Ese.' T_C.'],[-1 0 1 b.J+2],(b.J+2),(b.K+2)*(b.J+2));                
    % Jacobi of normal heat transfer with respect to T
        Ja.nht.T1=sparse(b.J+2,(b.K+2)*(b.J+2));
        Ja.nht.T1(2:end-1,:) = b.nht1.T(2:end-1,:) + spdiags(Ja_T_s.',1,b.J,(b.K+2)*(b.J+2));

    % equal temperature
        b.e1.T=speye(b.J+2,(b.K+2)*(b.J+2));

elseif o==2 && max(strncmp('s', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},1))
    index_surface = find(strncmp('s', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},1)==1);
    if strcmp(b.bc.z{index_bc_z+index_surface}(3),'c') == 1 && str2double(b.bc.z{o}(4)) == 0  % SFM
        l_lb  = diff(b.geom.z);
        r_i   = b.r.u(end);
        r_c   = b.r.u(1);
        T_d1l = str2func(b.bc.rhs.T.r{1}); T_d1l = T_d1l(r_i);
        T_d2l = str2func(b.bc.rhs.T.r{2}); T_d2l = T_d2l(r_i);
        T_ref = min(T_d1l,T_d2l);
        
        lambda = str2func(b.lambda);
        if strcmp(b.bc.rhs.T.z{index_bc_z+index_surface}(1),'@') % contant Biot number
            Bi = str2func(b.bc.rhs.T.z{index_bc_z+index_surface});
            h = Bi(b.Z.T(end,2:end-1)).*lambda(b.T.T(end,2:end-1))/l_lb; % heat transfer coefficient
        else
            subfolder = b.bc.rhs.T.z{index_bc_z+index_surface};

            Gamma  = l_lb/r_i;
            T_0    = (T_d1l+T_d2l)/2;
            rho    = str2func(b.rho);
            mu     = str2func(b.mu);
            cp     = str2func(b.cp);
            dsigma = str2func(b.dsigma);

            Re = abs(dsigma(T_0)*(T_d1l-T_d2l))*l_lb*rho(T_0)/mu(T_0)^2;
            Pr = mu(T_0)*cp(T_0)/lambda(T_0);
            Bi = Biot_func(Pr,flowopt.g,Gamma,Re,subfolder);
            h = Bi(0.5-b.Z.T(end,2:end-1)/l_lb).*lambda(b.T.T(end,2:end-1))/l_lb;
            
            % smoothing
            sm.bound  = 0.49;
            sm.z_orig = (l_lb/2-(b.Z.T(end,2:end-1)))/l_lb;
            sm.z      = sm.z_orig(sm.z_orig<=sm.bound & sm.z_orig>=-sm.bound);
            sm.h      = smooth(sm.z,h((sm.z_orig<=sm.bound & sm.z_orig>-sm.bound)),5,'moving');
            h = [h(sm.z_orig>sm.bound) sm.h' h(sm.z_orig<-sm.bound)];

            % constant Bi function
            %{
            %Bi = load(['BiotFunction_coeffs/' subfolder '/Bi_avg_Pr28.84_Gamma' num2str(Gamma,'%3.2f') '_g0.dat']);
            %Bi = load(['BiotFunction_coeffs/' subfolder '/Bi_avg_lin_Pr28.84_Gamma' num2str(Gamma,'%3.2f') '_g0.dat']);
            %Bi_avg = interp1(Bi(:,1),Bi(:,2),Re,'linear');
            Bi_avg = BiotAvgLin(Pr,flowopt.g,Gamma,Re,subfolder);
            h = Bi_avg*ones(size(b.Z.T(end,2:end-1))).*lambda(b.T.T(end,2:end-1))/l_lb;
            %}
            
            % linear Bi function
            %{
            Bi = load(['BiotFunction_coeffs/' subfolder '/Bi_lin_Pr28.84_Gamma' num2str(Gamma,'%3.2f') '_g0.dat']);            
            k = interp1(Bi(:,1),Bi(:,3),Re,'linear');
            d = interp1(Bi(:,1),Bi(:,2),Re,'linear');
            Bi_lin = @(z)k*z+d;
            h = Bi_lin(0.5-b.Z.T(end,2:end-1)/l_lb).*lambda(b.T.T(end,2:end-1))/l_lb;
            %}
            
        end
        T_s(2:end-1) = T_s(2:end-1)+h.*b.dz.c; % uncomment this line if it's not newton coolings law but arbitrary function
    end
    b.nhtend.T = spdiags([T_C.' T_Wsw.' T_s.' T_Ese.'],(b.K+1)*(b.J+2)+[-(b.J+2) -1 0 1],(b.J+2),(b.K+2)*(b.J+2));
    % Jacobi of normal heat transfer with respect to T
        Ja.nht.Tend=sparse(b.J+2,(b.K+2)*(b.J+2));
        Ja.nht.Tend(2:end-1,:) =  b.nhtend.T(2:end-1,:) + spdiags(Ja_T_s.',(b.K+1)*(b.J+2)+1,b.J,(b.K+2)*(b.J+2));

    % equal temperature
        b.eend.T=[sparse(b.J+2,(b.K+1)*(b.J+2)) speye(b.J+2)];

end