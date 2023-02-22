function[b]=block_grid(b, mesh, flowopt)
    ax=flowopt.ax;
    
    % segmentation of edges
        % r-direction
            b.r.u    = mesh.r.u(mesh.r.u>=b.geom.r(1) & mesh.r.u<=b.geom.r(end));
            %b.r.w    = [b.r.u(1) mesh.r.w(mesh.r.w>b.geom.r(1) & mesh.r.w<b.geom.r(end)) b.r.u(end)];
            b.r.w    = [b.r.u(1) (b.r.u(2:end)+b.r.u(1:end-1))/2 b.r.u(end)];
            b.dr.c   = diff(b.r.u);
            b.dr.rm  = diff(b.r.w);
            b.K=length(b.dr.c);
            
            % transformation
                b.delta_xi=1/b.K;
                
            %b.xi_r.u = mesh.xi_r.u(mesh.r.u>=b.geom.r(1) & mesh.r.u<=b.geom.r(end))*b.delta_xi;
            b.xi_r.u = diff(b.delta_xi*([0 0.5:(b.K-0.5) b.K]))./b.dr.rm;
            %b.xi_r.w = [b.xi_r.u(1)/b.delta_xi mesh.xi_r.w(mesh.r.w>b.geom.r(1) & mesh.r.w<b.geom.r(end)) b.xi_r.u(end)/b.delta_xi]*b.delta_xi;
            %b.xi_r.w = [b.xi_r.u(1) b.delta_xi./b.dr.c b.xi_r.u(end)];
            b.xi_r.w = [b.xi_r.u(1) (b.xi_r.u(2:end)+b.xi_r.u(1:end-1))/2 b.xi_r.u(end)];
            
        % z-direction
            b.z.w    = mesh.z.w(mesh.z.w>=b.geom.z(1) & mesh.z.w<=b.geom.z(end));
            %b.z.u    = [b.z.w(1) mesh.z.u(mesh.z.u>b.geom.z(1) & mesh.z.u<b.geom.z(end)) b.z.w(end)];
            b.z.u    = [b.z.w(1) (b.z.w(2:end)+b.z.w(1:end-1))/2 b.z.w(end)];
            b.dz.c   = diff(b.z.w);
            b.dz.zm  = diff(b.z.u);
            b.J=length(b.dz.c);
            
            % transformation
                b.delta_eta=1/b.J;
                
            %b.eta_z.w = mesh.eta_z.w(mesh.z.w>=b.geom.z(1) & mesh.z.w<=b.geom.z(end))*b.delta_eta;
            b.eta_z.w = diff(b.delta_eta*([0 0.5:(b.J-0.5) b.J]))./b.dz.zm;
            %b.eta_z.u = [b.eta_z.w(1)/b.delta_eta mesh.eta_z.u(mesh.z.u>b.geom.z(1) & mesh.z.u<b.geom.z(end)) b.eta_z.w(end)/b.delta_eta]*b.delta_eta;
            %b.eta_z.u = [b.eta_z.w(1) b.delta_eta./b.dz.c b.eta_z.w(end)];
            b.eta_z.u = [b.eta_z.w(1) (b.eta_z.w(2:end)+b.eta_z.w(1:end-1))/2 b.eta_z.w(end)];
            
        % xi-direction
            b.xi.u=b.delta_xi*(0:b.K);
            b.dxi.c=diff(b.xi.u);
            b.xi.w=b.delta_xi*([0 0.5:(b.K-0.5) b.K]);
            b.dxi.rm=diff(b.xi.w);

        % eta-direction
            b.eta.w=b.delta_eta*(0:b.J);
            b.deta.c=diff(b.eta.w);
            b.eta.u=b.delta_eta*([0 0.5:(b.J-0.5) b.J]);
            b.deta.zm=diff(b.eta.u);
    
        % radial shape
            b.R1.u=[];
            b.Rend.u=[];

            for n=1:length(b.geom.sh.z)/2
                J = length(b.z.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)));
                eval(['R1=' b.geom.sh.z{n} ';'])
                eval(['Rend=' b.geom.sh.z{length(b.geom.sh.z)/2+n} ';'])
                ph=0;
                if n==1
                    ph=ph+1;
                    b.R1.w=R1(b.geom.z(n));
                    b.Rend.w=Rend(b.geom.z(n));
                    
                end
                if n==length(b.geom.sh.z)/2
                    ph=ph+1;
                    
                end
                b.R1.u   = [b.R1.u R1(b.z.u(b.geom.z(n)<=b.z.u & b.geom.z(n+1)>=b.z.u)).*ones(1,J+ph)];
                b.R1.w   = [b.R1.w(1:end-1) [(b.R1.w(end)+R1(b.geom.z(n)))/2 R1(b.z.w(b.geom.z(n)<b.z.w & b.geom.z(n+1)>=b.z.w)).*ones(1,J)]];
                b.Rend.u = [b.Rend.u Rend(b.z.u(b.geom.z(n)<=b.z.u & b.geom.z(n+1)>=b.z.u)).*ones(1,J+ph)];
                b.Rend.w = [b.Rend.w(1:end-1) [(b.Rend.w(end)+Rend(b.geom.z(n)))/2 Rend(b.z.w(b.geom.z(n)<b.z.w & b.geom.z(n+1)>=b.z.w)).*ones(1,J)]];
                
            end
%             b.R1.w   = [b.R1.u(1) (b.R1.u(3:end-1)+b.R1.u(2:end-2))/2 b.R1.u(end)];
%             b.Rend.w = [b.Rend.u(1) (b.Rend.u(3:end-1)+b.Rend.u(2:end-2))/2 b.Rend.u(end)];
                
            %%{
            b.R1_z.w   = diff(b.R1.u)./b.dz.zm;
            b.R1_z.u   = [b.R1_z.w(1) diff(b.R1.w)./b.dz.c b.R1_z.w(end)];
            b.Rend_z.w = diff(b.Rend.u)./b.dz.zm;
            b.Rend_z.u = [b.Rend_z.w(1) diff(b.Rend.w)./b.dz.c b.Rend_z.w(end)];
            %}
            %{
            b.R1_z.w   = b.eta_z.w.*diff(b.R1.u)./b.deta.zm;
            b.R1_z.u   = [b.R1_z.w(1) b.eta_z.u(2:end-1).*diff(b.R1.w)./b.deta.c b.R1_z.w(end)];
            b.Rend_z.w = b.eta_z.w.*diff(b.Rend.u)./b.deta.zm;
            b.Rend_z.u = [b.Rend_z.w(1) b.eta_z.u(2:end-1).*diff(b.Rend.w)./b.deta.c b.Rend_z.w(end)];
            %}
            
            %%{
            b.R1_z_rlb1.w   = 1./b.dz.zm;
            b.R1_z_rlb1.u   = [b.R1_z_rlb1.w(1) 1./b.dz.c b.R1_z_rlb1.w(end)];
            b.Rend_z_rlbend.w = 1./b.dz.zm;
            b.Rend_z_rlbend.u = [b.Rend_z_rlbend.w(1) 1./b.dz.c b.Rend_z_rlbend.w(end)];
            %}
            %{
            b.R1_z_rlb1.w   = b.eta_z.w./b.deta.zm;
            b.R1_z_rlb1.u   = [b.R1_z_rlb1.w(1) b.eta_z.u(2:end-1)./b.deta.c b.R1_z_rlb1.w(end)];
            b.Rend_z_rlbend.w = b.eta_z.w./b.deta.zm;
            b.Rend_z_rlbend.u = [b.Rend_z_rlbend.w(1) b.eta_z.u(2:end-1)./b.deta.c b.Rend_z_rlbend.w(end)];
            %}
            
            %%{            
            b.R1_zz.w   = diff(b.R1_z.u)./b.dz.zm;
            b.R1_zz.u   = [b.R1_zz.w(1) diff(b.R1_z.w)./b.dz.c b.R1_zz.w(end)];
            b.Rend_zz.w = diff(b.Rend_z.u)./b.dz.zm;
            b.Rend_zz.u = [b.Rend_zz.w(1) diff(b.Rend_z.w)./b.dz.c b.Rend_zz.w(end)];
            %}
            %{
            b.R1_zz.w   = b.eta_z.w.*diff(b.R1_z.u)./b.deta.zm;
            b.R1_zz.u   = [b.R1_zz.w(1) b.eta_z.u(2:end-1).*diff(b.R1_z.w)./b.deta.c b.R1_zz.w(end)];
            b.Rend_zz.w = b.eta_z.w.*diff(b.Rend_z.u)./b.deta.zm;
            b.Rend_zz.u = [b.Rend_zz.w(1) b.eta_z.u(2:end-1).*diff(b.Rend_z.w)./b.deta.c b.Rend_zz.w(end)];
            %}
            
            b.R1.u   = ones(b.K+1,1)*b.R1.u;
            b.R1.w   = ones(b.K+2,1)*b.R1.w;
            b.R1.T   = [b.R1.u; b.R1.u(1,:)];
            b.R1.v   = b.R1.w(1:end-1,:);
            b.Rend.u = ones(b.K+1,1)*b.Rend.u;
            b.Rend.w = ones(b.K+2,1)*b.Rend.w;
            b.Rend.T = [b.Rend.u; b.Rend.u(1,:)];
            b.Rend.v = b.Rend.w(1:end-1,:);
            
            b.R1_z_rlb1.u     = ones(b.K+1,1)*b.R1_z_rlb1.u;
            b.R1_z_rlb1.w     = ones(b.K+2,1)*b.R1_z_rlb1.w;
            b.R1_z_rlb1.T     = [b.R1_z_rlb1.u; b.R1_z_rlb1.u(1,:)];
            b.R1_z_rlb1.v     = b.R1_z_rlb1.w(1:end-1,:);
            b.Rend_z_rlbend.u = ones(b.K+1,1)*b.Rend_z_rlbend.u;
            b.Rend_z_rlbend.w = ones(b.K+2,1)*b.Rend_z_rlbend.w;
            b.Rend_z_rlbend.T = [b.Rend_z_rlbend.u; b.Rend_z_rlbend.u(1,:)];
            b.Rend_z_rlbend.v = b.Rend_z_rlbend.w(1:end-1,:);
            
            b.R1_z.u   = ones(b.K+1,1)*b.R1_z.u;
            b.R1_z.w   = ones(b.K+2,1)*b.R1_z.w;
            b.R1_z.T   = [b.R1_z.u; b.R1_z.u(1,:)];
            b.R1_z.v   = b.R1_z.w(1:end-1,:);
            b.Rend_z.u   = ones(b.K+1,1)*b.Rend_z.u;
            b.Rend_z.w   = ones(b.K+2,1)*b.Rend_z.w;
            b.Rend_z.T   = [b.Rend_z.u; b.Rend_z.u(1,:)];
            b.Rend_z.v   = b.Rend_z.w(1:end-1,:);

            b.R1_zz.u   = ones(b.K+1,1)*b.R1_zz.u;
            b.R1_zz.w   = ones(b.K+2,1)*b.R1_zz.w;
            b.R1_zz.T   = [b.R1_zz.u; b.R1_zz.u(1,:)];
            b.R1_zz.v   = b.R1_zz.w(1:end-1,:);
            b.Rend_zz.u   = ones(b.K+1,1)*b.Rend_zz.u;
            b.Rend_zz.w   = ones(b.K+2,1)*b.Rend_zz.w;
            b.Rend_zz.T   = [b.Rend_zz.u; b.Rend_zz.u(1,:)];
            b.Rend_zz.v   = b.Rend_zz.w(1:end-1,:);
        
    % grid generation                
        % r,z-direction
            %[b.DZ.c,b.DR_cyl.c]       = meshgrid(b.dz.c, b.dr.c);
            %[b.DZ.rm,b.DR_cyl.rm]     = meshgrid(b.dz.c, b.dr.rm);
            %[b.DZ.zm,b.DR_cyl.zm]     = meshgrid(b.dz.zm, b.dr.c);
            [b.Z.u,b.R_cyl.u]         = meshgrid(b.z.u, b.r.u);
            [b.Z.w,b.R_cyl.w]         = meshgrid(b.z.w, b.r.w);
            [b.Z.p,b.R_cyl.p]         = meshgrid(b.z.u(2:end-1), b.r.w(2:end-1));
            [b.Z.T,b.R_cyl.T]         = meshgrid(b.z.u, b.r.w);
            [b.Z.v,b.R_cyl.v]         = meshgrid(b.z.w, b.r.u);
            
            % interpolation coefficients
                %{
                b.S.u=[ones(1,b.J+2); (b.R_cyl.T(3:end-1,:)-b.R_cyl.u(2:end-1,:))./(b.R_cyl.T(3:end-1,:)-b.R_cyl.T(2:end-2,:))];
                b.W.w=[ones(b.K+2,1) (b.Z.T(:,3:end-1)-b.Z.w(:,2:end-1))./(b.Z.T(:,3:end-1)-b.Z.T(:,2:end-2))];
                b.E.w=[(b.Z.w(:,2:end-1)-b.Z.T(:,2:end-2))./(b.Z.T(:,3:end-1)-b.Z.T(:,2:end-2)) ones(b.K+2,1)];
                b.N.u=[(b.R_cyl.u(2:end-1,:)-b.R_cyl.T(2:end-2,:))./(b.R_cyl.T(3:end-1,:)-b.R_cyl.T(2:end-2,:)); ones(1,b.J+2)];
                %}
                %%{
                b.S.u=[ones(1,b.J+2); 0.5*ones(b.K-1,b.J+2)];
                b.W.w=[ones(b.K+2,1) 0.5*ones(b.K+2,b.J-1)];
                b.E.w=[0.5*ones(b.K+2,b.J-1) ones(b.K+2,1)];
                b.N.u=[0.5*ones(b.K-1,b.J+2); ones(1,b.J+2)];
                %}
                
                b.SW.v=b.S.u(:,2:end-1).*b.W.w(2:end-1,:);
                b.SE.v=b.S.u(:,2:end-1).*b.E.w(2:end-1,:);
                b.NW.v=b.N.u(:,2:end-1).*b.W.w(2:end-1,:);
                b.NE.v=b.N.u(:,2:end-1).*b.E.w(2:end-1,:);
       
        % xi,eta-direction
            [b.DETA.c,b.DXI.c]        = meshgrid(b.deta.c, b.dxi.c);
            [b.DETA.rm,b.DXI.rm]      = meshgrid(b.deta.c, b.dxi.rm);
            [b.DETA.zm,b.DXI.zm]      = meshgrid(b.deta.zm, b.dxi.c);
            [b.ETA.u,b.XI.u]          = meshgrid(b.eta.u, b.xi.u);
            [b.ETA.w,b.XI.w]          = meshgrid(b.eta.w, b.xi.w);
            [b.ETA.p,b.XI.p]          = meshgrid(b.eta.u(2:end-1), b.xi.w(2:end-1));
            [b.ETA.T,b.XI.T]          = meshgrid(b.eta.u, b.xi.w);
            [b.ETA.v,b.XI.v]          = meshgrid(b.eta.w, b.xi.u);

            % del(ETA)/del(z), del(XI_cyl)/del(r)
                [b.ETA_z.u,b.XI_cyl_r.u]  = meshgrid(b.eta_z.u, b.xi_r.u);
                [b.ETA_z.w,b.XI_cyl_r.w]  = meshgrid(b.eta_z.w, b.xi_r.w);
                [b.ETA_z.p,b.XI_cyl_r.p]  = meshgrid(b.eta_z.u(2:end-1), b.xi_r.w(2:end-1));
                [b.ETA_z.T,b.XI_cyl_r.T]  = meshgrid(b.eta_z.u, b.xi_r.w);
                [b.ETA_z.v,b.XI_cyl_r.v]  = meshgrid(b.eta_z.w, b.xi_r.u);
                
            R           = @(R1, Rend, R_cyl) R1+(R_cyl-b.geom.r(1)).*(Rend-R1)/(b.geom.r(end)-b.geom.r(1));
            R_z         = @(R1_z, Rend_z, R_cyl) R1_z+(R_cyl-b.geom.r(1)).*(Rend_z-R1_z)/(b.geom.r(end)-b.geom.r(1));
            %R_zz        = @(R1_zz, Rend_zz, R_cyl) R1_zz+(R_cyl-b.geom.r(1)).*(Rend_zz-R1_zz)/(b.geom.r(end)-b.geom.r(1));
            F           = @(R1, Rend) (b.geom.r(end)-b.geom.r(1))./(Rend-R1);
            R_cyl_z     = @(R1, Rend, R1_z, Rend_z, R, XI_r) (XI_r.*((-R1_z).*(Rend-R1)-(Rend_z-R1_z).*(R-R1))./(Rend-R1));
    
            % R
                b.R.u=R(b.R1.u, b.Rend.u, b.R_cyl.u);
                b.R.w=R(b.R1.w, b.Rend.w, b.R_cyl.w);
                b.R.T=R(b.R1.T, b.Rend.T, b.R_cyl.T);
                b.R.p=b.R.T(2:end-1,2:end-1);
                b.R.v=R(b.R1.v, b.Rend.v, b.R_cyl.v);
            
            %%{
            % R_z
                b.R_z.u=R_z(b.R1_z.u, b.Rend_z.u, b.R_cyl.u);
                b.R_z.w=R_z(b.R1_z.w, b.Rend_z.w, b.R_cyl.w);
                b.R_z.T=R_z(b.R1_z.T, b.Rend_z.T, b.R_cyl.T);
                b.R_z.p=b.R_z.T(2:end-1,2:end-1);
                b.R_z.v=R_z(b.R1_z.v, b.Rend_z.v, b.R_cyl.v);
            %}
            %{
            % R_zz
                b.R_zz.u=R_zz(b.R1_zz.u, b.Rend_zz.u, b.R_cyl.u);
                b.R_zz.w=R_zz(b.R1_zz.w, b.Rend_zz.w, b.R_cyl.w);
                b.R_zz.T=R_zz(b.R1_zz.T, b.Rend_zz.T, b.R_cyl.T);
                b.R_zz.p=b.R_zz.T(2:end-1,2:end-1);
                b.R_zz.v=R_zz(b.R1_zz.v, b.Rend_zz.v, b.R_cyl.v);
            %}

            % del(XI)/del(r)
                b.XI_r.u   = b.XI_cyl_r.u.*F(b.R1.u, b.Rend.u);
                b.XI_r.w   = b.XI_cyl_r.w.*F(b.R1.w, b.Rend.w);
                b.XI_r.T   = b.XI_cyl_r.T.*F(b.R1.T, b.Rend.T);
                b.XI_r.p   = b.XI_r.T(2:end-1,2:end-1);
                b.XI_r.v   = b.XI_cyl_r.v.*F(b.R1.v, b.Rend.v);

            % Jacobi determinant
                b.JA.u=1./(b.XI_r.u.*b.ETA_z.u);
                b.JA.w=1./(b.XI_r.w.*b.ETA_z.w);
                b.JA.T=1./(b.XI_r.T.*b.ETA_z.T);
                b.JA.p=1./(b.XI_r.p.*b.ETA_z.p);
                b.JA.v=1./(b.XI_r.v.*b.ETA_z.v);
                
            % del(XI)/del(z)
                %%{
                b.XI_z.u   = R_cyl_z(b.R1.u, b.Rend.u, b.R1_z.u, b.Rend_z.u, b.R.u, b.XI_r.u);
                b.XI_z.w   = R_cyl_z(b.R1.w, b.Rend.w, b.R1_z.w, b.Rend_z.w, b.R.w, b.XI_r.w);
                b.XI_z.T   = R_cyl_z(b.R1.T, b.Rend.T, b.R1_z.T, b.Rend_z.T, b.R.T, b.XI_r.T);
                b.XI_z.p   = b.XI_z.T(2:end-1,2:end-1);
                b.XI_z.v   = R_cyl_z(b.R1.v, b.Rend.v, b.R1_z.v, b.Rend_z.v, b.R.v, b.XI_r.v);
                %}
                %{
                b.XI_z.w   = -(b.R.T(:,2:end)-b.R.T(:,1:end-1))./(b.ETA.T(:,2:end)-b.ETA.T(:,1:end-1))./b.JA.w;
                b.XI_z.T   = -[b.XI_z.w(:,1) (b.R.w(:,2:end)-b.R.w(:,1:end-1))./(b.ETA.w(:,2:end)-b.ETA.w(:,1:end-1))./b.JA.T(:,2:end-1) b.XI_z.w(:,end)];
                b.XI_z.p   = -b.XI_z.T(2:end-1,2:end-1);
                b.XI_z.v   = -(b.R.u(:,2:end)-b.R.u(:,1:end-1))./(b.ETA.u(:,2:end)-b.ETA.u(:,1:end-1))./b.JA.v;
                b.XI_z.u   = -[b.XI_z.v(:,1) (b.R.v(:,2:end)-b.R.v(:,1:end-1))./(b.ETA.v(:,2:end)-b.ETA.v(:,1:end-1))./b.JA.u(:,2:end-1) b.XI_z.v(:,end)];
                %}
                
            if max(strncmp('ss', {b.bc.z{1:end/2}},2))==1
                F_rlb_C       = @(R1, Rend, R1_z, Rend_z, R_rlb1, XI_r_rlb1, XI_z) XI_r_rlb1.*(R1_z.*R_rlb1 - Rend_z.*(R_rlb1-1))+2*XI_z./(Rend-R1);
                F_rlb_E       = @(R, Rend, R1_z_rlb1, XI_r_rlb1) XI_r_rlb1.*R1_z_rlb1.*(R-Rend);
                F_inv_rlb   = @(R1, Rend) (b.geom.r(end)-b.geom.r(1))./((Rend-R1).^2);

                % del(XI_r^-1)/del(r_lb)    
                    b.XI_r_inv_rlb1.u           = -1./b.XI_cyl_r.u/(b.geom.r(end)-b.geom.r(1));
                    b.XI_r_inv_rlb1.w           = -1./b.XI_cyl_r.w/(b.geom.r(end)-b.geom.r(1));
                    b.XI_r_inv_rlb1.T           = -1./b.XI_cyl_r.T/(b.geom.r(end)-b.geom.r(1));
                    b.XI_r_inv_rlb1.p           = b.XI_r_inv_rlb1.T(2:end-1,2:end-1);
                    b.XI_r_inv_rlb1.v           = -1./b.XI_cyl_r.v/(b.geom.r(end)-b.geom.r(1));

                % del(R)/del(r_lb)
                    b.R_rlb1.u                  = 1-(b.R_cyl.u-b.geom.r(1))/(b.geom.r(end)-b.geom.r(1));
                    b.R_rlb1.w                  = 1-(b.R_cyl.w-b.geom.r(1))/(b.geom.r(end)-b.geom.r(1));
                    b.R_rlb1.T                  = 1-(b.R_cyl.T-b.geom.r(1))/(b.geom.r(end)-b.geom.r(1));
                    b.R_rlb1.p                  = b.R_rlb1.T(2:end-1,2:end-1);
                    b.R_rlb1.v                  = 1-(b.R_cyl.v-b.geom.r(1))/(b.geom.r(end)-b.geom.r(1));

                % del(JA)/del(r_lb)
                    b.JA_rlb1.u                 = b.XI_r_inv_rlb1.u./b.ETA_z.u;
                    b.JA_rlb1.w                 = b.XI_r_inv_rlb1.w./b.ETA_z.w;
                    b.JA_rlb1.T                 = b.XI_r_inv_rlb1.T./b.ETA_z.T;
                    b.JA_rlb1.p                 = b.XI_r_inv_rlb1.p./b.ETA_z.p;
                    b.JA_rlb1.v                 = b.XI_r_inv_rlb1.v./b.ETA_z.v;

                % del(XI_r)/del(r_lb)
                    b.XI_r_rlb1.u   = b.XI_cyl_r.u.*F_inv_rlb(b.R1.u, b.Rend.u);
                    b.XI_r_rlb1.w   = b.XI_cyl_r.w.*F_inv_rlb(b.R1.w, b.Rend.w);
                    b.XI_r_rlb1.T   = b.XI_cyl_r.T.*F_inv_rlb(b.R1.T, b.Rend.T);
                    b.XI_r_rlb1.p   = b.XI_r_rlb1.T(2:end-1,2:end-1);
                    b.XI_r_rlb1.v   = b.XI_cyl_r.v.*F_inv_rlb(b.R1.v, b.Rend.v);

                % del(XI_z)/del(r_lb)
                    b.XI_z_rlb1_C.u   = F_rlb_C(b.R1.u, b.Rend.u, b.R1_z.u, b.Rend_z.u, b.R_rlb1.u, b.XI_r_rlb1.u, b.XI_z.u);
                    b.XI_z_rlb1_C.w   = F_rlb_C(b.R1.w, b.Rend.w, b.R1_z.w, b.Rend_z.w, b.R_rlb1.w, b.XI_r_rlb1.w, b.XI_z.w);
                    b.XI_z_rlb1_C.T   = F_rlb_C(b.R1.T, b.Rend.T, b.R1_z.T, b.Rend_z.T, b.R_rlb1.T, b.XI_r_rlb1.T, b.XI_z.T);
                    b.XI_z_rlb1_C.p   = b.XI_z_rlb1_C.T(2:end-1,2:end-1);
                    b.XI_z_rlb1_C.v   = F_rlb_C(b.R1.v, b.Rend.v, b.R1_z.v, b.Rend_z.v, b.R_rlb1.v, b.XI_r_rlb1.v, b.XI_z.v);

                    b.XI_z_rlb1_E.u   = F_rlb_E(b.R.u, b.Rend.u, b.R1_z_rlb1.u, b.XI_r_rlb1.u);
                    b.XI_z_rlb1_E.w   = F_rlb_E(b.R.w, b.Rend.w, b.R1_z_rlb1.w, b.XI_r_rlb1.w);
                    b.XI_z_rlb1_E.T   = F_rlb_E(b.R.T, b.Rend.T, b.R1_z_rlb1.T, b.XI_r_rlb1.T);
                    b.XI_z_rlb1_E.p   = b.XI_z_rlb1_E.T(2:end-1,2:end-1);
                    b.XI_z_rlb1_E.v   = F_rlb_E(b.R.v, b.Rend.v, b.R1_z_rlb1.v, b.XI_r_rlb1.v);

            end
            if max(strncmp('ss', {b.bc.z{end/2+1:end}},2))==1
                F_rlb_C       = @(R1, Rend, R1_z, Rend_z, R_rlbend, XI_r_rlbend, XI_z) XI_r_rlbend.*(Rend_z.*R_rlbend - R1_z.*(R_rlbend-1))-2*XI_z./(Rend-R1);
                F_rlb_E       = @(R, R1, Rend_z_rlbend, XI_r_rlbend) XI_r_rlbend.*Rend_z_rlbend.*(R-R1);
                F_inv_rlb   = @(R1, Rend) -(b.geom.r(end)-b.geom.r(1))./((Rend-R1).^2);

                % del(XI_r^-1)/del(r_lb)
                    b.XI_r_inv_rlbend.u           = 1./b.XI_cyl_r.u/(b.geom.r(end)-b.geom.r(1));
                    b.XI_r_inv_rlbend.w           = 1./b.XI_cyl_r.w/(b.geom.r(end)-b.geom.r(1));
                    b.XI_r_inv_rlbend.T           = 1./b.XI_cyl_r.T/(b.geom.r(end)-b.geom.r(1));
                    b.XI_r_inv_rlbend.p           = b.XI_r_inv_rlbend.T(2:end-1,2:end-1);
                    b.XI_r_inv_rlbend.v           = 1./b.XI_cyl_r.v/(b.geom.r(end)-b.geom.r(1));

                % del(R)/del(r_lb)
                    b.R_rlbend.u                  = (b.R_cyl.u-b.geom.r(1))/(b.geom.r(end)-b.geom.r(1));
                    b.R_rlbend.w                  = (b.R_cyl.w-b.geom.r(1))/(b.geom.r(end)-b.geom.r(1));
                    b.R_rlbend.T                  = (b.R_cyl.T-b.geom.r(1))/(b.geom.r(end)-b.geom.r(1));
                    b.R_rlbend.p                  = b.R_rlbend.T(2:end-1,2:end-1);
                    b.R_rlbend.v                  = (b.R_cyl.v-b.geom.r(1))/(b.geom.r(end)-b.geom.r(1));

                % del(JA)/del(r_lb)
                    b.JA_rlbend.u                 = b.XI_r_inv_rlbend.u./b.ETA_z.u;
                    b.JA_rlbend.w                 = b.XI_r_inv_rlbend.w./b.ETA_z.w;
                    b.JA_rlbend.T                 = b.XI_r_inv_rlbend.T./b.ETA_z.T;
                    b.JA_rlbend.p                 = b.XI_r_inv_rlbend.p./b.ETA_z.p;
                    b.JA_rlbend.v                 = b.XI_r_inv_rlbend.v./b.ETA_z.v;

                % del(XI_r)/del(r_lb)
                    b.XI_r_rlbend.u   = b.XI_cyl_r.u.*F_inv_rlb(b.R1.u, b.Rend.u);
                    b.XI_r_rlbend.w   = b.XI_cyl_r.w.*F_inv_rlb(b.R1.w, b.Rend.w);
                    b.XI_r_rlbend.T   = b.XI_cyl_r.T.*F_inv_rlb(b.R1.T, b.Rend.T);
                    b.XI_r_rlbend.p   = b.XI_r_rlbend.T(2:end-1,2:end-1);
                    b.XI_r_rlbend.v   = b.XI_cyl_r.v.*F_inv_rlb(b.R1.v, b.Rend.v);

                % del(XI_z)/del(r_lb)
                    b.XI_z_rlbend_C.u   = F_rlb_C(b.R1.u, b.Rend.u, b.R1_z.u, b.Rend_z.u, b.R_rlbend.u, b.XI_r_rlbend.u, b.XI_z.u);
                    b.XI_z_rlbend_C.w   = F_rlb_C(b.R1.w, b.Rend.w, b.R1_z.w, b.Rend_z.w, b.R_rlbend.w, b.XI_r_rlbend.w, b.XI_z.w);
                    b.XI_z_rlbend_C.T   = F_rlb_C(b.R1.T, b.Rend.T, b.R1_z.T, b.Rend_z.T, b.R_rlbend.T, b.XI_r_rlbend.T, b.XI_z.T);
                    b.XI_z_rlbend_C.p   = b.XI_z_rlbend_C.T(2:end-1,2:end-1);
                    b.XI_z_rlbend_C.v   = F_rlb_C(b.R1.v, b.Rend.v, b.R1_z.v, b.Rend_z.v, b.R_rlbend.v, b.XI_r_rlbend.v, b.XI_z.v);

                    b.XI_z_rlbend_E.u   = F_rlb_E(b.R.u, b.R1.u, b.Rend_z_rlbend.u, b.XI_r_rlbend.u);
                    b.XI_z_rlbend_E.w   = F_rlb_E(b.R.w, b.R1.w, b.Rend_z_rlbend.w, b.XI_r_rlbend.w);
                    b.XI_z_rlbend_E.T   = F_rlb_E(b.R.T, b.R1.T, b.Rend_z_rlbend.T, b.XI_r_rlbend.T);
                    b.XI_z_rlbend_E.p   = b.XI_z_rlbend_E.T(2:end-1,2:end-1);
                    b.XI_z_rlbend_E.v   = F_rlb_E(b.R.v, b.R1.v, b.Rend_z_rlbend.v, b.XI_r_rlbend.v);

            end
            
    if isfield(b, 'V_r')
        b.V_l=b.V_r*((pi*b.R.u(end,2:end-1)).^ax.*b.R.u(end,2:end-1) - (pi*b.R.u(1,2:end-1)).^ax.*b.R.u(1,2:end-1))*b.dz.c.';
        
    end