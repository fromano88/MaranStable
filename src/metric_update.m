function[b, oldb]=metric_update(b, x, Rs, old, oldb)
    
    for n=1:length(b.bc.z)/2
        if strcmp(b.bc.z{n}(1),'s')==1
            F_rlb_C     = @(R1, Rend, R1_z, Rend_z, R_rlb1, XI_r_rlb1, XI_z) XI_r_rlb1.*(R1_z.*R_rlb1 - Rend_z.*(R_rlb1-1))+2*XI_z./(Rend-R1);
            F_rlb_E     = @(R, Rend, R1_z_rlb1, XI_r_rlb1) XI_r_rlb1.*R1_z_rlb1.*(R-Rend);
            F_inv_rlb   = @(R1, Rend) (b.geom.r(end)-b.geom.r(1))./((Rend-R1).^2);

            if isfield(oldb,'Z')==0 || isequal(b.Z,oldb.Z)
                eval(['Rs1=x(Rs.' b.bc.z{n} '(1):Rs.' b.bc.z{n} '(2));'])
                
            else
                eval(['Rs1=interp1(oldb.ETA.u(1,(oldb.geom.z(n)<oldb.Z.u(1,:) & oldb.geom.z(n+1)>oldb.Z.u(1,:))).'',old.x(oldb.Rs.' b.bc.z{n} '(1):oldb.Rs.' b.bc.z{n} '(2)),b.ETA.u(1,(b.geom.z(n)<b.Z.u(1,:) & b.geom.z(n+1)>b.Z.u(1,:))).'',''spline'');'])
                
            end
            if eval(['b.' eval(['b.bc.z{' num2str(n) '}']) '.o==1'])
                b.Rs1=[b.Rs1; Rs1];
                
            end
                        
            b.R1.u(:,(b.geom.z(n)<b.z.u & b.geom.z(n+1)>b.z.u))=ones(b.K+1,1)*Rs1.';
            b.R1.w(:,(b.geom.z(n)<b.z.w & b.geom.z(n+1)>b.z.w))=ones(b.K+2,1)*(b.W.w(1,(b.geom.z(n)<b.z.w & b.geom.z(n+1)>b.z.w)).*Rs1(1:end-1).'+b.E.w(1,find(b.geom.z(n)<b.z.w & b.geom.z(n+1)>b.z.w)-1).*Rs1(2:end).');
            b.R1.T=[b.R1.u; b.R1.u(1,:)];
            b.R1.v=b.R1.w(1:end-1,:);
            %%{
            b.R1_z.w   = ones(b.K+2,1)*(diff(b.R1.u(1,:))./b.dz.zm);
            b.R1_z.u   = ones(b.K+1,1)*[b.R1_z.w(1) diff(b.R1.w(1,:))./b.dz.c b.R1_z.w(end)];
            b.R1_z.T   = [b.R1_z.u; b.R1_z.u(1,:)];
            b.R1_z.v   = b.R1_z.w(1:end-1,:);

            b.R1_zz.w   = ones(b.K+2,1)*(diff(b.R1_z.u(1,:))./b.dz.zm);
            b.R1_zz.u   = ones(b.K+1,1)*[b.R1_zz.w(1) diff(b.R1_z.w(1,:))./b.dz.c b.R1_zz.w(end)];
            b.R1_zz.T   = [b.R1_zz.u; b.R1_zz.u(1,:)];
            b.R1_zz.v   = b.R1_zz.w(1:end-1,:);
            %}
            
            %{
            b.R1_z.w   = ones(b.K+2,1)*(b.eta_z.w.*diff(b.R1.u(1,:))./b.deta.zm);
            b.R1_z.u   = ones(b.K+1,1)*[b.R1_z.w(1) b.eta_z.u(2:end-1).*diff(b.R1.w(1,:))./b.deta.c b.R1_z.w(end)];
            b.R1_z.T   = [b.R1_z.u; b.R1_z.u(1,:)];
            b.R1_z.v   = b.R1_z.w(1:end-1,:);

            b.R1_zz.w   = ones(b.K+2,1)*(b.eta_z.w.*diff(b.R1_z.u(1,:))./b.deta.zm);
            b.R1_zz.u   = ones(b.K+1,1)*[b.R1_zz.w(1) b.eta_z.u(2:end-1).*diff(b.R1_z.w(1,:))./b.deta.c b.R1_zz.w(end)];
            b.R1_zz.T   = [b.R1_zz.u; b.R1_zz.u(1,:)];
            b.R1_zz.v   = b.R1_zz.w(1:end-1,:);
            %}
            
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
        if strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'s')==1
            F_rlb_C     = @(R1, Rend, R1_z, Rend_z, R_rlbend, XI_r_rlbend, XI_z) XI_r_rlbend.*(Rend_z.*R_rlbend - R1_z.*(R_rlbend-1))-2*XI_z./(Rend-R1);
            F_rlb_E     = @(R, R1, Rend_z_rlbend, XI_r_rlbend) XI_r_rlbend.*Rend_z_rlbend.*(R-R1);
            F_inv_rlb   = @(R1, Rend) -(b.geom.r(end)-b.geom.r(1))./((Rend-R1).^2);

            if isfield(oldb,'Z')==0 || isequal(b.Z,oldb.Z)
                eval(['Rsend=x(Rs.' b.bc.z{length(b.bc.z)/2+n} '(1):Rs.' b.bc.z{length(b.bc.z)/2+n} '(2));'])
                
            else
                eval(['Rsend=interp1(oldb.ETA.u(1,(oldb.geom.z(n)<oldb.Z.u(1,:) & oldb.geom.z(n+1)>oldb.Z.u(1,:))).'',old.x(oldb.Rs.' b.bc.z{length(b.bc.z)/2+n} '(1):oldb.Rs.' b.bc.z{length(b.bc.z)/2+n} '(2)),b.ETA.u(1,(b.geom.z(n)<b.Z.u(1,:) & b.geom.z(n+1)>b.Z.u(1,:))).'',''spline'');'])
                
            end
            if eval(['b.' eval(['b.bc.z{' num2str(length(b.bc.z)/2+n) '}']) '.o==1'])
                b.Rsend=[b.Rsend; Rsend];
                
            end
            
            b.Rend.u(:,b.geom.z(n)<b.z.u & b.geom.z(n+1)>b.z.u)=ones(b.K+1,1)*Rsend.';
            b.Rend.w(:,b.geom.z(n)<b.z.w & b.geom.z(n+1)>b.z.w)=ones(b.K+2,1)*(b.W.w(1,(b.geom.z(n)<b.z.w & b.geom.z(n+1)>b.z.w)).*Rsend(1:end-1).'+b.E.w(1,find(b.geom.z(n)<b.z.w & b.geom.z(n+1)>b.z.w)-1).*Rsend(2:end).');
            b.Rend.T=[b.Rend.u; b.Rend.u(1,:)];
            b.Rend.v=b.Rend.w(1:end-1,:);
            %%{
            b.Rend_z.w   = ones(b.K+2,1)*(diff(b.Rend.u(1,:))./b.dz.zm);
            b.Rend_z.u   = ones(b.K+1,1)*[b.Rend_z.w(1) diff(b.Rend.w(1,:))./b.dz.c b.Rend_z.w(end)];
            b.Rend_z.T   = [b.Rend_z.u; b.Rend_z.u(1,:)];
            b.Rend_z.v   = b.Rend_z.w(1:end-1,:);

            b.Rend_zz.w   = ones(b.K+2,1)*(diff(b.Rend_z.u(1,:))./b.dz.zm);
            b.Rend_zz.u   = ones(b.K+1,1)*[b.Rend_zz.w(1) diff(b.Rend_z.w(1,:))./b.dz.c b.Rend_zz.w(end)];
            b.Rend_zz.T   = [b.Rend_zz.u; b.Rend_zz.u(1,:)];
            b.Rend_zz.v   = b.Rend_zz.w(1:end-1,:);
            %}
            
            %{
            b.Rend_z.w   = ones(b.K+2,1)*(b.eta_z.w.*diff(b.Rend.u(1,:))./b.deta.zm);
            b.Rend_z.u   = ones(b.K+1,1)*[b.Rend_z.w(1) b.eta_z.u(2:end-1).*diff(b.Rend.w(1,:))./b.deta.c b.Rend_z.w(end)];
            b.Rend_z.T   = [b.Rend_z.u; b.Rend_z.u(1,:)];
            b.Rend_z.v   = b.Rend_z.w(1:end-1,:);

            b.Rend_zz.w   = ones(b.K+2,1)*(b.eta_z.w.*diff(b.Rend_z.u(1,:))./b.deta.zm);
            b.Rend_zz.u   = ones(b.K+1,1)*[b.Rend_zz.w(1) b.eta_z.u(2:end-1).*diff(b.Rend_z.w(1,:))./b.deta.c b.Rend_zz.w(end)];
            b.Rend_zz.T   = [b.Rend_zz.u; b.Rend_zz.u(1,:)];
            b.Rend_zz.v   = b.Rend_zz.w(1:end-1,:);
            %}
            
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
        
    end
    
    oldb.Rs1=b.Rs1;
    oldb.Rsend=b.Rsend;
    oldb.Rs=Rs;
        
    R           = @(R1, Rend, R_cyl) R1+(R_cyl-b.geom.r(1)).*(Rend-R1)/(b.geom.r(end)-b.geom.r(1));
    R_z         = @(R1_z, Rend_z, R_cyl) R1_z+(R_cyl-b.geom.r(1)).*(Rend_z-R1_z)/(b.geom.r(end)-b.geom.r(1));
    %R_zz        = @(R1_zz, Rend_zz, R_cyl) R1_zz+(R_cyl-b.geom.r(1)).*(Rend_zz-R1_zz)/(b.geom.r(end)-b.geom.r(1));
    F           = @(R1, Rend) (b.geom.r(end)-b.geom.r(1))./(Rend-R1);
    R_cyl_z     = @(R1, Rend, R1_z, Rend_z, R, XI_r) (XI_r.*((-R1_z).*(Rend-R1)-(Rend_z-R1_z).*(R-R1))./(Rend-R1));
    
    b.R.u=R(b.R1.u, b.Rend.u, b.R_cyl.u);
    b.R.w=R(b.R1.w, b.Rend.w, b.R_cyl.w);
    b.R.T=R(b.R1.T, b.Rend.T, b.R_cyl.T);
    b.R.p=b.R.T(2:end-1,2:end-1);
    b.R.v=R(b.R1.v, b.Rend.v, b.R_cyl.v);
    
    %%{
    b.R_z.u=R_z(b.R1_z.u, b.Rend_z.u, b.R_cyl.u);
    b.R_z.w=R_z(b.R1_z.w, b.Rend_z.w, b.R_cyl.w);
    b.R_z.T=R_z(b.R1_z.T, b.Rend_z.T, b.R_cyl.T);
    b.R_z.p=b.R_z.T(2:end-1,2:end-1);
    b.R_z.v=R_z(b.R1_z.v, b.Rend_z.v, b.R_cyl.v);
    %}
    %{
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

    % del(XI)/del(z)
        b.XI_z.u   = R_cyl_z(b.R1.u, b.Rend.u, b.R1_z.u, b.Rend_z.u, b.R.u, b.XI_r.u);
        b.XI_z.w   = R_cyl_z(b.R1.w, b.Rend.w, b.R1_z.w, b.Rend_z.w, b.R.w, b.XI_r.w);
        b.XI_z.T   = R_cyl_z(b.R1.T, b.Rend.T, b.R1_z.T, b.Rend_z.T, b.R.T, b.XI_r.T);
        b.XI_z.p   = b.XI_z.T(2:end-1,2:end-1);
        b.XI_z.v   = R_cyl_z(b.R1.v, b.Rend.v, b.R1_z.v, b.Rend_z.v, b.R.v, b.XI_r.v);
        
    % Jacobi determinant
        b.JA.u=1./(b.XI_r.u.*b.ETA_z.u);
        b.JA.w=1./(b.XI_r.w.*b.ETA_z.w);
        b.JA.T=1./(b.XI_r.T.*b.ETA_z.T);
        b.JA.p=b.JA.T(2:end-1,2:end-1);
        b.JA.v=1./(b.XI_r.v.*b.ETA_z.v);