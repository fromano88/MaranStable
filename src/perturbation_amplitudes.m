vector_index = 0;
for n = 1:length(blocks)
    eval([blocks{n} ' = amplitudes(eigenvector,vector_index,flowopt,' blocks{n} ');'])
    vector_index = vector_index+(b1.K+1)*(b1.J+2)+(b1.K+2)*(b1.J+2)+(b1.K+2)*(b1.J+1)+(b1.K)*(b1.J)+(b1.K+2)*(b1.J+2);
end

for n=1:length(blocks)
    b = eval(blocks{n});

    if flowopt.energy == 1

        if flowopt.ax == 0
            time_T = 0;
            phi_T  = 0;                                % homogeneous direction
            z_T    = b1.geom.z(1)+0.5*diff(b1.geom.z); % center of "LB"
            block_T_max = plus_cc(b.T_hat.T,m,phi_T);
            block_T_max = unique(b.T_hat.T(abs(block_T_max)==max(max(abs(block_T_max)))));

            time_conv_energy = 0;
            phi_conv_energy  = 0;
            z_conv_energy    = b1.geom.z(1)+0.5*diff(b1.geom.z);
            block_conv_energy_max = plus_cc(b.conv_energy_hat,m,phi_conv_energy);
            block_conv_energy_max = unique(b.conv_energy_hat(abs(block_conv_energy_max)==max(max(abs(block_conv_energy_max)))));

        elseif m ~= 0
            time_T = 0;
            phi_T  = 1/m*atan(-imag(b.T_hat.T)./real(b.T_hat.T));
            block_T_max = plus_cc(b.T_hat.T,m,phi_T);
    
            phi_T = phi_T(abs(block_T_max)==max(max(abs(block_T_max))));
            z_T   = b.Z.T(abs(block_T_max)==max(max(abs(block_T_max))));
            block_T_max = b.T_hat.T(abs(block_T_max)==max(max(abs(block_T_max))));
            if imag(block_T_max(1))*sin(m*phi_T)-real(block_T_max(1))*cos(m*phi_T)>0
                phi_T = phi_T+pi/m;
            end
    
            time_conv_energy = 0;
            phi_conv_energy  = -1/m*atan(imag(b.conv_energy_hat)./real(b.conv_energy_hat));
            block_conv_energy_max = plus_cc(b.conv_energy_hat,m,phi_conv_energy);
            phi_conv_energy = phi_conv_energy(abs(block_conv_energy_max)==max(max(abs(block_conv_energy_max))));
            z_conv_energy   = b.Z.T(abs(block_conv_energy_max)==max(max(abs(block_conv_energy_max))));
            block_conv_energy_max=b.conv_energy_hat(abs(block_conv_energy_max)==max(max(abs(block_conv_energy_max))));
            if imag(block_conv_energy_max)*sin(m*phi_conv_energy)-real(block_conv_energy_max)*cos(m*phi_conv_energy)>0
                phi_conv_energy = phi_conv_energy+pi/m;
            end
    
        else
            phi_T  = 0;
            if imag(gamma(1,1)) == 0
                time_T = zeros(size(b.T.T));
            else
                time_T = -atan((real(b.T_hat.T)*real(gamma(1,1)) - imag(b.T_hat.T)*imag(gamma(1,1)))./(real(b.T_hat.T)*imag(gamma(1,1)) + imag(b.T_hat.T)*real(gamma(1,1))))/imag(gamma(1,1));
            end
            block_T_max = plus_cc(b.T_hat.T,0,0,gamma,time_T);
            time_T = time_T(abs(block_T_max)==max(max(abs(block_T_max))));
            z_T    = b.Z.T(abs(block_T_max)==max(max(abs(block_T_max))));
            block_T_max=b.T_hat.T(abs(block_T_max)==max(max(abs(block_T_max))));
            if -imag(gamma(1,1))*sign(real(block_T_max)*imag(gamma(1,1)) + imag(block_T_max)*real(gamma(1,1)))>0
                time_T = time_T+pi/imag(gamma(1,1));
            end
    
            phi_conv_energy  = 0;
            if imag(gamma(1,1)) == 0
                time_conv_energy = zeros(size(b.conv_energy_hat));
            else
                time_conv_energy = -atan((real(b.conv_energy_hat)*real(gamma(1,1)) - imag(b.conv_energy_hat)*imag(gamma(1,1)))./(real(b.conv_energy_hat)*imag(gamma(1,1)) + imag(b.conv_energy_hat)*real(gamma(1,1))))/imag(gamma(1,1));
            end
            block_conv_energy_max = plus_cc(b.conv_energy_hat,0,0,gamma,time_conv_energy);
            time_conv_energy = time_conv_energy(abs(block_conv_energy_max)==max(max(abs(block_conv_energy_max))));
            z_conv_energy    = b.Z.T(abs(block_conv_energy_max)==max(max(abs(block_conv_energy_max))));
            block_conv_energy_max = b.conv_energy_hat(abs(block_conv_energy_max)==max(max(abs(block_conv_energy_max))));
            if -imag(gamma(1,1))*sign(real(block_conv_energy_max)*imag(gamma(1,1)) + imag(block_conv_energy_max)*real(gamma(1,1)))>0
                time_conv_energy = time_conv_energy+pi/imag(gamma(1,1));
            end
    
        end

        if max(isinf(time_T)) == 1
            time_T = 0;
        end
        if max(isinf(phi_T)) == 1
            phi_T = 0;
        end
        if max(isinf(z_T)) == 1
            z_T = b1.geom.z(1)+0.5*diff(b1.geom.z);
        end
        if max(isinf(time_conv_energy)) == 1
            time_conv_energy = 0;
        end
        if max(isinf(phi_conv_energy)) == 1
            phi_conv_energy = 0;
        end
        if max(isinf(block_T_max)) == 1
            block_T_max = 0;
        end
    
        T_ex = plus_cc(block_T_max,m,phi_T,gamma,time_T);
        if isfield(b,'sigma') || T_ex>T_max
            T_max   = abs(T_ex);
            z_T_max = z_T;
            if m~=0
                phi_T_max  = phi_T;
                time_T_max = 0;
            else
                time_T_max = time_T;
                phi_T_max  = 0;
            end
        end
    
        conv_energy_ex = plus_cc(block_conv_energy_max,m,phi_conv_energy,gamma,time_conv_energy);
        b.conv_energy_max = conv_energy_ex;
        if isfield(b,'sigma') || conv_energy_ex>conv_energy_max
            conv_energy_max   = conv_energy_ex;
            z_conv_energy_max = z_conv_energy;
            if m~=0
                phi_conv_energy_max  = phi_conv_energy;
                time_conv_energy_max = 0;
    
            else
                time_conv_energy_max = time_conv_energy;
                phi_conv_energy_max  = 0;
    
            end
    
        end
    
    else
        z_T_max    = b1.geom.z(1)+0.5*diff(b1.geom.z);
        phi_T_max  = 0;
        time_T_max = 0;
    end

end

% find position of maximal thermal energy production
% rho_0*T'*(u'*d(c_p0 T_0)/dr + w'*d(c_p0 T_0)/dz)
if flowopt.energy == 1 %&& flowopt.ax == 1
    % thermal production in the whole LB
        Nphi = 100;
        [thermal_production, phi] = thermal_energy_production(b1,gamma,time_T_max,m,Nphi);

    % thermal production in the core region
        sx = size(b1.R.v,1); sy = size(b1.R.v,2); sz = length(phi);
        l_d1 = b1.geom.z(1);
        z_dimless = (l_d1+l_lb/2-b1.Z.v(1,:))/l_lb;
        idx_z1 = find(z_dimless<0.45,1);
        idx_z2 = find(z_dimless>-0.45,1,'last');
        idx_r  = find(b1.R.v(:,1)<0.9*r_i,1,'last');

        thermal_prod = reshape(thermal_production,sx,sy,sz);
        thermal_prod_core = thermal_prod(1:idx_r,idx_z1:idx_z2,:);
        thermal_prod_core_2D = reshape(thermal_prod_core,size(thermal_prod_core,1)*size(thermal_prod_core,2),sz);
        R_core = b1.R.v(1:idx_r,idx_z1:idx_z2,:);
        Z_core = b1.Z.v(1:idx_r,idx_z1:idx_z2,:);

    [production_max,idx_max] = max(thermal_prod_core_2D(:));
    idx_phi = ceil(idx_max/numel(R_core));
    if mod(idx_max,numel(R_core)) == 0
        idx_rz = numel(R_core);
    else
        idx_rz = mod(idx_max,numel(R_core));
    end
    idx_z = ceil(idx_rz/size(R_core,1));
    if mod(idx_rz,size(R_core,1)) == 0
        idx_r = size(R_core,1);
    else
        idx_r = mod(idx_rz,size(R_core,1));
    end

    r_T_max   = R_core(idx_r,idx_z);
    z_T_max   = Z_core(idx_r,idx_z);
    phi_T_max = phi(idx_phi);
    if sum(sum(plus_cc(b1.T_hat.v,m,phi_T_max,gamma,time_T_max))) < 0
        phi_T_max = phi_T_max+pi/m;
    end
end