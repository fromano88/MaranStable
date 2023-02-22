% enters if statement if file is not called from GUI
if isempty(get(gcf,'Name'))
    close gcf
    
    domain = [1 0]; dim = 0; unity = 'm'; Nphi = 100;
    export.T = 1; export.v = 0; export.p = 0; export.th_E = 1; % exporting variables
    export.bS = 0; export.pF = 1; % bS ... basic state, pF ... perturbation flow
    file_path = [folder 'newData.vtk'];
end

if domain(2) == 1
    b = b2;
else
    b = b1;
end

% unity
if dim == 0
    f = 1/l_lb;
else
    switch unity
        case 'm'
            f = 1;
        case 'dm'
            f = 10;
        case 'cm'
            f = 100;
        case 'mm'
            f = 1000;
    end
end

if isfield(b,'w_hat')
    if imag(gamma(1,1)) > 0
        gamma(1,1)  = conj(gamma(1,1));
        eigenvector = conj(eigenvector);
    end
else
    export.pF = 0;
end

delta_T = abs(T_d1l-T_d2l);
T_0     = (T_d1l+T_d2l)/2;
if flowopt.energy == 0
    export.th_E = 0;
    
    b.T.v = T_0*ones(size(b.R.v));
    b.T_hat.v = zeros(size(b.R.v));
end

% for non-dimensionalization use always the properties of the liquid phase
mu     = str2func(b1.mu);
rho    = str2func(b1.rho);
dsigma = str2func(b1.dsigma);
nu  = mu(T_0)/rho(T_0);

R_2D  = b.R.v*f;
if length(blocks) == 1
    l_d1 = 0;
end
Z_2D  = (l_d1+l_lb/2-b.Z.v)*f;
if dim == 1
    Ur_2D = (1./b.XI_r.v.*b.u.v-b.JA.v.*b.XI_z.v.*b.w.v)*f;
    Uz_2D = -1./b.ETA_z.v.*b.w.v*f;
    T_2D  = b.T.v;
    p_2D  = b.p.v;
else
    Ur_2D = (1./b.XI_r.v.*b.u.v-b.JA.v.*b.XI_z.v.*b.w.v)*l_lb/nu;
    Uz_2D = -1./b.ETA_z.v.*b.w.v*l_lb/nu;
    T_2D  = (b.T.v-T_0)/delta_T;
    p_2D  = b.p.v*l_lb/(abs(dsigma(T_0))*delta_T);
end
%% Reconstruct Perturbation
N_t = 1;
phi = 0:2*pi/Nphi:2*pi;
R_3D    = repmat(R_2D(:),1,length(phi));
Phi_3D  = repmat(phi,length(R_2D(:)),1);
Z_3D    = repmat(Z_2D(:),1,length(phi));
T_3Dbs  = T_2D(:).*exp(0*phi);
p_3Dbs  = p_2D(:).*exp(0*phi);
Ur_3Dbs = Ur_2D(:).*exp(0*phi);
Uz_3Dbs = Uz_2D(:).*exp(0*phi);
if export.pF == 1
    A_pert = 3000;
    if dim == 1
        omega = imag(gamma);
    else
        omega = imag(gamma)*l_lb^2/nu;
    end
    if N_t == 1
        t = 0;
    else
        t = (0:1/N_t:1)*2*pi/omega;
    end

    % since the perturbations are defined only up to a constant factor, we can leave it in dimensional form
    T_pert    = b.T_hat.v(:).*exp(-gamma*t+1i*m*phi)+conj(b.T_hat.v(:).*exp(-gamma*t+1i*m*phi));
    p_pert    = b.p_hat.v(:).*exp(-gamma*t+1i*m*phi)+conj(b.p_hat.v(:).*exp(-gamma*t+1i*m*phi));
    Ur_pert   = 1./b.XI_r.v(:).*(b.u_hat.v(:).*exp(-gamma*t+1i*m*phi)+conj(b.u_hat.v(:).*exp(-gamma*t+1i*m*phi))) - b.JA.v(:).*b.XI_z.v(:).*(b.w_hat.v(:).*exp(-gamma*t+1i*m*phi)+conj(b.w_hat.v(:).*exp(-gamma*t+1i*m*phi)));
    Uphi_pert = b.v_hat.v(:)*exp(-gamma*t+1i*m*phi)+conj(b.v_hat.v(:)*exp(-gamma*t+1i*m*phi));
    Uz_pert   = -1./b.ETA_z.v(:).*(b.w_hat.v(:).*exp(-gamma*t+1i*m*phi)+conj(b.w_hat.v(:).*exp(-gamma*t+1i*m*phi)));
    
    T_tot     = T_3Dbs  + A_pert*T_pert;
    p_tot     = p_3Dbs  + A_pert*p_pert;
    Ur_tot    = Ur_3Dbs + A_pert*Ur_pert;
    Uphi_tot  =           A_pert*Uphi_pert-omega*R_3D; % corotating with the wave
    Uz_tot    = Uz_3Dbs + A_pert*Uz_pert;
end
%% Reconstruct Cartesian Grid Components
% repmat on R and Z and match them with phi to create a 3D grid [R_3D,Z_3D]
Ux_3Dbs = Ur_3Dbs.*cos(Phi_3D);
Uy_3Dbs = Ur_3Dbs.*sin(Phi_3D);
if export.pF == 1
    Ux_pert = Ur_pert.*cos(Phi_3D) - Uphi_pert.*sin(Phi_3D);
    Uy_pert = Ur_pert.*sin(Phi_3D) + Uphi_pert.*cos(Phi_3D);
    Ux_tot  = Ur_tot.*cos(Phi_3D)  - Uphi_tot.*sin(Phi_3D);
    Uy_tot  = Ur_tot.*sin(Phi_3D)  + Uphi_tot.*cos(Phi_3D);
end
[X_3D,Y_3D] = pol2cart(Phi_3D,R_3D);
%% Reshape matrices
x = size(R_2D,1); y = size(R_2D,2); z = length(phi);
X_3D    = reshape(X_3D,x,y,z);
Y_3D    = reshape(Y_3D,x,y,z);
Z_3D    = reshape(Z_3D,x,y,z);
Ur_3Dbs = reshape(Ur_3Dbs,x,y,z);
Ux_3Dbs = reshape(Ux_3Dbs,x,y,z);
Uy_3Dbs = reshape(Uy_3Dbs,x,y,z);
Uz_3Dbs = reshape(Uz_3Dbs,x,y,z);
T_3Dbs  = reshape(T_3Dbs,x,y,z);
p_3Dbs  = reshape(p_3Dbs,x,y,z);
if export.pF == 1
    Ur_pert   = reshape(Ur_pert,x,y,z);
    Uphi_pert = reshape(Uphi_pert,x,y,z);
    Ur_tot    = reshape(Ur_tot,x,y,z);
    Uphi_tot  = reshape(Uphi_tot,x,y,z);
    Ux_pert   = reshape(Ux_pert,x,y,z);
    Uy_pert   = reshape(Uy_pert,x,y,z);
    Uz_pert   = reshape(Uz_pert,x,y,z);
    Ux_tot    = reshape(Ux_tot,x,y,z);
    Uy_tot    = reshape(Uy_tot,x,y,z);
    Uz_tot    = reshape(Uz_tot,x,y,z);
    T_pert    = reshape(T_pert,x,y,z);
    T_tot     = reshape(T_tot,x,y,z);
    p_pert    = reshape(p_pert,x,y,z);
    p_tot     = reshape(p_tot,x,y,z);
    if export.th_E == 1
        thermal_production = thermal_energy_production(b,gamma,time_T_max,m,Nphi);
        thermal_production = reshape(thermal_production,x,y,z);
    end
end
%% Convert Data from Matlab to ParaView
nr_of_elements=numel(X_3D);
fid = fopen(file_path, 'w'); 
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'BINARY\n\n');
fprintf(fid, 'DATASET STRUCTURED_GRID\n');
fprintf(fid, ['DIMENSIONS ' num2str(size(X_3D,1)) ' ' num2str(size(X_3D,2)) ' ' num2str(size(X_3D,3)) '\n']);
fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
fclose(fid);
fid = fopen(file_path, 'a'); 
fwrite(fid, [reshape(X_3D,1,nr_of_elements);  reshape(Y_3D,1,nr_of_elements); reshape(Z_3D,1,nr_of_elements)],'float','b');
fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
if export.v == 1
    if export.bS == 1
        fprintf(fid, 'VECTORS velocity_bS_cartesian float\n');
        fwrite(fid, [reshape(Ux_3Dbs,1,nr_of_elements);  reshape(Uy_3Dbs,1,nr_of_elements); reshape(Uz_3Dbs,1,nr_of_elements)],'float','b');
    end
    if export.pF == 1
        fprintf(fid, '\nVECTORS velocity_pF_cartesian float\n');
        fwrite(fid, [reshape(Ux_pert,1,nr_of_elements);  reshape(Uy_pert,1,nr_of_elements); reshape(Uz_pert,1,nr_of_elements)],'float','b');
    end
end
if export.T == 1
    if export.bS == 1
        fprintf(fid, '\nSCALARS T_bS float\n'); 
        fprintf(fid, 'LOOKUP_TABLE default\n'); 
        fwrite (fid, reshape(T_3Dbs,1,nr_of_elements),'float','b');
    end
    if export.pF == 1
        fprintf(fid, '\nSCALARS T_pF float\n'); 
        fprintf(fid, 'LOOKUP_TABLE default\n'); 
        fwrite (fid, reshape(T_pert,1,nr_of_elements),'float','b');
    end
end
if export.p == 1
    if export.bS == 1
        fprintf(fid, '\nSCALARS p_bS float\n'); 
        fprintf(fid, 'LOOKUP_TABLE default\n'); 
        fwrite (fid, reshape(p_3Dbs,1,nr_of_elements),'float','b');
    end
    if export.pF == 1
        fprintf(fid, '\nSCALARS p_pF float\n'); 
        fprintf(fid, 'LOOKUP_TABLE default\n'); 
        fwrite (fid, reshape(p_pert,1,nr_of_elements),'float','b');
    end
end
if export.pF == 1 && export.th_E == 1
    fprintf(fid, '\nSCALARS thermal_production float\n'); 
    fprintf(fid, 'LOOKUP_TABLE default\n'); 
    fwrite (fid, reshape(thermal_production,1,nr_of_elements),'float','b');
end
fclose(fid);