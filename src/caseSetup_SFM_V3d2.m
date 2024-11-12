% flowoptions --> don't change because for inviscid flows no solution could be found
    flowopt.inviscid = 0;

% stability properties
    % critical curve settings
        % sense of rotation for finding the critical point
            rotdir = 1;
            
        % imaginary eigenvalue problem
            flowopt.iv = 1;

        % set independent parameters from the critical curve
            independent = 'xxx';                   % only needed for critical curve
            eval([independent '_start = 0;'])      % only needed for critical curve
            eval([independent '_end = 0.1;'])      % only needed for critical curve
            eval([independent '_resolution = 1;']) % only needed for critical curve

        % set dependent parameters from the critical curve
            eval([dependent '_start = ' num2str(T_init) ';'])
            eval([dependent '_resolution = 1;'])   % only needed for critical curve

    
% geometry
    blocks={'b1'};

    b1.geom.r=[r_c r_i];
    b1.geom.z=[0 l_lb];

% fluid pressure
    if isfile([tempdir 'my_fluids.mat'])
        my_fluids = load([tempdir 'my_fluids.mat']);
    end
    b1.p_fluid=00000;
    b2.p_fluid=00000;

    % location
        b1.rp=(b1.geom.r(1) + b1.geom.r(end))/2;
        b1.zp=(b1.geom.z(1) + b1.geom.z(end))/2;

% fluid properties
    T_0 = mean([T_d1l T_d2l]);
    if flowopt.boussinesq == 1 || NS == 1 % this is equivalent
        flowopt.boussinesq = 1; b1.T_0 = T_0;
        [b1]=fluid_properties_GUI_V3d2(selectedLiquid,'liquid', [], p_fluid, b1);

        % find out beta
        rho_l   = str2func(b1.rho);
        drho_l  = str2func(b1.drho);
        b1.beta = -drho_l(T_0)/rho_l(T_0);
        [b1]=fluid_properties_GUI_V3d2(selectedLiquid,'liquid', T_0, p_fluid, b1);

    elseif NS == 2
        [b1]=fluid_properties_GUI_V3d2(selectedLiquid,'liquid', [], p_fluid, b1);

        % enforce linear dependence
        rho_l    = str2func(b1.rho);    drho_l    = str2func(b1.drho);
        mu_l     = str2func(b1.mu);     dmu_l     = str2func(b1.dmu);
        lambda_l = str2func(b1.lambda); dlambda_l = str2func(b1.dlambda);
        cp_l     = str2func(b1.cp);     dcp_l     = str2func(b1.dcp);
        sigma_l  = str2func(b1.sigma);  dsigma_l = str2func(b1.dsigma);
        b1.rho    = ['@(theta) ' num2str(rho_l(T_0),16)    '+' num2str(drho_l(T_0),16)    '*(theta-' num2str(T_0) ')'];
        b1.mu     = ['@(theta) ' num2str(mu_l(T_0),16)     '+' num2str(dmu_l(T_0),16)     '*(theta-' num2str(T_0) ')'];
        b1.lambda = ['@(theta) ' num2str(lambda_l(T_0),16) '+' num2str(dlambda_l(T_0),16) '*(theta-' num2str(T_0) ')'];
        b1.cp     = ['@(theta) ' num2str(cp_l(T_0),16)     '+' num2str(dcp_l(T_0),16)     '*(theta-' num2str(T_0) ')'];
        b1.sigma  = ['@(theta) ' num2str(sigma_l(T_0),16)  '+' num2str(dsigma_l(T_0),16)  '*(theta-' num2str(T_0) ')'];
        b1.drho    = ['@(theta) ' num2str(drho_l(T_0),16)];
        b1.dmu     = ['@(theta) ' num2str(dmu_l(T_0),16)];
        b1.dlambda = ['@(theta) ' num2str(dlambda_l(T_0),16)];
        b1.dcp     = ['@(theta) ' num2str(dcp_l(T_0),16)];
        b1.dsigma  = ['@(theta) ' num2str(dsigma_l(T_0),16)];

    else
        [b1]=fluid_properties_GUI_V3d2(selectedLiquid,'liquid', T_fluid, p_fluid, b1);
    end
    [b2]=fluid_properties_GUI_V3d2(selectedGas,'gas', T_fluid, p_fluid, b2);


% initial conditions
    b1.u_init.u='@(Z,R) 0*Z + 0*R';
    b1.w_init.w='@(Z,R) 0*Z + 0*R';
    b1.p_init.p=['@(Z,R) ones(size(Z))*' num2str(b1.p_fluid)];
    b1.T_init.T=['@(Z,R) ones(size(Z))*' num2str((T_d1l+T_d2l)/2)];

% Volume ratio V_r=Volume_liquid/Volume_cylinder
    b1.V_r=V_r;

%static surface shape
    if strcmp(b1.bc.z{2}(1),'s')==1 && strcmp(b1.bc.z{2}(2),'i')==1
        if exist('oldb1','var')==0 || isfield(oldb1, 'dp_lb')==0
            oldb1.r_lb_l=[];
            oldb1.dp_lb=[];
        end
            rho_l=str2func(b1.rho);
            rho_g=str2func(b2.rho);
            sigma=str2func(b1.sigma);
    
            [b1.r_lb_l, b1.dp_lb, convergence]=surface_shape(oldb1.r_lb_l, oldb1.dp_lb, b1.V_r, -flowopt.g, rho_l((T_d1l+T_d2l)/2), rho_g((T_d1l+T_d2l)/2), sigma((T_d1l+T_d2l)/2), r_c, r_i, l_lb, 0, 0.0005*l_lb);
            oldb1.r_lb_l=b1.r_lb_l;
            oldb1.dp_lb=b1.dp_lb;
            shape='b.r_lb_l';
    else
        shape=['@(z)' num2str(b1.geom.r(2))];
        convergence = 1;
    end

% functions on boundaries
    b1.bc.rhs.u.z = {'@(z)0'; '@(z)0'};
    b1.bc.rhs.w.z = {'@(z)0'; '@(z)0'};
    b1.geom.sh.z  = {['@(z)' num2str(r_c)]; shape};

    b1.bc.rhs.u.r={'@(r)0'; '@(r)0'};
    b1.bc.rhs.w.r={'@(r)0'; '@(r)0'};

% grid parameters
    geom.r=[r_c r_i];
    geom.z=[0 l_lb];
    
    [mesh.r.u, mesh.r.w, mesh.xi_r.u, mesh.xi_r.w] = segmentation(geom.r, mesh.r); 
    [mesh.z.w, mesh.z.u, mesh.eta_z.w, mesh.eta_z.u] = segmentation(geom.z, mesh.z);