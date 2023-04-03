% SFM
    h.blocks = {'b1'};

% flowOptions
    % --- do not change this default values ---
    h.flowopt.ax                  = 1;
    h.flowopt.energy              = 1;
    h.flowopt.thermcapcon         = 1;
    h.flowopt.creeping            = 0;
    h.flowopt.boussinesq          = 0;
    % -----------------------------------------
    h.flowopt.tolerance.residuals = 1e-6;
    h.flowopt.tolerance.newton    = 1e-6;
    h.flowopt.tolerance.growth    = 1e-5;
    h.flowopt.g                   = 0*9.81;
    h.flowopt.eigs.n              = 12;
    h.flowopt.eigs.n_cayley       = 5;
    h.flowopt.eigs.tol            = 1e-14;
    h.flowopt.eigs.maxit          = 300;
    h.flowopt.eigs.krylov         = 100;
    % -----------------------------------------
    h.NS = 3; % 1...Boussinesq, 2...linear on T, 3...fully temperature dep.

% geometry in m
    h.r_c  = 0;
    h.r_i  = 0.004;
    h.l_lb = 0.004;

% temperature
    h.T_d1l   = 40;
    h.T_d2l   = 10;
    h.delta_T = h.T_d1l-h.T_d2l;
    h.T_0     = (h.T_d1l+h.T_d2l)/2;
    
% stability
    h.m_start = 0;
    h.m_delta = 1;
    h.m_end   = 1;
    h.T_init  = h.T_d1l-h.T_d2l;
    h.dependent = 'delta_T';
    h.equal_T_d1 = {'T_d1l'};
    h.equal_T_d2 = {'T_d2l'};

% volume ratio
    h.V_r = 1;

% fluids
    h.selectedLiquid = 'KF-96L-5cs';
    h.selectedGas    = 'Air';

    % possible liquids
    h.liquid_list = {'KF-96-5cs','KF-96L-2cs','KF-96L-20cs','H2O','Br2','HF',...
    'HCN','HBr','PCl3','SiCl4','CS2','C2H6S','C4H4S','CH2Cl2',...
    'CHCl3','CCl4','CH2Br2','CHBr3','CCl3F','C2H5Br','C2H4Cl2',...
    'C2H4Br2','C2H3Cl3','C3H7Cl','C4H9Cl','C5H11Cl','C2HCl3','C2Cl4',...
    'C6H5F','C6H5Cl','C6H5Br','C6H5I','C7H7Cl','C6H14','C7H16','C8H18',...
    'C9H20','C10H22','C11H24','C12H26','C6H12','C7H14','C8H16','C9H18',...
    'C10H20','C11H22','C12H24','C6H10','C6H6','C7H8','C8H8','C13H12',...
    'CH4O','C2H6O','C3H8O','C4H10O','C5H12O','C6H14O','C7H16O','C8H18O',...
    'C2H6O2','C3H8O2','C3H8O3','C6H12O','C7H8O','C6H6O','C2H4O2',...
    'C3H6O2','C4H8O2','C5H10O2','C6H12O2','C4H6O3','C6H10O3','C2H2Cl2O2',...
    'C3H6O','C4H8O','C5H10O','C7H14O','C8H8O','C6H12O3','C5H4O2','C7H6O',...
    'C7H6O2','C8H8O2','C9H10O2','C8H8O3','C3H9N','C5H11N','C7H9N',...
    'C8H11N','C10H15N','C6H8N2','C3H5N','C4H7N','C7H5N','C6H5NO2',...
    'CH2O2','C2H3N','CH3NO','CH3NO2'};

% grid parameters
    h.mesh.r.delta.start = 0.1;
    h.mesh.r.delta.fit   = 0.75;
    h.mesh.r.delta.end   = 0.005;
    h.mesh.r.f           = [1.15 1.15];
    h.mesh.r.sf          = {'tanh'};
    
    h.mesh.z.delta.start = 0.005;
    h.mesh.z.delta.fit   = 0.75;
    h.mesh.z.delta.end   = 0.005;
    h.mesh.z.f           = [1.15 1.15];
    h.mesh.z.sf          = {'tanh'};

% boundary conditions
    h.Bi = 0;
    h.b1.bc.rhs.T.z = {'@(z)0'; ['@(z)' num2str(h.Bi)]};
    h.b1.bc.rhs.T.r = {['@(r)' num2str(h.T_d1l)]; ['@(r)' num2str(h.T_d2l)]};

    h.b1.bc.z(1,1) = 'a';
    %------------------------%
    h.b1.bc.z(2,1) = 's';
    h.b1.bc.z(2,2) = 'i';
    h.b1.bc.z(2,3) = 'a';
    h.b1.bc.z(2,4) = '0';
    %------------------------%
    h.b1.bc.r(1,1) = 'w';
    h.b1.bc.r(1,2) = 'n';
    h.b1.bc.r(1,3) = 'c';
    %------------------------%
    h.b1.bc.r(2,1) = 'w';
    h.b1.bc.r(2,2) = 'n';
    h.b1.bc.r(2,3) = 'c';
    
% initialization
    h.initialize = 'standard';

% number of iteration steps
    h.steps = 100;
    
% optical ray options
    h.r_0 = 0.8;
    h.z_p = -0.4;
    h.N_coeff(1) = 1.4094; % N(T) = N_coeff(1)-N_coeff(2)*(T-N_coeff(3))
    h.N_coeff(2) = 3.9e-4;
    h.N_coeff(3) = 20;

% temperatures at the corners
    h.temperature.rc_d1 = h.T_d1l;
    h.temperature.rc_d2 = h.T_d2l;
    h.temperature.ri_d1 = h.T_d1l;
    h.temperature.ri_d2 = h.T_d2l;
    
% list of boundary conditions for temperature and velocity
    h.boundaryNames = {'radialLB.rc', 'axialLB.d1', 'axialLB.d2'};
    
% text position (left, right)
    h.tp_l = -0.09; h.tp_r = 1.1;
    if ispc
        h.tp_y = 0.5;
    else
        h.tp_y = 0.47;
    end
    
% Font sizes
    h.latex_l = {'interpreter','Latex','FontSize',12,'HorizontalAlignment','right'};
    h.latex_r = {'interpreter','Latex','FontSize',12,'HorizontalAlignment','left'};
    h.Latex = {'interpreter','Latex','FontSize',22};
    h.LATEX = {'interpreter','Latex','FontSize',26};
    
% Tab Colors
    h.unselectedTabColor = [0 0.4470 0.7410];
    h.selectedTabColor   = [0.8500 0.3250 0.0980];
    h.blue   = [0 0.4470 0.7410];
    h.red    = [0.8500 0.3250 0.0980];
    h.ocher  = [0.9290 0.6940 0.1250];
    h.violet = [0.4940 0.1840 0.5560];
    h.green  = [0.4660 0.6740 0.1880];
    h.colors = [h.blue; h.red; h.ocher; h.violet; h.green; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];
    
% Windows style
    h.opts = struct('WindowStyle','non-modal','Interpreter','tex');
    
% every message box which is listed here will be shown only once
    % edit texts - save writings
        % inputValues
            sw = ["g","V_r","l_lb","r_c","r_i","spacing","stretchingC","temperature","residuals","steps","r_0","z_p","N","n_eig","m_start","m_delta","m_end","m","n_cayley","tolerance"];
            for i = 1:length(sw)
                h.show.(sw(i)).error   = 1;
                h.show.(sw(i)).warning = 1;
            end
        % additional edit texts
            h.show.stretching = 1;
            
    % button groups
        sw = ["NS1","NS2","stretchingF_lin","stretchingF","surface"];
        for i = 1:length(sw)
            h.show.(sw(i)) = 1;
        end

% visibilty of subpanels
    h.mh.Visible = 'axial';
    h.bc.Visible = 'axial';