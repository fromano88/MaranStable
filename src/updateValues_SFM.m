assignin('base','tmp',h.flowopt.ax);
evalin('base','flowopt.ax=tmp;');
% ------------------------------------------------------------------- %
assignin('base','tmp',h.flowopt.energy);
evalin('base','flowopt.energy=tmp;');
% ------------------------------------------------------------------- %
assignin('base','tmp',h.flowopt.thermcapcon);
evalin('base','flowopt.thermcapcon=tmp;');
% ------------------------------------------------------------------- %
assignin('base','tmp',h.flowopt.creeping);
evalin('base','flowopt.creeping=tmp;');
% ------------------------------------------------------------------- %
assignin('base','tmp',h.flowopt.boussinesq);
evalin('base','flowopt.boussinesq=tmp;');
% ------------------------------------------------------------------- %
assignin('base','r_c',h.r_c);
assignin('base','d_c',h.r_c*2);
% ------------------------------------------------------------------- %
assignin('base','r_i',h.r_i);
assignin('base','d_i',h.r_i*2);
% ------------------------------------------------------------------- %
assignin('base','l_lb',h.l_lb);
% ------------------------------------------------------------------- %
assignin('base','p_fluid',101325);

% stability
    assignin('base','tmp',h.flowopt.eigs.n);
    evalin('base','flowopt.eigs.n=tmp;');
    % ------------------------------------------------------------------- %
    assignin('base','tmp',h.flowopt.eigs.n_cayley);
    evalin('base','flowopt.eigs.n_cayley=tmp;');
    % ------------------------------------------------------------------- %
    assignin('base','m_start',h.m_start);
    assignin('base','m_delta',h.m_delta);
    assignin('base','m_end',h.m_end);
    % ------------------------------------------------------------------- %
    assignin('base','dependent',h.dependent);
    assignin('base','T_0',h.T_0);
    assignin('base','T_init',h.T_init);
    % ------------------------------------------------------------------- %
    assignin('base','equal_T_d1',h.equal_T_d1);
    assignin('base','equal_T_d2',h.equal_T_d2);

% convergence criteria
    assignin('base','tmp',h.flowopt.tolerance.residuals);
    evalin('base','flowopt.tolerance.residuals=tmp;');
    % --------------------------------------------------------------- %
    assignin('base','tmp',h.flowopt.tolerance.newton);
    evalin('base','flowopt.tolerance.newton=tmp;');
    
% gravity
    assignin('base','tmp',h.flowopt.g);
    evalin('base','flowopt.g=tmp;');
    
% volume ratio
    assignin('base','V_r',h.V_r);
    
% (mean) temperatures
    assignin('base','T_d1l',h.T_d1l);
    assignin('base','T_d2l',h.T_d2l);
    assignin('base','T_initial',(h.T_d1l+h.T_d2l)/2);
    assignin('base','delta_T',h.delta_T);

% boundary conditions - general
    if h.r_c == 0
        assignin('base','tmp',{h.b1.bc.z(1,1); h.b1.bc.z(2,:)});
    else
        assignin('base','tmp',{h.b1.bc.z(1,1:3); h.b1.bc.z(2,:)});
    end
    evalin('base','b1.bc.z=tmp;');
    
    assignin('base','tmp',{h.b1.bc.r(1,1:3); h.b1.bc.r(2,1:3)});
    evalin('base','b1.bc.r=tmp;');

% fluid properies
    if h.NS == 1
        h.T_fluid = (h.T_d1l+h.T_d2l)/2;
    else
        h.T_fluid = [];
    end
    assignin('base','T_fluid',h.T_fluid);
    % --------------------------------------------------------------- %
    assignin('base','NS',h.NS);
    
% boundary conditions - functions on boundaries
    if strcmp(h.b1.bc.z(2,3),'a')
        h.Bi = 0;
    end
    assignin('base','Bi',h.Bi)
    if h.r_c == 0
        assignin('base','tmp',{'@(z)0'; ['@(z)' num2str(h.Bi)]});
    else
        assignin('base','tmp',{h.b1.bc.rhs.T.z{1}, ['@(z)' num2str(h.Bi)]});
    end
    evalin('base','b1.bc.rhs.T.z=tmp;');
    assignin('base','tmp',h.b1.bc.rhs.T.r);
    evalin('base','b1.bc.rhs.T.r=tmp;');
    
% grid parameters
    % changing from nondimensional inputs to dimensional outputs
        mesh.r.delta.start  = h.mesh.r.delta.start*h.l_lb/100;
        mesh.r.delta.fit = h.mesh.r.delta.fit*h.l_lb/100;
        mesh.r.delta.end = h.mesh.r.delta.end*h.l_lb/100;
        % -------------------------------------------------------------- %
        mesh.z.delta.start = h.mesh.z.delta.start*h.l_lb/100;
        mesh.z.delta.fit   = h.mesh.z.delta.fit*h.l_lb/100;
        mesh.z.delta.end   = h.mesh.z.delta.end*h.l_lb/100;
        
    % if linear distribution change some spacings to avoid errors
        if h.mesh.r.f(1) == 1
            mesh.r.delta.start = mesh.r.delta.fit;
        end
        if h.mesh.r.f(2) == 1
            mesh.r.delta.end   = mesh.r.delta.fit;
        end
        % --------------------------------------------------------- %
        if h.mesh.z.f(1) == 1
            mesh.z.delta.start = mesh.z.delta.fit;
        end
        if h.mesh.z.f(2) == 1
            mesh.z.delta.end   = mesh.z.delta.fit;
        end
        
    assignin('base','tmp',mesh.r.delta.start);
    evalin('base','mesh.r.delta_start=tmp;')
    assignin('base','tmp',mesh.r.delta.fit);
    evalin('base','mesh.r.delta_fit=tmp;')
    assignin('base','tmp',mesh.r.delta.end);
    evalin('base','mesh.r.delta_end=tmp;')
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.mesh.r.f(1));
    evalin('base','mesh.r.f_start=tmp;')
    assignin('base','tmp',h.mesh.r.f(2));
    evalin('base','mesh.r.f_end=tmp;')
    % ------------------------------------------------------------------ %
    assignin('base','tmp',[h.mesh.r.sf; h.mesh.r.sf]);
    evalin('base','mesh.r.sf=tmp;');
    % ------------------------------------------------------------------ %
    % ------------------------------------------------------------------ %
    assignin('base','tmp',mesh.z.delta.start);
    evalin('base','mesh.z.delta_start=tmp;')
    assignin('base','tmp',mesh.z.delta.fit);
    evalin('base','mesh.z.delta_fit=tmp;')
    assignin('base','tmp',mesh.z.delta.end);
    evalin('base','mesh.z.delta_end=tmp;')
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.mesh.z.f(1));
    evalin('base','mesh.z.f_start=tmp;')
    assignin('base','tmp',h.mesh.z.f(2));
    evalin('base','mesh.z.f_end=tmp;')
    % ------------------------------------------------------------------ %
    assignin('base','tmp',[h.mesh.z.sf; h.mesh.z.sf]);
    evalin('base','mesh.z.sf=tmp;');

% fluids
    assignin('base','selectedLiquid',h.selectedLiquid);
    assignin('base','selectedGas',h.selectedGas);
    
% optical ray
    assignin('base','r_0',h.r_0*h.r_i/h.l_lb)
    assignin('base','z_p',h.z_p+0.5)
    assignin('base','N_coeff',[h.N_coeff(1) h.N_coeff(2) h.N_coeff(3)])
    
% define colors for residual plot
    evalin('base','mycolor.blue = [0 0.4470 0.7410];')
    evalin('base','mycolor.red = [0.8500 0.3250 0.0980];')
    evalin('base','mycolor.ocher = [0.9290 0.6940 0.1250];')
    evalin('base','mycolor.violet = [0.4940 0.1840 0.5560];')
    evalin('base','mycolor.green = [0.4660 0.6740 0.1880];')
    evalin('base','h = zeros(1,5);')

% clear intermediate communication variables
    evalin('base','clear tmp');