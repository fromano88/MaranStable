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
assignin('base','r_o',h.r_o);
assignin('base','d_o',h.r_o*2);
% ------------------------------------------------------------------- %
assignin('base','l_lb',h.l_lb);
assignin('base','l_d1',h.l_d1);
assignin('base','l_d2',h.l_d2);
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
    assignin('base','tmp',h.flowopt.eigs.krylov);
    evalin('base','flowopt.eigs.krylov=tmp;');
    % ------------------------------------------------------------------- %
    assignin('base','tmp',h.flowopt.eigs.tol);
    evalin('base','flowopt.eigs.tol=tmp;');
    % ------------------------------------------------------------------- %
    assignin('base','tmp',h.flowopt.eigs.maxit);
    evalin('base','flowopt.eigs.maxit=tmp;');
    % ------------------------------------------------------------------- %
    assignin('base','tmp',h.flowopt.tolerance.growth);
    evalin('base','flowopt.tolerance.growth=tmp;');
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
    assignin('base','T_d1g',h.T_d1g);
    assignin('base','T_d2g',h.T_d2g);
    assignin('base','T_p1',h.T_p1);
    assignin('base','T_p2',h.T_p2);
    assignin('base','T_initial',(h.T_d1l+h.T_d2l)/2);
    assignin('base','delta_T',h.delta_T);
    
% (mean) inflow velocity
    assignin('base','w_in',h.w_in);

% boundary conditions - general
    if h.r_c == 0
        assignin('base','tmp',{h.b1.bc.z(1,1); h.b1.bc.z(2,:)});
    else
        assignin('base','tmp',{h.b1.bc.z(1,1:3); h.b1.bc.z(2,:)});
    end
    evalin('base','b1.bc.z=tmp;');
    
    assignin('base','tmp',{h.b1.bc.r(1,1:3); h.b1.bc.r(2,1:3)});
    evalin('base','b1.bc.r=tmp;');

    assignin('base','tmp',{h.b2.bc.z(1,1:3); h.b1.bc.z(2,:); h.b2.bc.z(3,1:3); h.b2.bc.z(4,1:3); h.b2.bc.z(4,1:3); h.b2.bc.z(4,1:3);});
    evalin('base','b2.bc.z=tmp;');
    
    assignin('base','tmp',{h.b2.bc.r(1,1:3); h.b2.bc.r(2,1:3)});
    evalin('base','b2.bc.r=tmp;');

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
    if h.r_c == 0
        assignin('base','tmp',{'@(z)0'; '@(z)0'});
    else
        assignin('base','tmp',h.b1.bc.rhs.T.z);
    end
    evalin('base','b1.bc.rhs.T.z=tmp;');
    
    assignin('base','tmp',h.b1.bc.rhs.T.r);
    evalin('base','b1.bc.rhs.T.r=tmp;');
    
    h.b2.bc.rhs.T.z{5} = h.b2.bc.rhs.T.z{4}; h.b2.bc.rhs.T.z{6} = h.b2.bc.rhs.T.z{4};
    assignin('base','tmp',h.b2.bc.rhs.T.z);
    evalin('base','b2.bc.rhs.T.z=tmp;');
    
    assignin('base','tmp',h.b2.bc.rhs.T.r);
    evalin('base','b2.bc.rhs.T.r=tmp;');
    
    assignin('base','tmp',h.b2.bc.rhs.w.r);
    evalin('base','b2.bc.rhs.w.r=tmp;');
    
% grid parameters
    % changing from nondimensional inputs to dimensional outputs
        mesh.r.delta.start  = h.mesh.r.delta.start*h.l_lb/100;
        mesh.r.delta.fit(1) = h.mesh.r.delta.fit(1)*h.l_lb/100;
        mesh.r.delta.fit(2) = h.mesh.r.delta.fit(2)*(h.r_o-h.r_i)/100;
        mesh.r.delta.end(1) = h.mesh.r.delta.end(1)*h.l_lb/100;
        mesh.r.delta.end(2) = h.mesh.r.delta.end(2)*(h.r_o-h.r_i)/100;
        % -------------------------------------------------------------- %
        mesh.z.delta.start(1)   = h.mesh.z.delta.start(1)*h.l_d1/100;
        mesh.z.delta.start(2:3) = h.mesh.z.delta.start(2:3)*h.l_lb/100;
        mesh.z.delta.fit(1)     = h.mesh.z.delta.fit(1)*h.l_d1/100;
        mesh.z.delta.fit(2)     = h.mesh.z.delta.fit(2)*h.l_lb/100;
        mesh.z.delta.fit(3)     = h.mesh.z.delta.fit(3)*h.l_d2/100;
        mesh.z.delta.end(1:2)   = h.mesh.z.delta.end(1:2)*h.l_lb/100;
        mesh.z.delta.end(3)     = h.mesh.z.delta.end(3)*h.l_d2/100;
        
    % if linear distribution change some spacings to avoid errors
        if h.mesh.r.f(1) == 1
            mesh.r.delta.start(1) = mesh.r.delta.fit(1);
        end
        if h.mesh.r.f(2) == 1
            mesh.r.delta.end(1)   = mesh.r.delta.fit(1);
        end
        if h.mesh.r.f(3) == 1
            mesh.r.delta.start(2) = mesh.r.delta.fit(2);
        end
        if h.mesh.r.f(4) == 1
            mesh.r.delta.end(2) = mesh.r.delta.fit(2);
        end
        % --------------------------------------------------------- %
        if h.mesh.z.f(1) == 1
            if h.mesh.z.f(2) == 1
                mesh.z.delta.start(1) = mesh.z.delta.fit(2);
            else
                mesh.z.delta.start(1) = mesh.z.delta.fit(1);
            end
        end
        if h.mesh.z.f(2) == 1
            mesh.z.delta.fit(1)   = mesh.z.delta.fit(2);
            mesh.z.delta.end(1)   = mesh.z.delta.fit(2);
            mesh.z.delta.start(2) = mesh.z.delta.fit(2);
        end
        if h.mesh.z.f(3) == 1
            mesh.z.delta.end(2)   = mesh.z.delta.fit(2);
            mesh.z.delta.start(3) = mesh.z.delta.fit(2);
            mesh.z.delta.fit(3)   = mesh.z.delta.fit(2);
        end
        if h.mesh.z.f(4) == 1
            if h.mesh.z.f(3) == 1
                mesh.z.delta.end(3) = mesh.z.delta.fit(2);
            else
                mesh.z.delta.end(3) = mesh.z.delta.fit(3);
            end
        end
        
    assignin('base','tmp',mesh.r.delta.start);
    evalin('base','mesh.r.delta_start=tmp;')
    assignin('base','tmp',mesh.r.delta.fit);
    evalin('base','mesh.r.delta_fit=tmp;')
    assignin('base','tmp',mesh.r.delta.end);
    evalin('base','mesh.r.delta_end=tmp;')
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.mesh.r.f([1 3]));
    evalin('base','mesh.r.f_start=tmp;')
    assignin('base','tmp',h.mesh.r.f([2 4]));
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
    assignin('base','tmp',h.mesh.z.f(1:3));
    evalin('base','mesh.z.f_start=tmp;')
    assignin('base','tmp',h.mesh.z.f(2:4));
    evalin('base','mesh.z.f_end=tmp;')
    % ------------------------------------------------------------------ %
    assignin('base','tmp',[h.mesh.z.sf; h.mesh.z.sf; h.mesh.z.sf]);
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