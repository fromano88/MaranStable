function SingleFluidModel_V3d2

evalin('base','clearvars'); close all; clc

f_main = figure('units','pixels','menubar','none','name',...
    'SINGLE FLUID MODEL','numbertitle','off','resize','off','Visible','Off');

if ispc % windows
    f_size = [1100 680];
else    % linux
    f_size = [1131 760];
end
f_main.Position(3:4) = f_size;
movegui(f_main,'center')

% delete user defined fluids
if isfile([tempdir 'my_fluids.mat'])
    delete([tempdir 'my_fluids.mat'])
end

defaultValues_SFM

createFig_SFM

% default sketch
    updateSketch_SFM
    
% default text
    defaultText_SFM_V3d2
    
% logo
    axes(h.as.logo);
    logo_MS(0)

set(f_main,'WindowButtonMotionFcn',{@ce_button,h})
set([h.tabs{:}],'callback',{@ce_changeMenu,h})
set([h.editTextBoxes{:}],'call',{@ce_getEditValue,h})
set(h.gl.bg.formulation,'SelectionChangedFcn',{@gl_bg_formulation,h})
set(h.gl.bg.model,'SelectionChangedFcn',{@gl_bg_model,h})
set(h.gl.pm.liquid,'call',{@gl_pm_fluid,h})
set(h.gl.pb.create,'call',{@gl_pb_create,h})
set(h.mh.pm.headline,'call',{@mh_pm_headline,h})
set([h.mh.buttonGroups.radial{:}],'SelectionChangedFcn',{@mh_bg_radial,h})
set([h.mh.buttonGroups.axial{:}],'SelectionChangedFcn',{@mh_bg_axial,h})
set([h.mh.pb.checkMesh h.mh.pb.drawMesh],'call',{@mh_pbs,h})
set(h.es.cb.energy,'call',{@es_cb_energy,h})
set([h.es.cb.stokes h.es.cb.marangoni],'call',{@es_cbs,h})
set(h.bc.pm.headline,'call',{@bc_pm_headline,h})
set([h.bc.buttonGroups.velocity{:}],'SelectionChangedFcn',{@bc_bg_velocity,h})
set([h.bc.buttonGroups.temperature{:}],'SelectionChangedFcn',{@bc_bg_temperature,h})
set(h.bc.bg.radialLB.rc(1),'SelectionChangedFcn',{@bc_bg_axis,h})
set(h.bc.bg.radialLB.ri(1),'SelectionChangedFcn',{@bc_bg_surface,h})
set([h.bc.popupMenus.temperature{:}],'call',{@bc_pm_temperature,h})
set([h.bc.et.axialLB.d1.temperature h.bc.et.axialLB.d2.temperature h.bc.et.radialLB.rc.temperature],'ButtonDownFcn',{@bc_pm_temperature,h})
set(h.sn.bg.initialization,'SelectionChangedFcn',{@sn_bg_initialization,h})
set(h.sn.pb.run,'call',{@sn_pb_run,h})
set([h.sn.pb.stop h.sn.pb.load h.sn.pb.save h.sn.pb.plot h.sy.pb.load h.sy.pb.save h.sy.pb.plot],'call',{@sn_pbs,h})
set(h.or.pb.run,'call',{@or_pb_run,h})
set(h.sy.bg.parameter,'SelectionChangedFcn',{@sy_bg_parameter,h})
set(h.sy.pb.run(1),'call',{@sy_pb_most_d_mode,h})
set(h.sy.pb.run(2),'call',{@sy_pb_critical_mode,h})
set(h.sy.pb.stop,'call',{@sy_pb_stop,h})

% correct figure size
    f_main.Position(3:4) = f_size;
    set(f_main,'Visible','On')

guidata(f_main,h)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%                  ____    ____     ____                                  %
%                 |       |    |   |    |    ___                          %
%                 |       |    |   |  __|   |___                          %
%                 |____   |____|   |  \     |___                          %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ce_button(hObject, ~, ~)
h = guidata(hObject);

pos = get(hObject, 'currentpoint'); % get mouse location on figure
x = pos(1); y = pos(2);

if isunix && x > 1091 && x < 1125 && y > 730 && y < 756
    set(h.pl.info(1),'visible','on')
elseif ispc && x > 1060 && x < 1094 && y > 653 && y < 678
    set(h.pl.info(1),'visible','on')
else
    set(h.pl.info(1),'visible','off')
end

function ce_changeMenu(hObject, ~, ~)
h = guidata(hObject);

% checking if b.c.s are well defined
boundaryList = {h.b1.bc.rhs.T.z{1} h.b1.bc.rhs.T.r{1} h.b1.bc.rhs.T.r{2}};
if max(cellfun(@isempty,boundaryList)) == 1 && hObject ~= h.pb.bcs
    waitfor(warndlg({'Boundary condition not defined!','', 'Set boundary condition before going on!'},'Warning 158'))
    return
end

set([h.mainPanels{:}],'visible','off')
set([h.tabs{:}],'ForegroundColor',h.unselectedTabColor)

switch hObject
    case h.pb.general
        set(h.pl.general,'visible','on')
        set(h.pb.general,'ForegroundColor',h.selectedTabColor)
        
    case h.pb.geometry
        set(h.pl.geometry,'visible','on')
        set(h.pb.geometry,'ForegroundColor',h.selectedTabColor)
        
    case h.pb.mesh
        set(h.pl.mesh,'visible','on')
        set(h.pb.mesh,'ForegroundColor',h.selectedTabColor)
        set([h.mh.pl.radialLB.main h.mh.pl.axialLB.main],'visible','off')
        if strcmp(h.mh.Visible,'radial')
            set(h.mh.pm.headline,'value',1)
            set(h.mh.pl.radialLB.main,'visible','on')
        else
            set(h.mh.pm.headline,'value',2)
            set(h.mh.pl.axialLB.main,'visible','on')
        end
        
    case h.pb.equations
        set(h.pl.equations,'visible','on')
        set(h.pb.equations,'ForegroundColor',h.selectedTabColor)
        
    case h.pb.bcs
        set(h.pl.bcs,'visible','on')
        set(h.pb.bcs,'ForegroundColor',h.selectedTabColor)
        set([h.bc.pl.radialLB.main h.bc.pl.axialLB.main],'visible','off')
        if strcmp(h.bc.Visible,'radial')
            set(h.bc.pm.headline,'value',1)
            set(h.bc.pl.radialLB.main,'visible','on')
        else
            set(h.bc.pm.headline,'value',2)
            set(h.bc.pl.axialLB.main,'visible','on')
        end
        
    case h.pb.simulation
        set(h.pl.simulation,'visible','on')
        set(h.pb.simulation,'ForegroundColor',h.selectedTabColor)
        
    case h.pb.ray
        set(h.pl.ray,'visible','on')
        set(h.pb.ray,'ForegroundColor',h.selectedTabColor)
        
    case h.pb.lsa
        set(h.pl.lsa,'visible','on')
        set(h.pb.lsa,'ForegroundColor',h.selectedTabColor)
        
    case h.pb.about
        set(h.pl.about,'visible','on')
        set(h.pb.about,'ForegroundColor',h.selectedTabColor)
end

guidata(hObject,h)

function ce_getEditValue(hObject, ~, ~)
h = guidata(hObject);

p = h.flowopt.ax + 1; rx = 'xr';

% define limits for the values
switch hObject
    case h.gl.et.g
        h.lowerLimit=-inf; h.upperLimit=inf; name = 'g';
        callbackName = 'gl_et_g'; valueName = 'flowopt.g';

    case h.gl.et.Vr
        h.lowerLimit=0; h.upperLimit=inf; name = 'V_r';
        errormsg = 'Value must be positive.'; %#ok<NASGU> 
        callbackName = 'gl_et_Vr';

    case h.gy.et.llb
        h.lowerLimit=-inf; h.upperLimit=inf; name = 'l_lb';
        errormsg = 'Value must be positive.'; default = h.l_lb*1000; %#ok<NASGU> 
        callbackName = 'gy_et_d';

    case h.gy.et.rc
        h.lowerLimit=0; h.upperLimit=h.r_i*1000; name = 'r_c';
        callbackName = 'gy_et_rc'; opts = h.opts; default = h.r_c*1000; %#ok<NASGU> 
        errormsg = ['Value must be greater or equal zero and less than ' rx(p) '_i.']; %#ok<NASGU> 

    case h.gy.et.ri
        h.lowerLimit=h.r_c*1000; h.upperLimit=inf; name = 'r_i';
        callbackName = 'gy_et_ri'; opts = h.opts; default = h.r_i*1000;  %#ok<NASGU> 
        errormsg = ['Value must be greater than ' rx(p) '_c and less than ' rx(p) '_o.']; %#ok<NASGU> 

    case {h.mh.et.radialLB.spacing(1),h.mh.et.radialLB.spacing(2),...
            h.mh.et.radialLB.spacing(3),h.mh.et.axialLB.spacing(1),...
            h.mh.et.axialLB.spacing(2),h.mh.et.axialLB.spacing(3)}
        h.lowerLimit=1e-14; h.upperLimit=100; name = 'spacing';
        switch hObject
            case h.mh.et.radialLB.spacing(1)
                valueName = 'mesh.r.delta.start';
            case h.mh.et.radialLB.spacing(2)
                valueName = 'mesh.r.delta.fit';
            case h.mh.et.radialLB.spacing(3)
                valueName = 'mesh.r.delta.end';
            case h.mh.et.axialLB.spacing(1)
                valueName = 'mesh.z.delta.start';
            case h.mh.et.axialLB.spacing(2)
                valueName = 'mesh.z.delta.fit';
            case h.mh.et.axialLB.spacing(3)
                valueName = 'mesh.z.delta.end';
        end

    case {h.mh.et.radialLB.stretchingC(1),h.mh.et.radialLB.stretchingC(3),...
            h.mh.et.axialLB.stretchingC(1),h.mh.et.axialLB.stretchingC(3)}
        h.lowerLimit=1; h.upperLimit=inf; name = 'stretchingC';
        callbackName = 'mh_et_stretchingC';
        errormsg = 'Value must be greater than 1.'; %#ok<NASGU> 
        switch hObject
            case h.mh.et.radialLB.stretchingC(1)
                valueName = 'mesh.r.f(1)';
            case h.mh.et.radialLB.stretchingC(3)
                valueName = 'mesh.r.f(2)';
            case h.mh.et.axialLB.stretchingC(1)
                valueName = 'mesh.z.f(1)';
            case h.mh.et.axialLB.stretchingC(3)
                valueName = 'mesh.z.f(2)';
        end
        strchngC = ['h.' valueName];
        
    case {h.bc.et.axialLB.d1.temperature h.bc.et.axialLB.d2.temperature ...
            h.bc.et.radialLB.rc.temperature h.bc.et.radialLB.ri.Bi}
        h.lowerLimit=-inf; h.upperLimit=inf; name = 'temperature';
        callbackName = 'bc_et_temperature';
        switch hObject
            case h.bc.et.axialLB.d1.temperature
                valueName = 'b1.bc.rhs.T.r{1}';
            case h.bc.et.axialLB.d2.temperature
                valueName = 'b1.bc.rhs.T.r{2}';
            case h.bc.et.radialLB.rc.temperature
                valueName = 'b1.bc.rhs.T.z{1}';
            case h.bc.et.radialLB.ri.Bi
                valueName = 'b1.bc.rhs.T.z{2}';
        end
        h.vlNm = ['h.' valueName];
        default = eval(h.vlNm);
        if ~isempty(default)
            default = str2double(default(5:end)); %#ok<NASGU> 
        else
            default = 25; %#ok<NASGU> 
        end
        
    case {h.sn.et.residuals(1),h.sn.et.residuals(2)}
        h.lowerLimit=1.0e-15; h.upperLimit=inf; name = 'residuals';
        errormsg = 'Value must be greater or equal 1.0e-15!'; %#ok<NASGU> 
        switch hObject
            case h.sn.et.residuals(1)
                valueName = 'flowopt.tolerance.residuals';
            case h.sn.et.residuals(2)
                valueName = 'flowopt.tolerance.newton';
        end
        
    case h.sn.et.steps
        h.lowerLimit=-inf; h.upperLimit=inf; name = 'steps';
        callbackName = 'sn_et_steps';
        
    case h.or.et.ray
        h.lowerLimit=0; h.upperLimit=1; name = 'r_0';
        callbackName = 'or_et_ray';
        
    case h.or.et.particle
        h.lowerLimit=-0.5; h.upperLimit=0.5; name = 'z_p';
        callbackName = 'or_et_particle';

    case {h.or.et.N(1),h.or.et.N(2),h.or.et.N(3)}
        h.lowerLimit=-inf; h.upperLimit=inf; name = 'N';
        callbackName = 'or_et_N';
        switch hObject
            case h.or.et.N(1)
                valueName = 'N_coeff(1)';
            case h.or.et.N(2)
                valueName = 'N_coeff(2)';
            case h.or.et.N(3)
                valueName = 'N_coeff(3)';
        end
        
    case {h.sy.et.n_eig,h.sy.et.n_eig_c,h.sy.et.kryl,h.sy.et.maxit}
        h.lowerLimit=-inf; h.upperLimit=inf; name = 'n_eig';
        callbackName = 'sy_et_eigs';
        switch hObject
            case h.sy.et.n_eig
                valueName = 'flowopt.eigs.n';
            case h.sy.et.n_eig_c
                valueName = 'flowopt.eigs.n_cayley';
            case h.sy.et.kryl
                valueName = 'flowopt.eigs.krylov';
            case h.sy.et.maxit
                valueName = 'flowopt.eigs.maxit';
        end

    case {h.sy.et.tol_eigs,h.sy.et.conv}
        h.lowerLimit=1.0e-15; h.upperLimit=inf; name = 'tolerance';
        errormsg = 'Value must be greater or equal 1.0e-15!'; %#ok<NASGU> 
        switch hObject
            case h.sy.et.tol_eigs
                valueName = 'flowopt.eigs.tol';
            case h.sy.et.conv
                valueName = 'flowopt.tolerance.growth';
        end
        
    case {h.sy.et.m_start,h.sy.et.m_delta,h.sy.et.m_end}
        h.lowerLimit=-inf; h.upperLimit=inf;
        callbackName = 'sy_et_m'; opts = h.opts; %#ok<NASGU>
        switch hObject
            case h.sy.et.m_start
                name = 'm_start';
            case h.sy.et.m_delta
                name = 'm_delta';
            case h.sy.et.m_end
                name = 'm_end';
        end
        
end

if ~exist('valueName','var')
    valueName = name;
end
if ~exist('default','var')
    eval(['default = h.' valueName ';']) %#ok<EVLDOT> 
end
showOption = h.show.(name); %#ok<NASGU> 
if ~exist('errormsg','var')
    errormsg = ['Choose a value between ' num2str(h.lowerLimit) ' and ' num2str(h.upperLimit) '.']; %#ok<NASGU> 
end
if ~exist('opts','var')
    opts = struct('WindowStyle','non-modal','Interpreter','none'); %#ok<NASGU> 
end

% get value
    eval(['[h.' valueName ', showMessage] = getValue(h.lowerLimit,h.upperLimit,default,hObject,errormsg,showOption,opts);']) %#ok<EVLDOT> 

% update showoptions
    if isempty(showMessage.error) == 0
        h.show.(name).error = showMessage.error;
    end
    if isempty(showMessage.warning) == 0
        h.show.(name).warning = showMessage.warning;
    end
    
if exist('strchngC','var')
    h.strchngC = eval(strchngC); %#ok<EVLDOT> 
end

guidata(hObject,h)

if exist('callbackName','var')
    eval([callbackName '(hObject, h)'])
end

function h = ce_updateFunctions(h)

% update every temperature profile in case the profile contains this geometry parameter
    if strcmp(h.b1.bc.z(1,3),'c')
        h.b1.bc.rhs.T.z{1} = get(h.bc.et.radialLB.rc.temperature,'string');
        h.b1.bc.rhs.T.z{1} = convertTemperature(h,h.b1.bc.rhs.T.z{1},'z');
    end
    % ------------------------------------------------------------------- %
    if strcmp(h.b1.bc.r(1,3),'c')
        h.b1.bc.rhs.T.r{1} = get(h.bc.et.axialLB.d1.temperature,'string');
        h.b1.bc.rhs.T.r{1} = convertTemperature(h,h.b1.bc.rhs.T.r{1},'r');
    end
    % --------------------------------------------------------------------%
    if strcmp(h.b1.bc.r(2,3),'c')
        h.b1.bc.rhs.T.r{2} = get(h.bc.et.axialLB.d2.temperature,'string');
        h.b1.bc.rhs.T.r{2} = convertTemperature(h,h.b1.bc.rhs.T.r{2},'r');
    end

function ce_basic_state(h)

axes(h.sy.as.eigs); cla(h.sy.as.eigs,'reset');
ax = gca; ax.YScale = 'log'; ax.XLim = [1 2]; ax.XTick = [1 2]; ax.Box = 'on'; xlabel('Iterations'); ylabel('Residuals');
ce_display_dependent(h,0); axes(h.sy.as.eigs)
updateValues_SFM; evalin('base','caseSetup_SFM_V3d2')
evalin('base','flowopt.stability=0; grid_generation; block_assembly_V3d2');

try
    for iteration = 1:h.steps
        assignin('base','iteration',iteration)
        evalin('base','dx = J\(A*x-b);')
        evalin('base','x = x - dx;')
        evalin('base','block_assembly_V3d2')
        evalin('base','plotResiduals')
        residual = evalin('base', 'residual');
        dx = evalin('base', 'dx');
        if (residual <= h.flowopt.tolerance.residuals && norm(dx) <= h.flowopt.tolerance.newton) || isnan(residual) || residual>10^15  || get(h.sy.pb.stop,'userdata'), break, end
    end
catch
    waitfor(errordlg('Simulation stopped!','Error 201'))
    return
end
evalin('base','clear A J b convergence')
evalin('base','old.x=x;')

if residual > 10^15 || isnan(residual)
    errordlg('Simulation stopped!','Error 201')
elseif get(h.sy.pb.stop,'userdata')
    set(h.sy.st.run,'String','running ...','visible','off')
    set(h.sy.pb.run(2),'Visible','on')
    set(h.sy.pb.stop,'Visible','off')
    helpdlg('Simulation stopped!','Notice')
    set([h.tabs{:} h.sy.et.n_eig h.sy.et.n_eig_c h.sy.et.m_start h.sy.et.m_delta h.sy.et.m_end h.sy.rb.delta_T h.sy.rb.T_d1 h.sy.rb.T_d2 h.sy.pb.run(1) h.sy.pb.run(2) h.sy.pb.load h.sy.pb.save h.sy.pb.plot],'Enable','on')
end

function h = sy_update_parameter(h)

bc_list = {'h.b1.bc.rhs.T.r{1}', 'h.b1.bc.rhs.T.r{2}', 'h.b1.bc.rhs.T.z{1}'};

switch h.dependent
    case 'delta_T'
        h.T_0 = (h.T_d1l+h.T_d2l)/2;
        if h.T_d1l>h.T_d2l
            h.T_d1l = h.T_0 + h.delta_T/2;
            h.T_d2l = h.T_0 - h.delta_T/2;
        else
            h.T_d1l = h.T_0 - h.delta_T/2;
            h.T_d2l = h.T_0 + h.delta_T/2;
        end
        for i = 1:length(h.equal_T_d1)
            h.(h.equal_T_d1{i}) = h.T_d1l;
            eval([bc_list{h.idx_T_d1(i)} ' = ''@(' bc_list{h.idx_T_d1(i)}(end-3) ')' num2str(h.T_d1l) ''';'])
        end
        for i = 1:length(h.equal_T_d2)
            h.(h.equal_T_d2{i}) = h.T_d2l;
            eval([bc_list{h.idx_T_d2(i)} ' = ''@(' bc_list{h.idx_T_d2(i)}(end-3) ')' num2str(h.T_d2l) ''';'])
        end

    case 'T_d1l'
        for i = 1:length(h.equal_T_d1)
            h.(h.equal_T_d1{i}) = h.T_d1l;
            eval([bc_list{h.idx_T_d1(i)} ' = ''@(' bc_list{h.idx_T_d1(i)}(end-3) ')' num2str(h.T_d1l) ''';'])
        end

    case 'T_d2l'
        for i = 1:length(h.equal_T_d2)
            h.(h.equal_T_d2{i}) = h.T_d2l;
            eval([bc_list{h.idx_T_d2(i)} ' = ''@(' bc_list{h.idx_T_d2(i)}(end-3) ')' num2str(h.T_d2l) ''';'])
        end

end


function ce_display_dependent(h,disp_g)

if nargin == 1, disp_g = 1; end

axes(h.sy.as.dependent); cla(h.sy.as.dependent); %#ok<*LAXES>
switch h.dependent
    case 'delta_T'
        if disp_g==1 && isfield(h,'gamma')
            text(0,0,['$\Delta T = ' num2str(h.delta_T) '\,^\circ \mathrm{C} \; \to \;\gamma_c=' num2str(-real(h.gamma(1))) '\pm ' num2str(abs(imag(h.gamma(1)))) '\,i$ [1/s]'],h.latex_r{:});
        else
            text(0,0,['$\Delta T = ' num2str(h.delta_T) '\,^\circ$C'],h.latex_r{:});
        end

    case 'T_d1l'
        if disp_g==1 && isfield(h,'gamma')
            text(0,0,['$T_{d1} = ' num2str(h.T_d1l) '\,^\circ \mathrm{C} \; \to \;\gamma_c=' num2str(-real(h.gamma(1))) '\pm ' num2str(abs(imag(h.gamma(1)))) '\,i$ [1/s]'],h.latex_r{:});
        else
            text(0,0,['$T_{d1} = ' num2str(h.T_d1l) '\,^\circ$C'],h.latex_r{:});
        end

    case 'T_d2l'
        if disp_g==1 && isfield(h,'gamma')
            text(0,0,['$T_{d2} = ' num2str(h.T_d2l) '\,^\circ \mathrm{C} \; \to \;\gamma_c=' num2str(-real(h.gamma(1))) '\pm ' num2str(abs(imag(h.gamma(1)))) '\,i$ [1/s]'],h.latex_r{:});
        else
            text(0,0,['$T_{d2} = ' num2str(h.T_d2l) '\,^\circ$C'],h.latex_r{:});
        end

end
% make a little break such that user can see the new value of gamma before
% the new basic state is computed
if disp_g==1 && isfield(h,'gamma')
    pause(7)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%           ____                           ____                           %
%          |         ___    |\  |   ___   |    |    /\    |               %
%          |  __    |___    | \ |  |___   |  __|   /__\   |               %
%          |____|   |___    |  \|  |___   |  \    /    \  |___            %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gl_bg_formulation(hObject, ~, ~)
h = guidata(hObject);

selection = get(get(hObject,'SelectedObject'),'String');

switch selection
    case 'Planar'
        h.flowopt.ax=0;
        set(h.sy.et.m_delta,'enable','on')
        if ispc
            txt_wavenumber = 'Wave Numbers along z';
        else
            txt_wavenumber = 'Wave Numbers in z';
        end

        h.b1.bc.z(1,1) = 'w';
        set([h.bc.rb.radialLB.rc.noSlip h.bc.rb.radialLB.rc.slip],'Enable','on')
        set(h.bc.rb.radialLB.rc.wall,'Value',1)
        if strcmp(h.b1.bc.z(1,2),'s')
            set(h.bc.rb.radialLB.rc.slip,'Value',1)
            h.b1.bc.z(1,2) = 's';
        else
            set(h.bc.rb.radialLB.rc.noSlip,'Value',1)
            h.b1.bc.z(1,2) = 'n';
        end
    
        % conductive or adiabatic if enery equation turned on
        if get(h.es.cb.energy,'value')
            h.b1.bc.z(1,3) = 'c';
            set(h.bc.rb.radialLB.rc.adiabatic,'Value',0,'Enable','on')
            set(h.bc.rb.radialLB.rc.conductive,'Value',1,'Enable','on')
            set([h.bc.et.radialLB.rc.temperature h.bc.pm.radialLB.rc.temperature h.bc.st.radialLB.rc.temperature],'visible','on')
        else
            h.b1.bc.z(1,3) = 'a';
            set([h.bc.rb.radialLB.rc.adiabatic h.bc.rb.radialLB.rc.conductive],'Value',0,'Enable','off')
        end

    case 'Axisymmetric'
        h.flowopt.ax=1;
        h.m_start = 0; h.m_delta = 1; h.m_end = 1;
        set(h.sy.et.m_start,'string',num2str(h.m_start))
        set(h.sy.et.m_delta,'string',num2str(h.m_delta),'enable','off')
        set(h.sy.et.m_end,'string',num2str(h.m_end))
        if ispc
            txt_wavenumber = 'Wave Numbers along φ';
        else
            txt_wavenumber = 'Wave Numbers in φ';
        end

        if h.r_c == 0
            h.b1.bc.z(1,1:3) = 'a  ';
            % no further b.c.s can be imposed at r=0 except symmetry condition
            set(h.bc.rb.radialLB.rc.axis,'Value',1)
            set([h.bc.rb.radialLB.rc.noSlip h.bc.rb.radialLB.rc.slip h.bc.rb.radialLB.rc.adiabatic h.bc.rb.radialLB.rc.conductive],'Value',0,'Enable','off')
        
             % cancel visibility of imposed temperature profile at r=0
            cla(h.bc.as.radialLB.rc(2),'reset');
            set([h.bc.as.radialLB.rc(2) h.bc.et.radialLB.rc.temperature h.bc.pm.radialLB.rc.temperature h.bc.st.radialLB.rc.temperature],'visible','off')
        end
end

set(h.sy.pl.m,'title',txt_wavenumber);

h = ce_updateFunctions(h);

guidata(hObject, h);

updateSketch_SFM
updateText_SFM_V3d2

function gl_et_g(hObject, ~, ~) %#ok<*DEFNU>
h = guidata(hObject);

h.flowopt.g = h.flowopt.g*9.81;
    
guidata(hObject,h)

% update sketch
updated = 'g'; %#ok<NASGU>
updateSketch_SFM

function gl_et_Vr(hObject, ~, ~) %#ok<*DEFNU>
h = guidata(hObject);

% changing value if V_r <= 0
if str2double(get(hObject,'string')) == 0
    if h.show.V_r.warning == 1
        waitfor(warndlg('Value must be positive.','Attention'))
        h.show.V_r.warning = 0;
    end
    h.V_r = 0.1;
    set(hObject,'string',num2str(h.V_r));
end

guidata(hObject,h)

updateSketch_SFM

function gl_bg_model(hObject, ~, ~)
h = guidata(hObject);

selection = get(get(hObject,'SelectedObject'),'String');

switch selection
    case 'Oberbeck-Boussinesq Approximation'
        if h.show.NS1 == 1
            waitfor(helpdlg({'All fluid properties are assumed to be constant at the mean reference temperature except for the density in the buoyancy term.'},'Notice'))
            h.show.NS1 = 0;
        end
        h.NS = 1; h.flowopt.boussinesq = 1;

    case 'Linearly Temperature Dependent'
        if h.show.NS2 == 1
            waitfor(helpdlg({'All fluid properties are linearized around the mean reference temperature.'},'Notice'))
            h.show.NS2 = 0;
        end
        h.NS = 2; h.flowopt.boussinesq = 0;

    case 'Fully Temperature Dependent'
        h.NS = 3; h.flowopt.boussinesq = 0;
        
end

guidata(hObject,h)

updateText_SFM_V3d2

function gl_pm_fluid(hObject, ~, ~)
h = guidata(hObject);

str=get(hObject,'String'); val=get(hObject,'Value');

switch str{val}
    case 'Others ...'
        selectedFluid = getFluid('liquid',h);
        if ~isempty(selectedFluid)
            h.selectedLiquid = selectedFluid;
        else
            selectedFluid = h.selectedLiquid;
        end
        % if the selected Liquid is not part of the popup menu, then it will be added to the list
        if max(contains(str,selectedFluid))
            idx = find(ismember(str,selectedFluid));
            set(hObject,'Value',idx)                
        else
            str=[str(1:end-1); selectedFluid; str(end)];
            set(hObject,'String',str);
        end
        
    otherwise
        h.selectedLiquid = str{val};
end

if h.flowopt.energy == 1
    set(h.gl.rb.NS3,'enable','on')
    if isfile([tempdir 'my_fluids.mat'])
        myfluids = who('-file',[tempdir 'my_fluids.mat']);
        % if user selects a customized fluid --> switch off NS3
        if max(ismember(myfluids,h.selectedLiquid))
            if h.NS == 3
                set(h.gl.rb.NS2,'value',1)
                h.NS = 2;
            end
            set(h.gl.rb.NS3,'enable','off')
        end
    end
end
    
guidata(hObject,h)

function gl_pb_create(hObject, ~, ~)
h = guidata(hObject);

[name, prop] = createFluid('SFM',h); %#ok<ASGLU> 
if isempty(name)
    return
end

h.liquid_list = {name, h.liquid_list{1:end}};
h.selectedLiquid = name;
str = get(h.gl.pm.liquid,'String');
str = [name; str];
set(h.gl.pm.liquid,'String',str,'Value',1)

eval([name '= prop;'])
if isfile([tempdir 'my_fluids.mat'])
    save([tempdir 'my_fluids.mat'],name,'-append')
else
    save([tempdir 'my_fluids.mat'],name)
end

if h.NS == 3 && h.flowopt.energy == 1
    waitfor(warndlg('Customoized fluids are only available for the Oberbeck-Boussinesq Approximation and for the linearized temperature model.','Attention'))
    set(h.gl.rb.NS2,'value',1)
    h.NS = 2;
end
set(h.gl.rb.NS3,'enable','off')

guidata(hObject,h)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%        ____             ___                 _____    ____               %
%       |         ___    |   |   |\  /|   ___   |     |    |   \  /       %
%       |  __    |___    |   |   | \/ |  |___   |     |  __|    \/        %
%       |____|   |___    |___|   |    |  |___   |     |  \      /        %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gy_et_d(hObject, ~, ~)
h = guidata(hObject);

h.l_lb = h.l_lb/1000;

% changing value if l <= 0
if str2double(get(hObject,'string')) <= 0
    if h.show.l_lb.warning == 1
        waitfor(warndlg('Value must be positive.','Attention'))
        h.show.l_lb.warning = 0;
    end
    h.l_lb = 0.1; set(hObject,'string',num2str(vah.l_lbl));
end

h = ce_updateFunctions(h);

guidata(hObject,h)

% update sketch
updated = 'l_lb'; %#ok<NASGU>
updateSketch_SFM

function gy_et_rc(hObject, ~, ~)
h = guidata(hObject);

h.r_c = h.r_c/1000;

% changing value if r_c == r_i
if h.r_c == h.r_i
    if h.r_i > 0.1/1000
        h.r_c = h.r_c-0.1/1000;
    else
        h.r_c = 0;
        h.r_i = 0.1/1000;
    end
    set(hObject,'string',num2str(h.r_c*1000));
    set(h.gy.et.ri,'String',num2str(h.r_i*1000));
end
    
% changing from anular to pure liquid bridge and vice versa
if h.r_c == 0
    if h.flowopt.ax == 1
        h.b1.bc.z(1,1:3) = 'a  ';
         % no further b.c.s can be imposed at r=0 except symmetry condition
        set(h.bc.rb.radialLB.rc.axis,'Value',1)
        set([h.bc.rb.radialLB.rc.noSlip h.bc.rb.radialLB.rc.slip h.bc.rb.radialLB.rc.adiabatic h.bc.rb.radialLB.rc.conductive],'Value',0,'Enable','off')
    
         % cancel visibility of imposed temperature profile at r=0
        cla(h.bc.as.radialLB.rc(2),'reset');
        set([h.bc.as.radialLB.rc(2) h.bc.et.radialLB.rc.temperature h.bc.pm.radialLB.rc.temperature h.bc.st.radialLB.rc.temperature],'visible','off')

    end
else
    h.b1.bc.z(1,1) = 'w';
    % wall at r=r_c
        set([h.bc.rb.radialLB.rc.noSlip h.bc.rb.radialLB.rc.slip],'Enable','on')
        if strcmp(h.b1.bc.z(1,2),'s')
            set(h.bc.rb.radialLB.rc.slip,'Value',1)
            h.b1.bc.z(1,2) = 's';
        else
            set(h.bc.rb.radialLB.rc.noSlip,'Value',1)
            h.b1.bc.z(1,2) = 'n';
        end
        set(h.bc.rb.radialLB.rc.wall,'Value',1)

    % conductive or adiabatic if enery equation turned on
        if get(h.es.cb.energy,'value')
            h.b1.bc.z(1,3) = 'c';
            set(h.bc.rb.radialLB.rc.adiabatic,'Value',0,'Enable','on')
            set(h.bc.rb.radialLB.rc.conductive,'Value',1,'Enable','on')
        else
            h.b1.bc.z(1,3) = 'a';
            set([h.bc.rb.radialLB.rc.adiabatic h.bc.rb.radialLB.rc.conductive],'Value',0,'Enable','off')
        end

    % make imposed temperature profile visible if energy eq. turned on and if conductive
        if get(h.es.cb.energy,'value') && strcmp(h.b1.bc.z(1,3),'c')
            set([h.bc.et.radialLB.rc.temperature h.bc.pm.radialLB.rc.temperature h.bc.st.radialLB.rc.temperature],'visible','on')
        end
end
    
h = ce_updateFunctions(h);

guidata(hObject,h)

% update sketch
if h.r_c ~= 0
    updated = 'r_c'; %#ok<NASGU>
end
updateSketch_SFM
updateText_SFM_V3d2

function gy_et_ri(hObject, ~, ~)
h = guidata(hObject);

h.r_i = h.r_i/1000;

% changing value if r_c == r_i
if h.r_c == h.r_i
    if h.r_c < h.upperLimit - 0.1/1000
        h.r_i = h.r_i + 0.1/1000;
    else
        h.r_i = h.r_c + 0.1/1000;
        h.r_o = h.r_i + 0.1/1000;
    end
    set(hObject,'string',num2str(h.r_i*1000,16));
end

h = ce_updateFunctions(h);

guidata(hObject,h)

% update sketch
updated = 'r_i'; %#ok<NASGU>
updateSketch_SFM


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%                                                                         %
%                      |\  /|   ___   ____                                %
%                      | \/ |  |___  |____   |__|                         %
%                      |    |  |___   ____|  |  |                         %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mh_pm_headline(hObject, ~, ~)
h = guidata(hObject);

% here use value in switch function because strings are changing in terms
% of planar or axisymmetric
val=get(hObject, 'Value');

set([h.mh.pl.radialLB.main h.mh.pl.axialLB.main],'visible','off')

switch val
    case 1
        set(h.mh.pl.radialLB.main, 'visible', 'on')
        h.mh.Visible = 'radial';
    case 2
        set(h.mh.pl.axialLB.main, 'visible', 'on')
        h.mh.Visible = 'axial';
end

guidata(hObject,h)

function mh_et_stretchingC(hObject, ~, ~)
h = guidata(hObject);

% if stretching factor == 1 --> linear distribution
% if stretching factor  > 1 --> distribution function of corresponding panels
radial_gp = get(h.mh.rb.radialLB(1).gp,'value') + get(h.mh.rb.radialLB(3).gp,'value');
axial_gp = get(h.mh.rb.axialLB(1).gp,'value') + get(h.mh.rb.axialLB(3).gp,'value');
msg_txt = 'Corresponding Stretching Function f set to LINEAR!';
if h.strchngC == 1
    if h.show.stretching == 1
        waitfor(helpdlg(msg_txt,'Notice'))
        h.show.stretching = 0;
    end
end
switch hObject
    case h.mh.et.radialLB.stretchingC(1)
        if h.mesh.r.f(1) == 1
            set(h.mh.rb.radialLB(1).lin,'Value',1);
            set(h.mh.et.radialLB.spacing(1),'enable','off')
        else
            if radial_gp > 0
                set(h.mh.rb.radialLB(1).gp,'Value',1);
            else
                set(h.mh.rb.radialLB(1).tanh,'Value',1);
            end
            set(h.mh.et.radialLB.spacing(1),'enable','on')
        end
        
    case h.mh.et.radialLB.stretchingC(3)
        if h.mesh.r.f(2) == 1
            set(h.mh.rb.radialLB(3).lin,'Value',1);
            set(h.mh.et.radialLB.spacing(3),'enable','off')
        else
            if radial_gp > 0
                set(h.mh.rb.radialLB(3).gp,'Value',1);
            else
                set(h.mh.rb.radialLB(3).tanh,'Value',1);
            end
            set(h.mh.et.radialLB.spacing(3),'enable','on')
        end
        
    case h.mh.et.axialLB.stretchingC(1)
        if h.mesh.z.f(1) == 1
            set(h.mh.rb.axialLB(1).lin,'Value',1);
            set(h.mh.et.axialLB.spacing(1),'enable','off')
        else
            if axial_gp > 0
                set(h.mh.rb.axialLB(1).gp,'Value',1);
            else
                set(h.mh.rb.axialLB(1).tanh,'Value',1);
            end
            set(h.mh.et.axialLB.spacing(1),'enable','on')
        end
        
    case h.mh.et.axialLB.stretchingC(3)
        if h.mesh.z.f(2) == 1
            set(h.mh.rb.axialLB(3).lin,'Value',1);
            set(h.mh.et.axialLB.spacing(3),'enable','off')
        else
            if axial_gp > 0
                set(h.mh.rb.axialLB(3).gp,'Value',1);
            else
                set(h.mh.rb.axialLB(3).tanh,'Value',1);
            end
            set(h.mh.et.axialLB.spacing(3),'enable','on')
        end
    
end

guidata(hObject,h)

function mh_bg_radial(hObject, ~, ~)
h = guidata(hObject);

% - since there is only one nonlinear distribution function for each ...
%   direction in each phase, the other button group must change accordingly
% - change value of stretching factor if changing from linear to nonlinear
default = 1.15;
msg_txt = 'Corresponding Stretching Factor f set to default value!';
selection = get(get(hObject,'SelectedObject'),'String');

% general for any button group in radial direction
switch selection
    case 'Hyperbolic Tangent'
        if get(h.mh.rb.radialLB(1).lin,'Value') == 0
            set(h.mh.rb.radialLB(1).tanh,'Value',1);
        end
        if get(h.mh.rb.radialLB(3).lin,'Value') == 0
            set(h.mh.rb.radialLB(3).tanh,'Value',1);
        end
        h.mesh.r.sf={'tanh'};
        
    case 'Geometric Progression'
        if get(h.mh.rb.radialLB(1).lin,'Value') == 0
            set(h.mh.rb.radialLB(1).gp,'Value',1);
        end
        if get(h.mh.rb.radialLB(3).lin,'Value') == 0
            set(h.mh.rb.radialLB(3).gp,'Value',1);
        end
        h.mesh.r.sf={'gp'};
        
    case 'Linear'
        if h.show.stretchingF_lin == 1
            waitfor(helpdlg('Corresponding Stretching Factor f set to 1!','Notice'))
            h.show.stretchingF_lin = 0;
        end
end

% for specific groups
switch hObject
    case h.mh.bg.radialLB(1)
        switch selection
            case {'Hyperbolic Tangent','Geometric Progression'}
                if h.mesh.r.f(1) == 1
                    if h.show.stretchingF == 1
                        waitfor(helpdlg(msg_txt,'Notice'))
                        h.show.stretchingF = 0;
                    end
                    h.mesh.r.f(1) = default;
                    set(h.mh.et.radialLB.stretchingC(1),'String',num2str(default))
                    set(h.mh.et.radialLB.spacing(1),'enable','on')
                end
                
            case 'Linear'
                h.mesh.r.f(1) = 1;
                set(h.mh.et.radialLB.stretchingC(1),'String','1')
                set(h.mh.et.radialLB.spacing(1),'enable','off')
        end
        
    case h.mh.bg.radialLB(3)
        switch selection
            case 'Hyperbolic Tangent'
                set(h.mh.rb.radialLB(3).tanh,'Value',1);
                if h.mesh.r.f(2) == 1
                    if h.show.stretchingF == 1
                        waitfor(helpdlg(msg_txt,'Notice'))
                        h.show.stretchingF = 0;
                    end
                    h.mesh.r.f(2) = default;
                    set(h.mh.et.radialLB.stretchingC(3),'String',num2str(default))
                    set(h.mh.et.radialLB.spacing(3),'enable','on')
                end
                
            case 'Geometric Progression'
                set(h.mh.rb.radialLB(3).gp,'Value',1);
                if h.mesh.r.f(2) == 1
                    if h.show.stretchingF == 1
                        waitfor(helpdlg(msg_txt,'Notice'))
                        h.show.stretchingF = 0;
                    end
                    h.mesh.r.f(2) = default;
                    set(h.mh.et.radialLB.stretchingC(3),'String',num2str(default))
                    set(h.mh.et.radialLB.spacing(3),'enable','on')
                end

            case 'Linear'
                h.mesh.r.f(2) = 1;
                set(h.mh.et.radialLB.stretchingC(3),'String','1')
                set(h.mh.rb.radialLB(3).lin,'Value',1);
                set(h.mh.et.radialLB.spacing(3),'enable','off')
        end
end

guidata(hObject,h)

function mh_bg_axial(hObject, ~, ~)
h = guidata(hObject);

% see description in function mh_bg_radial
default = 1.15;
msg_txt = 'Corresponding Stretching Factor f set to default value!';
selection = get(get(hObject,'SelectedObject'),'String');

% general for any button group in axial direction
switch selection
    case 'Hyperbolic Tangent'
        if get(h.mh.rb.axialLB(1).lin,'Value') == 0
            set(h.mh.rb.axialLB(1).tanh,'Value',1);
        end
        if get(h.mh.rb.axialLB(3).lin,'Value') == 0
            set(h.mh.rb.axialLB(3).tanh,'Value',1);
        end
        h.mesh.z.sf={'tanh'};
        
    case 'Geometric Progression'
        if get(h.mh.rb.axialLB(1).lin,'Value') == 0
            set(h.mh.rb.axialLB(1).gp,'Value',1);
        end
        if get(h.mh.rb.axialLB(3).lin,'Value') == 0
            set(h.mh.rb.axialLB(3).gp,'Value',1);
        end
        h.mesh.z.sf={'gp'};
        
    case 'Linear'
        if h.show.stretchingF_lin == 1
            waitfor(helpdlg('Corresponding Stretching Factor f set to 1!','Notice'))
            h.show.stretchingF_lin = 0;
        end
end

% for specific groups
switch hObject      
    case h.mh.bg.axialLB(1)
        switch selection
            case 'Hyperbolic Tangent'
                if h.mesh.z.f(1) == 1
                    if h.show.stretchingF == 1
                        waitfor(helpdlg(msg_txt,'Notice'))
                        h.show.stretchingF = 0;
                    end
                    h.mesh.z.f(1) = default;
                    set(h.mh.et.axialLB.stretchingC(1),'String',num2str(default))
                    set(h.mh.et.axialLB.spacing(1),'enable','on')
                end

            case 'Geometric Progression'
                if h.mesh.z.f(1) == 1
                    if h.show.stretchingF == 1
                        waitfor(helpdlg(msg_txt,'Notice'))
                        h.show.stretchingF = 0;
                    end
                    h.mesh.z.f(1) = default;
                    set(h.mh.et.axialLB.stretchingC(1),'String',num2str(default))
                    set(h.mh.et.axialLB.spacing(1),'enable','on')
                end

            case 'Linear'
                h.mesh.z.f(1) = 1;
                set(h.mh.et.axialLB.stretchingC(1),'String','1')
                set(h.mh.rb.axialLB(1).lin,'Value',1);
                set(h.mh.et.axialLB.spacing(1),'enable','off')
        end
        
    case h.mh.bg.axialLB(3)
        switch selection
            case 'Hyperbolic Tangent'
                if h.mesh.z.f(2) == 1
                    if h.show.stretchingF == 1
                        waitfor(helpdlg(msg_txt,'Notice'))
                        h.show.stretchingF = 0;
                    end
                    h.mesh.z.f(2) = default;
                    set(h.mh.et.axialLB.stretchingC(3),'String',num2str(default))
                    set(h.mh.et.axialLB.spacing(3),'enable','on')
                end

            case 'Geometric Progression'
                set(h.mh.rb.axialLB(3).gp,'Value',1);
                if h.mesh.z.f(2) == 1
                    if h.show.stretchingF == 1
                        waitfor(helpdlg(msg_txt,'Notice'))
                        h.show.stretchingF = 0;
                    end
                    h.mesh.z.f(2) = default;
                    set(h.mh.et.axialLB.stretchingC(3),'String',num2str(default))
                    set(h.mh.et.axialLB.spacing(3),'enable','on')
                end

            case 'Linear'
                h.mesh.z.f(2) = 1;
                set(h.mh.et.axialLB.stretchingC(3),'String','1')
                set(h.mh.rb.axialLB(3).lin,'Value',1);
                set(h.mh.et.axialLB.spacing(3),'enable','off')
        end
end

guidata(hObject,h)

function mh_pbs(hObject, ~, ~)
h = guidata(hObject);

% create waitbar
    wb = waitbar(0,'Please wait...');

% transfer variables from gui to workspace
    updateValues_SFM
    waitbar(1/4)

% run the setup needed for all computations
    try
        evalin('base','caseSetup_SFM_V3d2')
    catch
        close(wb)
        errordlg('Something went wrong! Check your mesh parameters!','Error 148')
        return
    end
    waitbar(2/4)
    convergence = evalin('base','convergence');
    if convergence == 0
        close(wb)
        waitfor(errordlg('Failed to compute the static surface shape. Check your general and/or geometry parameters!','Error 173'))
        guidata(hObject,h)
        evalin('base','clearvars')
        return
    end

% grid generation
    evalin('base','grid_generation')
    waitbar(3/4)

switch hObject
    case h.mh.pb.drawMesh
        evalin('base','draw_grid')
        if exist('wb','var')
            waitbar(4/4)
            close(wb)
        end
        
    case h.mh.pb.checkMesh
        evalin('base','[mesh.quality.nCells, mesh.quality.aspectRatio, mesh.quality.skewness]=gridQuality(blocks);')
        if exist('wb','var')
            waitbar(4/4)
            close(wb)
        end
        h.nCells      = evalin('base','mesh.quality.nCells');
        h.aspectRatio = evalin('base','mesh.quality.aspectRatio');
        h.skewness    = evalin('base','mesh.quality.skewness');

        msgbox([strcat("Number of Cells = ",num2str(h.nCells));"";strcat("Maximal Aspect Ratio = ",num2str(h.aspectRatio));"";strcat("Maximal Skewness = ",num2str(h.skewness))],'Mesh Info',"help")
    
end

evalin('base','clear convergence')

guidata(hObject,h)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%            ____                    _____        ____                    %
%     ___   |    |   |   |     /\      |     |   |    |   |\  |    ___    %
%    |___   |    |   |   |    /__\     |     |   |    |   | \ |   |___    %
%    |___   |___\|   |___|   /    \    |     |   |____|   |  \|    ___|   %
%                      \                                                  %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function es_cb_energy(hObject, ~, ~)
h = guidata(hObject);

h.flowopt.energy = get(hObject,'Value');

% set to lowest NS model because energy eq. either got switched off or switched on from being switched off
    set(h.gl.rb.NS1,'value',1)
    h.NS = 1; h.flowopt.boussinesq = 1;

if h.flowopt.energy == 1
    set(hObject,'string','On')

    % Enable NS models
        set([h.gl.rb.NS1 h.gl.rb.NS2 h.gl.rb.NS3],'Enable','On')
        % if user selects a customized fluid --> switch off NS3
        if isfile([tempdir 'my_fluids.mat'])
            myfluids = who('-file',[tempdir 'my_fluids.mat']);
            if max(ismember(myfluids,h.selectedLiquid))
                set(h.gl.rb.NS3,'enable','off')
            end
        end
            
    % set visibility on of everything that is connected to temperature
        % Liquid Bridge - radial Direction
        if (h.r_c>0 || h.flowopt.ax==0) && strcmp(h.b1.bc.z(1,3),'c')
            set([h.bc.et.radialLB.rc.temperature h.bc.pm.radialLB.rc.temperature h.bc.st.radialLB.rc.temperature],'visible','on')
        end
        % ---------------------------------------------------------- %
        if strcmp(h.b1.bc.z(2,3),'c')
            set(h.bc.et.radialLB.ri.Bi,'visible','on')
        end
        % Liquid Bridge - axial Direction
        if strcmp(h.b1.bc.r(1,3),'c')
            set([h.bc.et.axialLB.d1.temperature h.bc.pm.axialLB.d1.temperature h.bc.st.axialLB.d1.temperature],'visible','on')
        end
        % ---------------------------------------------------------- %
        if strcmp(h.b1.bc.r(2,3),'c')
            set([h.bc.et.axialLB.d2.temperature h.bc.pm.axialLB.d2.temperature h.bc.st.axialLB.d2.temperature],'visible','on')
        end
            
    % change radio buttons that are connected to temperature
        % Liquid Bridge - radial Direction
        if h.r_c>0 || h.flowopt.ax==0
            if strcmp(h.b1.bc.z(1,3),'c')
                set(h.bc.rb.radialLB.rc.adiabatic,'Value',0,'Enable','on')
                set(h.bc.rb.radialLB.rc.conductive,'Value',1,'Enable','on')
            else
                set(h.bc.rb.radialLB.rc.adiabatic,'Value',1,'Enable','on')
                set(h.bc.rb.radialLB.rc.conductive,'Value',0,'Enable','on')
            end
        end
        % ---------------------------------------------------------- %
        set([h.bc.rb.radialLB.ri.adiabatic h.bc.rb.radialLB.ri.conductive h.bc.rb.axialLB.d1.adiabatic h.bc.rb.axialLB.d1.conductive h.bc.rb.axialLB.d2.adiabatic h.bc.rb.axialLB.d2.conductive],'Enable','on')
        if strcmp(h.b1.bc.z(2,3),'c')
            set(h.bc.rb.radialLB.ri.conductive,'Value',1)
        else
            set(h.bc.rb.radialLB.ri.adiabatic,'Value',1)
        end
        % Liquid Bridge - axial Direction
        if strcmp(h.b1.bc.r(1,3),'c')
            set(h.bc.rb.axialLB.d1.conductive,'Value',1)
        else
            set(h.bc.rb.axialLB.d1.adiabatic,'Value',1)
        end
        % ---------------------------------------------------------- %
        if strcmp(h.b1.bc.r(2,3),'c')
            set(h.bc.rb.axialLB.d2.conductive,'Value',1)
        else
            set(h.bc.rb.axialLB.d2.adiabatic,'Value',1)
        end

    % critical mode on
        set([h.sy.pb.run(2) h.sy.rb.delta_T h.sy.rb.T_d1 h.sy.rb.T_d2],'Enable','On')

else
    set(hObject,'string','Off')

    % Disable NS models
        set([h.gl.rb.NS1 h.gl.rb.NS2 h.gl.rb.NS3],'Enable','Off')
            
    % set visibility off of everything that is connected to temperature
        set([h.bc.et.radialLB.rc.temperature h.bc.pm.radialLB.rc.temperature h.bc.st.radialLB.rc.temperature h.bc.et.radialLB.ri.Bi h.bc.et.axialLB.d1.temperature h.bc.pm.axialLB.d1.temperature h.bc.st.axialLB.d1.temperature h.bc.et.axialLB.d2.temperature h.bc.pm.axialLB.d2.temperature h.bc.st.axialLB.d2.temperature],'visible','off')
        cla(h.bc.as.radialLB.rc(2),'reset'); set(h.bc.as.radialLB.rc(2),'visible','off')
        cla(h.bc.as.radialLB.Bi,'reset'); set(h.bc.as.radialLB.Bi,'visible','off')
        cla(h.bc.as.axialLB.d1(2),'reset'); set(h.bc.as.axialLB.d1(2),'visible','off')
        cla(h.bc.as.axialLB.d2(2),'reset'); set(h.bc.as.axialLB.d2(2),'visible','off')
            
    % change radio buttons that are connected to temperature
        set([h.bc.rb.radialLB.rc.adiabatic h.bc.rb.radialLB.rc.conductive h.bc.rb.radialLB.ri.adiabatic h.bc.rb.axialLB.d1.adiabatic h.bc.rb.axialLB.d1.conductive h.bc.rb.radialLB.ri.conductive h.bc.rb.axialLB.d2.adiabatic h.bc.rb.axialLB.d2.conductive],'Value',0,'Enable','off')

    % critical mode off
        set([h.sy.pb.run(2) h.sy.rb.delta_T h.sy.rb.T_d1 h.sy.rb.T_d2],'Enable','Off')
end

set(h.es.cb.marangoni,'Value',h.flowopt.energy,'Enable',get(hObject,'String'))

guidata(hObject,h)

updateText_SFM_V3d2
updateSketch_SFM

function es_cbs(hObject, ~, ~)
h = guidata(hObject);

switch hObject
    case h.es.cb.stokes
        h.flowopt.creeping = get(hObject,'Value');
        updateText_SFM_V3d2
        
    case h.es.cb.marangoni
        h.flowopt.thermcapcon = get(hObject,'Value');

end

guidata(hObject,h)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%                  ____           _____                                   %
%                 |    |         |                                        %
%                 |  __|         |               ___                      %
%                 |    |    _    |         _    |___     _                %
%                 |____|   |_|   |_____   |_|    ___|   |_|               %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bc_pm_headline(hObject, ~, ~)
h = guidata(hObject);

% here use value in switch function because strings are changing in terms
% of planar or axisymmetric
val=get(hObject, 'Value');

set([h.bc.pl.radialLB.main h.bc.pl.axialLB.main],'visible','off')

switch val
    case 1
        set(h.bc.pl.radialLB.main,'visible','on')
        h.bc.Visible = 'radial';
    case 2
        set(h.bc.pl.axialLB.main,'visible','on')
        h.bc.Visible = 'axial';
end

guidata(hObject,h)

function bc_bg_velocity(hObject, ~, ~)
h = guidata(hObject);

selection = get(get(hObject,'SelectedObject'),'String');

switch selection
    case 'No Slip'
        switch hObject
            case h.bc.bg.axialLB.d1(2)
                h.b1.bc.r(1,2) = 'n';

            case h.bc.bg.axialLB.d2(2)
                h.b1.bc.r(2,2) = 'n';

            case h.bc.bg.radialLB.rc(2)
                h.b1.bc.z(1,2) = 'n';
        end
        
    case 'Slip'
        switch hObject
            case h.bc.bg.axialLB.d1(2)
                h.b1.bc.r(1,2) = 's';
                
            case h.bc.bg.axialLB.d2(2)
                h.b1.bc.r(2,2) = 's';
                
            case h.bc.bg.radialLB.rc(2)
                h.b1.bc.z(1,2) = 's';
        end
end

guidata(hObject,h)

function bc_bg_temperature(hObject, ~, ~)
h = guidata(hObject);

p = h.flowopt.ax + 1; rx = 'xr'; zy = 'yz';
selection = get(get(hObject,'SelectedObject'),'String');

switch selection
    case 'Adiabatic'
        switch hObject
            case h.bc.bg.axialLB.d1(3)
                h.b1.bc.r(1,3) = 'a'; h.b1.bc.rhs.T.r{1} = '@(r)0';
                cla(h.bc.as.axialLB.d1(2),'reset'); set(h.bc.as.axialLB.d1(2),'visible','off')
                set([h.bc.et.axialLB.d1.temperature h.bc.pm.axialLB.d1.temperature h.bc.st.axialLB.d1.temperature],'visible','off')
                    
            case h.bc.bg.axialLB.d2(3)
                h.b1.bc.r(2,3) = 'a'; h.b1.bc.rhs.T.r{2} = '@(r)0';
                cla(h.bc.as.axialLB.d2(2),'reset'); set(h.bc.as.axialLB.d2(2),'visible','off')
                set([h.bc.et.axialLB.d2.temperature h.bc.pm.axialLB.d2.temperature h.bc.st.axialLB.d2.temperature],'visible','off')
                    
            case h.bc.bg.radialLB.rc(3)
                h.b1.bc.z(1,3) = 'a'; h.b1.bc.rhs.T.z{1} = '@(z)0';
                cla(h.bc.as.radialLB.rc(2),'reset'); set(h.bc.as.radialLB.rc(2),'visible','off')
                set([h.bc.et.radialLB.rc.temperature h.bc.pm.radialLB.rc.temperature h.bc.st.radialLB.rc.temperature],'visible','off')
                
            case h.bc.bg.radialLB.ri(2)
                h.b1.bc.z(2,3) = 'a';
                set(h.bc.et.radialLB.ri.Bi,'visible','off')
                cla(h.bc.as.radialLB.Bi)
        end
        
    case 'Conductive'
        switch hObject
            case h.bc.bg.axialLB.d1(3)
                h.b1.bc.r(1,3) = 'c';
                if isempty(get(h.bc.et.axialLB.d1.temperature,'String'))
                    h.b1.bc.rhs.T.r{1} = [];
                else
                    h.b1.bc.rhs.T.r{1} = get(h.bc.et.axialLB.d1.temperature,'string');
                    h.b1.bc.rhs.T.r{1} = convertTemperature(h,h.b1.bc.rhs.T.r{1},'r');
                end

                set([h.bc.et.axialLB.d1.temperature h.bc.pm.axialLB.d1.temperature h.bc.st.axialLB.d1.temperature],'visible','on')
                axes(h.bc.as.axialLB.d1(2));
                text(0.94*h.tp_r,0.47,'$^\circ$C',h.latex_r{:})
                if get(h.bc.pm.axialLB.d1.temperature,'value') == 2
                    text(0.4*h.tp_l,0.47,['$T_{d1}(' rx(p) ') =$'],h.latex_l{:})
                else
                    text(0.4*h.tp_l,0.47,'$T_{d1} =$',h.latex_l{:})
                end
                displayTemperature = get(h.bc.et.axialLB.d1.temperature,'String');
                idx = strfind(displayTemperature,rx(3-p));
                if ~isempty(idx)
                    displayTemperature(idx) = rx(p);
                    set(h.bc.et.axialLB.d1.temperature,'String',displayTemperature);
                end
                
            case h.bc.bg.axialLB.d2(3)
                h.b1.bc.r(2,3) = 'c';
                if isempty(get(h.bc.et.axialLB.d1.temperature,'String'))
                    h.b1.bc.rhs.T.r{2} = [];
                else
                    h.b1.bc.rhs.T.r{2} = get(h.bc.et.axialLB.d2.temperature,'string');
                    h.b1.bc.rhs.T.r{2} = convertTemperature(h,h.b1.bc.rhs.T.r{2},'r');
                end

                set([h.bc.et.axialLB.d2.temperature h.bc.pm.axialLB.d2.temperature h.bc.st.axialLB.d2.temperature],'visible','on')
                axes(h.bc.as.axialLB.d2(2));
                text(0.94*h.tp_r,0.47,'$^\circ$C',h.latex_r{:})
                if get(h.bc.pm.axialLB.d2.temperature,'value') == 2
                    text(0.4*h.tp_l,0.47,['$T_{d2}(' rx(p) ') =$'],h.latex_l{:})
                else
                    text(0.4*h.tp_l,0.47,'$T_{d2} =$',h.latex_l{:})
                end
                displayTemperature = get(h.bc.et.axialLB.d2.temperature,'String');
                idx = strfind(displayTemperature,rx(3-p));
                if ~isempty(idx)
                    displayTemperature(idx) = rx(p);
                    set(h.bc.et.axialLB.d2.temperature,'String',displayTemperature);
                end
                
            case h.bc.bg.radialLB.rc(3)
                h.b1.bc.z(1,3) = 'c';
                if isempty(get(h.bc.et.radialLB.rc.temperature,'String'))
                    h.b1.bc.rhs.T.z{1} = [];
                else
                    h.b1.bc.rhs.T.z{1} = get(h.bc.et.radialLB.rc.temperature,'string');
                    h.b1.bc.rhs.T.z{1} = convertTemperature(h,h.b1.bc.rhs.T.z{1},'z');
                end

                set([h.bc.et.radialLB.rc.temperature h.bc.pm.radialLB.rc.temperature h.bc.st.radialLB.rc.temperature],'visible','on')
                axes(h.bc.as.radialLB.rc(2))
                text(0.94*h.tp_r,0.47,'$^\circ$C',h.latex_r{:})
                if get(h.bc.pm.radialLB.rc.temperature,'value') == 2
                    text(0.4*h.tp_l,0.47,['$T_{' rx(p) '_c}(' zy(p) ') =$'],h.latex_l{:})
                else
                    text(0.4*h.tp_l,0.47,['$T_{' rx(p) '_c} =$'],h.latex_l{:})
                end
                displayTemperature = get(h.bc.et.radialLB.rc.temperature,'String');
                idx = strfind(displayTemperature,zy(3-p));
                if ~isempty(idx)
                    displayTemperature(idx) = zy(p);
                    set(h.bc.et.radialLB.rc.temperature,'String',displayTemperature);
                end
                
            case h.bc.bg.radialLB.ri(2)
                h.b1.bc.z(2,3) = 'c';
        end
        
    case 'Newton''s law'
        if h.Bi == 0, h.Bi = 0.2; end
        set(h.bc.et.radialLB.ri.Bi,'string',num2str(h.Bi),'visible','on')
        h.b1.bc.z(2,3) = 'c';
        axes(h.bc.as.radialLB.Bi); cla(h.bc.as.radialLB.Bi)
        text(h.tp_l,h.tp_y,'$\mathrm{Bi}=$',h.latex_l{:})
end

guidata(hObject,h)

function bc_bg_axis(hObject, ~, ~)
h = guidata(hObject);

selection = get(get(hObject,'SelectedObject'),'String');
switch selection
    case 'Symmetry'
        if h.r_c > 0
            rx = 'xr'; p = h.flowopt.ax + 1;
            waitfor(warndlg(['Please, first set ' rx(p) '_c to 0 before using this boundary condition!'],'Warning',h.opts))
            set(h.bc.rb.radialLB.rc.wall,'value',1)
            return
        end

        h.b1.bc.z(1,1:3) = 'a  ';
        set([h.bc.rb.radialLB.rc.noSlip h.bc.rb.radialLB.rc.slip h.bc.rb.radialLB.rc.adiabatic h.bc.rb.radialLB.rc.conductive],'Value',0,'Enable','off')
    
         % cancel visibility of imposed temperature profile at r=0
        cla(h.bc.as.radialLB.rc(2),'reset');
        set([h.bc.as.radialLB.rc(2) h.bc.et.radialLB.rc.temperature h.bc.pm.radialLB.rc.temperature h.bc.st.radialLB.rc.temperature],'visible','off')
        
    case 'No Penetration'
        if h.flowopt.ax == 1
            waitfor(warndlg('Please, first select a non-zero value for r_c before using this boundary condition.','Warning',h.opts))
            set(h.bc.rb.radialLB.rc.axis,'value',1)
            return
        end

        h.b1.bc.z(1,1:2) = 'wn';
        set(h.bc.rb.radialLB.rc.slip,'Value',0,'Enable','on')
        set(h.bc.rb.radialLB.rc.noSlip,'Value',1,'Enable','on')

        if h.flowopt.energy == 1
            h.b1.bc.z(1,3) = 'c';
            set(h.bc.rb.radialLB.rc.adiabatic,'Value',0,'Enable','on')
            set(h.bc.rb.radialLB.rc.conductive,'Value',1,'Enable','on')
            set([h.bc.et.radialLB.rc.temperature h.bc.pm.radialLB.rc.temperature h.bc.st.radialLB.rc.temperature],'visible','on')
        else
            h.b1.bc.z(1,3) = 'a';
            set([h.bc.rb.radialLB.rc.adiabatic h.bc.rb.radialLB.rc.conductive],'Value',0,'Enable','off')
        end

end

updateSketch_SFM
updateText_SFM_V3d2

guidata(hObject,h)

function bc_bg_surface(hObject, ~, ~)
h = guidata(hObject);

selection = get(get(hObject,'SelectedObject'),'String');
switch selection
    case 'Dynamically deformed surface shape'
        h.b1.bc.z(2,2) = 's'; set(h.gl.et.Vr,'Enable','on')
        
    case 'Indeformable hydrostatic surface shape'
        h.b1.bc.z(2,2) = 'i'; set(h.gl.et.Vr,'Enable','on')
        
    case 'Straight indeformable surface shape'
        h.b1.bc.z(2,2) = 'r'; set(h.gl.et.Vr,'String','1','Enable','off');
        if h.show.surface == 1 && h.V_r ~= 1
            h.V_r = 1;
            waitfor(helpdlg('Volume ratio set to 1.','Attention'))
            h.show.surface = 0;
        end           
end

updateSketch_SFM

guidata(hObject,h)

function bc_et_temperature(hObject, ~, ~)
h = guidata(hObject);

p = h.flowopt.ax + 1; zy = 'yz';

eval([h.vlNm '=''@(' h.vlNm(end-3) ')' num2str(eval(h.vlNm),16) ''';'])

switch h.vlNm
    case {'h.b1.bc.rhs.T.r{1}','h.b1.bc.rhs.T.r{2}'}
        if strcmp(h.vlNm,'h.b1.bc.rhs.T.r{1}')
            h.temperature.ri_d1 = str2double(h.b1.bc.rhs.T.r{1}(5:end)); h.temperature.rc_d1 = h.temperature.ri_d1; h.T_d1l = h.temperature.ri_d1;
        else
            h.temperature.ri_d2 = str2double(h.b1.bc.rhs.T.r{2}(5:end)); h.temperature.rc_d2 = h.temperature.ri_d2; h.T_d2l = h.temperature.ri_d2;
        end
        
        updateSketch_SFM

        % change temperature profile at central rod in order to get a linear profile with continous temperatures in the corners
        if get(h.bc.pm.radialLB.rc.temperature,'value') == 2
            if h.temperature.rc_d1 >= h.temperature.rc_d2
                h.b1.bc.rhs.T.z{1} = [num2str((h.temperature.rc_d1+h.temperature.rc_d2)*0.5) '+' num2str(h.temperature.rc_d1-h.temperature.rc_d2) '*' zy(p) '/d'];
            else
                h.b1.bc.rhs.T.z{1} = [num2str((h.temperature.rc_d1+h.temperature.rc_d2)*0.5) '-' num2str(h.temperature.rc_d2-h.temperature.rc_d1) '*' zy(p) '/d'];
            end
            set(h.bc.et.radialLB.rc.temperature,'string',h.b1.bc.rhs.T.z{1})
            h.b1.bc.rhs.T.z{1} = convertTemperature(h,h.b1.bc.rhs.T.z{1},'z');
        end
        
    case 'h.b1.bc.rhs.T.z{2}'
        h.Bi = str2double(h.b1.bc.rhs.T.z{2}(5:end));
        
end

guidata(hObject,h)

function bc_pm_temperature(hObject, ~, ~)
h = guidata(hObject);

p = h.flowopt.ax + 1; rx = 'xr'; zy = 'yz'; ax = h.flowopt.ax;
val=get(hObject, 'Value');
if val == 0, val = 2; end

switch hObject
    case {h.bc.pm.axialLB.d1.temperature,h.bc.et.axialLB.d1.temperature}
        axes(h.bc.as.axialLB.d1(2)); cla(h.bc.as.axialLB.d1(2));
        text(0.94*h.tp_r,0.47,'$^\circ$C',h.latex_r{:})
        switch val
            case 1
                text(0.4*h.tp_l,0.47,'$T_{d1} =$',h.latex_l{:})
                set(h.bc.et.axialLB.d1.temperature,'Enable','on')
                
                % if users selects constant --> program checks whether edit box is a number
                % if yes --> leave it; else delete the string and set temperature empty
                if ~isempty(get(h.bc.et.axialLB.d1.temperature,'String'))
                    [~, status]=str2num(get(h.bc.et.axialLB.d1.temperature,'string'));
                    if status == 0
                        set(h.bc.et.axialLB.d1.temperature,'String','')
                        h.b1.bc.rhs.T.r{1} = [];
                    end
                end
            case 2
                text(0.4*h.tp_l,0.47,['$T_{d1}(' rx(p) ') =$'],h.latex_l{:})
                set(h.bc.et.axialLB.d1.temperature,'Enable','inactive')
                if contains(get(h.bc.et.axialLB.d1.temperature,'string'),rx(p)) == 0
                    set(h.bc.et.axialLB.d1.temperature,'string','')
                end
                
                h.b1.bc.rhs.T.r{1} = inputTemperature_SFM(rx(p),h.r_c,h.flowopt.ax);
                if isempty(h.b1.bc.rhs.T.r{1}) && contains(get(h.bc.et.axialLB.d1.temperature,'string'),rx(p))
                    h.b1.bc.rhs.T.r{1} = get(h.bc.et.axialLB.d1.temperature,'string');
                end
                set(h.bc.et.axialLB.d1.temperature,'String',h.b1.bc.rhs.T.r{1})
                h.b1.bc.rhs.T.r{1} = convertTemperature(h,h.b1.bc.rhs.T.r{1},'r');
        end
        % change temperature profile at central rod in order to get a linear profile with continous temperatures in the corners
        if ~isempty(h.b1.bc.rhs.T.r{1})
            tempProf = str2func(h.b1.bc.rhs.T.r{1});
            h.temperature.rc_d1 = tempProf(h.r_c);
            if get(h.bc.pm.radialLB.rc.temperature,'value') == 2
                if h.temperature.rc_d1 >= h.temperature.rc_d2
                    h.b1.bc.rhs.T.z{1} = [num2str((h.temperature.rc_d1+h.temperature.rc_d2)*0.5) '+' num2str(h.temperature.rc_d1-h.temperature.rc_d2) '*' zy(p) '/d'];
                else
                    h.b1.bc.rhs.T.z{1} = [num2str((h.temperature.rc_d1+h.temperature.rc_d2)*0.5) '-' num2str(h.temperature.rc_d2-h.temperature.rc_d1) '*' zy(p) '/d'];
                end
                set(h.bc.et.radialLB.rc.temperature,'string',h.b1.bc.rhs.T.z{1})
                h.b1.bc.rhs.T.z{1} = convertTemperature(h,h.b1.bc.rhs.T.z{1},'z');
            end
            h.T_d1l = integral(@(r)tempProf(r).*r.^ax,h.r_c,h.r_i)*2/(h.r_i^2-h.r_c^2);
            h.temperature.ri_d1 = tempProf(h.r_i);
        end
        
    case {h.bc.pm.axialLB.d2.temperature,h.bc.et.axialLB.d2.temperature}
        axes(h.bc.as.axialLB.d2(2)); cla(h.bc.as.axialLB.d2(2));
        text(0.94*h.tp_r,0.47,'$^\circ$C',h.latex_r{:})
        switch val
            case 1
                text(0.4*h.tp_l,0.47,'$T_{d1} =$',h.latex_l{:})
                set(h.bc.et.axialLB.d2.temperature,'Enable','on')

                if ~isempty(get(h.bc.et.axialLB.d2.temperature,'String'))
                    [~, status]=str2num(get(h.bc.et.axialLB.d2.temperature,'string'));
                    if status == 0
                        set(h.bc.et.axialLB.d2.temperature,'String','')
                        h.b1.bc.rhs.T.r{2} = [];
                    end
                end
            case 2
                text(0.4*h.tp_l,0.47,['$T_{d2}(' rx(p) ') =$'],h.latex_l{:})
                set(h.bc.et.axialLB.d2.temperature,'Enable','inactive')
                if contains(get(h.bc.et.axialLB.d2.temperature,'string'),rx(p)) == 0
                    set(h.bc.et.axialLB.d2.temperature,'string','')
                end
                
                h.b1.bc.rhs.T.r{2} = inputTemperature_SFM(rx(p),h.r_c,h.flowopt.ax);
                if isempty(h.b1.bc.rhs.T.r{2}) && contains(get(h.bc.et.axialLB.d2.temperature,'string'),rx(p))
                    h.b1.bc.rhs.T.r{2} = get(h.bc.et.axialLB.d2.temperature,'string');
                end
                set(h.bc.et.axialLB.d2.temperature,'String',h.b1.bc.rhs.T.r{2})
                h.b1.bc.rhs.T.r{2} = convertTemperature(h,h.b1.bc.rhs.T.r{2},'r');
        end
        % change temperature profile at central rod in order to get a linear profile with continous temperatures in the corners
        if ~isempty(h.b1.bc.rhs.T.r{2})
            tempProf = str2func(h.b1.bc.rhs.T.r{2});
            h.temperature.rc_d2 = tempProf(h.r_c);
            if get(h.bc.pm.radialLB.rc.temperature,'value') == 2
                if h.temperature.rc_d1 >= h.temperature.rc_d2
                    h.b1.bc.rhs.T.z{1} = [num2str((h.temperature.rc_d1+h.temperature.rc_d2)*0.5) '+' num2str(h.temperature.rc_d1-h.temperature.rc_d2) '*' zy(p) '/d'];
                else
                    h.b1.bc.rhs.T.z{1} = [num2str((h.temperature.rc_d1+h.temperature.rc_d2)*0.5) '-' num2str(h.temperature.rc_d2-h.temperature.rc_d1) '*' zy(p) '/d'];
                end
                set(h.bc.et.radialLB.rc.temperature,'string',h.b1.bc.rhs.T.z{1})
                h.b1.bc.rhs.T.z{1} = convertTemperature(h,h.b1.bc.rhs.T.z{1},'z');
            end
            h.T_d2l = integral(@(r)tempProf(r).*r.^ax,h.r_c,h.r_i)*2/(h.r_i^2-h.r_c^2);
            h.temperature.ri_d2 = tempProf(h.r_i);
        end
        
    case {h.bc.pm.radialLB.rc.temperature,h.bc.et.radialLB.rc.temperature}
        axes(h.bc.as.radialLB.rc(2)); cla(h.bc.as.radialLB.rc(2));
        text(0.94*h.tp_r,0.47,'$^\circ$C',h.latex_r{:})
        switch val
            case 1
                text(0.4*h.tp_l,0.47,['$T_{' rx(p) '_c} =$'],h.latex_l{:})
                set(h.bc.et.radialLB.rc.temperature,'Enable','on')
                
                if ~isempty(get(h.bc.et.radialLB.rc.temperature,'String'))
                    [~, status]=str2num(get(h.bc.et.radialLB.rc.temperature,'string'));
                    if status == 0
                        set(h.bc.et.radialLB.rc.temperature,'String','')
                        h.b1.bc.rhs.T.z{1} = [];
                    end
                end
            case 2
                text(0.4*h.tp_l,0.47,['$T_{' rx(p) '_c}(' zy(p) ') =$'],h.latex_l{:})
                set(h.bc.et.radialLB.rc.temperature,'Enable','inactive')
                if contains(get(h.bc.et.radialLB.rc.temperature,'string'),zy(p)) == 0
                    set(h.bc.et.radialLB.rc.temperature,'string','')
                end
                
                h.b1.bc.rhs.T.z{1} = inputTemperature_SFM(zy(p),h.r_c,h.flowopt.ax);
                if isempty(h.b1.bc.rhs.T.z{1}) && contains(get(h.bc.et.radialLB.rc.temperature,'string'),zy(p))
                    h.b1.bc.rhs.T.z{1} = get(h.bc.et.radialLB.rc.temperature,'string');
                end
                set(h.bc.et.radialLB.rc.temperature,'String',h.b1.bc.rhs.T.z{1})
                h.b1.bc.rhs.T.z{1} = convertTemperature(h,h.b1.bc.rhs.T.z{1},'z');
        end
end

updateSketch_SFM

guidata(hObject,h)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%                                                           ____          %
%  __    |    |\  /|   |    |   |        /\    _____   |   |    |   |\  | %
% |__    |    | \/ |   |    |   |       /__\     |     |   |    |   | \ | %
%  __|   |    |    |   |____|   |___   /    \    |     |   |____|   |  \| %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sn_et_steps(hObject, ~, ~)
h = guidata(hObject);

if floor(h.steps) ~= h.steps || h.steps <= 0
    if h.show.steps.warning == 1
        waitfor(warndlg('Value must be a positive integer.','Warning'))
    end
    h.show.steps.warning = 0;
    h.steps = floor(abs(h.steps));
    set(h.sn.et.steps,'string',num2str(h.steps))
end

guidata(hObject,h)

function sn_bg_initialization(hObject, ~, ~)
h = guidata(hObject);

selection = get(get(hObject,'SelectedObject'),'String');

switch selection
    case 'Fluid at rest. Constant mean temperature.'
        h.initialize = 'standard';
    case 'Use previous simulation.'
        if evalin('base','exist(''b1'',''var'') && isfield(b1,''w'') && residual<1e15') == 0
            waitfor(errordlg('Basic state solution not found!','Error 404'))
            set(h.sn.rb.standard,'value',1)
            return
        end
        h.initialize = 'previous';
end

guidata(hObject,h)

function sn_pb_run(hObject, ~, ~)
h = guidata(hObject);

% checking if everything is well defined
    set(h.sn.st.run,'visible','off')

    wb = waitbar(0,'Checking case ...');
    waitbar(1/4)

    % transfer variables from gui to workspace
    if strcmp(h.initialize,'standard')
        evalin('base','clearvars; clc;')
    end

    updateValues_SFM
    waitbar(2/4)

    % run the setup needed for all computations
    try
        evalin('base','caseSetup_SFM_V3d2')
    catch
        close(wb)
        errordlg('Something went wrong! Check your mesh parameters!','Error 148')
        return
    end
    waitbar(3/4)

    % delete perturbation flow if it exists
    if isfield(h,'old'), h=rmfield(h,'old'); end
    evalin('base','if isfield(b1,''w_hat''), b1=rmfield(b1,''w_hat''), end;')

    convergence = evalin('base','convergence');
    if convergence == 0
        close(wb)
        waitfor(errordlg('Failed to compute the surface shape. Check your general and/or geometry parameters!','Error 173'))
        return
    end

    waitbar(3/3)
    close(wb)

% changing button options
    set(h.sn.st.run,'String','running...')
    set([h.sn.st.run h.sn.pb.stop],'visible','on')
    set(h.sn.pb.run,'Visible','off')
    axes(h.sy.as.eigs); cla(h.sy.as.eigs,'reset'); set(h.sy.as.eigs,'visible','off')
    axes(h.sy.as.dependent); cla(h.sy.as.dependent,'reset'); set(h.sy.as.dependent,'visible','off')
    set([h.tabs{:} h.sn.et.residuals(1) h.sn.et.residuals(2) h.sn.rb.standard h.sn.rb.previous h.sn.et.steps h.sn.pb.load h.sn.pb.save h.sn.pb.plot],'Enable','inactive')

    axes(h.sn.as.residuals)
    if strcmp(h.initialize,'standard')
        cla(h.sn.as.residuals,'reset');
        ax = gca; ax.YScale = 'log'; ax.XLim = [1 2]; ax.XTick = [1 2];
    end
    ax = gca; ax.Box = 'on'; xlabel('Iterations'); ylabel('Residuals');
    
    
% run simulation
    % transfer variables from gui to workspace
    if strcmp(h.initialize,'standard')
        evalin('base','clearvars; clc;')
    end
    updateValues_SFM
    evalin('base','caseSetup_SFM_V3d2')
    try
        evalin('base','flowopt.stability=0; grid_generation; block_assembly_V3d2');
    catch
        errordlg('Initialization failed.','Error 158')
        set([h.sn.st.run h.sn.pb.stop],'visible','off')
        set(h.sn.pb.run,'Visible','on')
        set([h.tabs{:} h.sn.et.residuals(1) h.sn.et.residuals(2) h.sn.rb.standard h.sn.rb.previous h.sn.et.steps h.sn.pb.load h.sn.pb.save h.sn.pb.plot],'Enable','on')
        return
    end

    if strcmp(h.initialize,'standard')
        evalin('base','res.u_vec = []; res.w_vec = []; res.T_vec = []; res.conti_vec = []; res.newton_vec = [];')
    end

    set(h.sn.pb.stop,'userdata',0);
    try
        for iteration = 1:h.steps
            assignin('base','iteration',iteration)
            evalin('base','dx = J\(A*x-b);')
            evalin('base','x = x - dx;')
            evalin('base','block_assembly_V3d2')

            evalin('base','plotResiduals')
            residual = evalin('base', 'residual');
            dx = evalin('base', 'dx');
            if (residual <= h.flowopt.tolerance.residuals && norm(dx) <= h.flowopt.tolerance.newton) || isnan(residual) || residual>10^15 || get(h.sn.pb.stop,'userdata'), break, end
        end
    catch
        errordlg('Simulation stopped!','Error 201')
        set([h.tabs{:} h.sn.et.residuals(1) h.sn.et.residuals(2) h.sn.rb.standard h.sn.rb.previous h.sn.et.steps h.sn.pb.load h.sn.pb.save h.sn.pb.plot],'Enable','on')
        set([h.sn.st.run h.sn.pb.stop],'visible','off')
        set(h.sn.pb.run,'Visible','on')
        return
    end
    evalin('base','clear A J b convergence')
    evalin('base','old.x=x;')
    
% analyze simulation
    h.initialize = 'previous';
    if residual <= h.flowopt.tolerance.residuals && norm(dx) <= h.flowopt.tolerance.newton
        helpdlg('Solution converged!')
    elseif iteration == h.steps
        warndlg('Basic State Solution not yet converged!','Attention')
    elseif residual > 10^15 || isnan(residual)
        errordlg('Simulation stopped!','Error 201')
        h.initialize = 'standard';
    elseif get(h.sn.pb.stop,'userdata')
        helpdlg('Simulation stopped!','Notice')
    end

    if strcmp(h.initialize,'previous')
        set(h.sn.rb.previous,'Value',1)
    else
        set(h.sn.rb.standard,'Value',1)
    end

    set([h.tabs{:} h.sn.et.residuals(1) h.sn.et.residuals(2) h.sn.rb.standard h.sn.rb.previous h.sn.et.steps h.sn.pb.load h.sn.pb.save h.sn.pb.plot],'Enable','on')
    set([h.sn.st.run h.sn.pb.stop],'visible','off')
    set(h.sn.pb.run,'Visible','on')

guidata(hObject,h)

function sn_pbs(hObject, ~, ~)
h = guidata(hObject);

switch hObject
    case h.sn.pb.stop
        set(h.sn.pb.stop,'userdata',1)
        
    case {h.sn.pb.load,h.sy.pb.load}
        set([h.sn.st.run h.sy.st.run],'String','loading...','visible','on')
        try
            [file,path] = uigetfile('*.mat','Open');
            if ischar(file)
                evalin('base','clear all; clc')
                evalin('base',['load(''' [path file] ''')'])
                loadSetup_SFM_V3d2
            end
        catch
            set([h.sn.st.run h.sy.st.run],'visible','off');
            errordlg('Invalid setup file.','Error 315')
        end
        set([h.sn.st.run h.sy.st.run],'visible','off');
        
    case {h.sn.pb.save,h.sy.pb.save}
        set([h.sn.st.run h.sy.st.run],'String','saving...','visible','on')
        
        % save setup if there is no basic state
        flg = evalin('base','exist(''b1'',''var'') && isfield(b1,''w'')');
        if flg == 0
            updateValues_SFM
            evalin('base','caseSetup_SFM_V3d2')
        end

        evalin('base','[file,path] = uiputfile(''newData.mat'');')
        evalin('base','if max(file~=0) || max(path~=0), save(fullfile(path, file)), end')

        set([h.sn.st.run h.sy.st.run],'visible','off')
        
    case h.sn.pb.plot
        flg = evalin('base','exist(''b1'',''var'') && isfield(b1,''w'') && residual<1e15');
        if flg == 0
            waitfor(errordlg('Basic state solution not found!','Error 404'))
            return
        else
            plot_bS
        end
        
    case h.sy.pb.plot
        flg = evalin('base','exist(''b1'',''var'') && isfield(b1,''w_hat'') && residual<1e15');
        if flg == 0
            waitfor(errordlg('Perturbation flow not found!','Error 404'))
            return
        else
            plot_pF
        end
        
end

guidata(hObject,h)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%                                                                         %
%                            optical ray                                  %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function or_et_ray(hObject, ~, ~)
h = guidata(hObject);

r_0 = h.r_0*h.r_i/h.l_lb;

assignin('base','r_0',r_0);

guidata(hObject,h)

function or_et_particle(hObject, ~, ~)
h = guidata(hObject);

z_p = h.z_p+0.5;

assignin('base','z_p',z_p);

guidata(hObject,h)

function or_et_N(hObject, ~, ~)
h = guidata(hObject);

assignin('base','N_coeff',[h.N_coeff(1) h.N_coeff(2) h.N_coeff(3)]);

guidata(hObject,h)

function or_pb_run(hObject, ~, ~)
h = guidata(hObject);

p = h.flowopt.ax + 1; % pointer
rx = 'xr'; zy = 'yz';

flg = evalin('base','exist(''b1'',''var'') && isfield(b1,''w'') && residual<1e15');
if flg == 0
    waitfor(errordlg('Basic state solution not found!','Error 404'))
    return
end

T = evalin('base','b1.T.T');
R = evalin('base','b1.R.T');
Z = evalin('base','b1.Z.T');
Z = h.l_lb/2-Z;

% index of refraction
    N = h.N_coeff(1) - h.N_coeff(2)*(T-h.N_coeff(3));

% compute optical ray
    ray = evalin('base','optical_ray(b1,r_0,z_p,N_coeff)');
    assignin('base','ray',ray);

% plot the index of refraction
    axes(h.or.as.plot(1))
    cla(h.or.as.plot(1),'reset');
    [~, c] = contourf(R*1000,Z*1000,N,256);
    set(c,'Edgecolor','none');
    colormap(h.or.as.plot(1),'pink')
    axis equal
    cb = colorbar; set(get(cb,'ylabel'),'String','Index of refraction $\mathcal{N}(T)$','Interpreter','latex','FontSize',12)
    xlabel(['$' rx(p) '$ [mm]'],'Interpreter','latex','FontSize',12)
    ylabel(['$' zy(p) '$ [mm]'],'Interpreter','latex','FontSize',12)
    hold on
    plot(ray.path_real(:,1)*h.l_lb*1000,(ray.path_real(:,2)-0.5)*h.l_lb*1000,'color',h.selectedTabColor)
    plot(ray.path_fict(:,1)*h.l_lb*1000,(ray.path_fict(:,2)-0.5)*h.l_lb*1000,'LineStyle','--','Color',h.unselectedTabColor)
    
    
% plot the real and the fictional path of the ray
    axes(h.or.as.plot(2))
    cla(h.or.as.plot(2),'reset');
    plot(ray.path_real(:,1)*h.l_lb*1000,(ray.path_real(:,2)-0.5)*h.l_lb*1000,'Color',h.selectedTabColor)
    hold on
    plot(ray.path_fict(:,1)*h.l_lb*1000,(ray.path_fict(:,2)-0.5)*h.l_lb*1000,'LineStyle','--','Color',h.unselectedTabColor)
    xlabel(['$' rx(p) '$ [mm]'],'Interpreter','latex','FontSize',12)
    ylabel(['$' zy(p) '$ [mm]'],'Interpreter','latex','FontSize',12)
    % adjust xlim because otherwise the vertical line is always at the border
    x_lim = xlim;
    delta_r = (ray.path_real(end,1)-ray.r_0)*h.l_lb*1000;
    if delta_r < 0
        xlim([x_lim(1) x_lim(2)+diff(x_lim)*0.05])
    else
        xlim([x_lim(1)-diff(x_lim)*0.05 x_lim(2)])
    end
    % adjust ylim such that upper bound = 0.5*l_lb
    y_lim = ylim;
    ylim([y_lim(1) h.l_lb*1000/2])
    legend({'$\mathcal{N}=\mathcal{N}(T)$','$\mathcal{N}=\mathrm{const}.$'},'location','best','Interpreter','latex','FontSize',12)
    
% print the actual position of the particle
    axes(h.or.as.particle(1)); cla(h.or.as.particle(1))
    if isnan(ray.r_p)
        text(-0.53,h.tp_y,'$r_\mathrm{p}=$',h.latex_r{:})
        helpdlg('The ray reached the free surface.','Info')
    else
        text(-0.53,h.tp_y,['$' rx(p) '_\mathrm{p}=' num2str(ray.r_p*h.l_lb/h.r_i) '\cdot\, ' rx(p) '_i$'],h.latex_r{:})
    end

guidata(hObject,h)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%             _____          ____                _____                    %
%         __    |     /\    |    |  |  |      |    |     \  /             %
%        |__    |    /__\   |  __|  |  |      |    |      \/              %
%         __|   |   /    \  |    |  |  |      |    |      |               %
%               |  /      \ |____|  |  |____  |    |      |               %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sy_et_eigs(hObject, ~, ~)
h = guidata(hObject);

switch hObject
    case h.sy.et.n_eig
        if floor(h.flowopt.eigs.n) ~= h.flowopt.eigs.n || h.flowopt.eigs.n <= 0
            if h.show.n_eig.warning == 1
                waitfor(warndlg('Value must be a positive integer.','Warning'))
            end
            h.show.n_eig.warning = 0;
            h.flowopt.eigs.n = floor(abs(h.flowopt.eigs.n));
            if h.flowopt.eigs.n == 0, h.flowopt.eigs.n = 1; end
            set(h.sy.et.n_eig,'string',num2str(h.flowopt.eigs.n))
        end
    case h.sy.et.n_eig_c
        if floor(h.flowopt.eigs.n_cayley) ~= h.flowopt.eigs.n_cayley || h.flowopt.eigs.n_cayley < 0
            if h.show.n_eig.warning == 1
                waitfor(warndlg('Value must be a positive integer (or 0).','Warning'))
            end
            h.show.n_eig.warning = 0;
            h.flowopt.eigs.n_cayley = floor(abs(h.flowopt.eigs.n_cayley));
            set(h.sy.et.n_eig_c,'string',num2str(h.flowopt.eigs.n_cayley))
        end
        if h.flowopt.eigs.n_cayley > 5 % default value
            if h.show.n_cayley.warning == 1
                waitfor(warndlg({'Increasing n_{cay} beyond its default value may cause a disproportional increase in computational time.','','Please reconsider changing to the default value!'},'Attention',h.opts))
            end
            h.show.n_cayley.warning = 0;
        end
    case h.sy.et.kryl
        if floor(h.flowopt.eigs.krylov) ~= h.flowopt.eigs.krylov || h.flowopt.eigs.krylov <= 0
            if h.show.n_eig.warning == 1
                waitfor(warndlg('Value must be a positive integer.','Warning'))
            end
            h.show.n_eig.warning = 0;
            h.flowopt.eigs.krylov = floor(abs(h.flowopt.eigs.krylov));
            if h.flowopt.eigs.krylov == 0, h.flowopt.eigs.krylov = 1; end
            set(h.sy.et.kryl,'string',num2str(h.flowopt.eigs.krylov))
        end
    case h.sy.et.maxit
        if floor(h.flowopt.eigs.maxit) ~= h.flowopt.eigs.maxit || h.flowopt.eigs.maxit <= 0
            if h.show.n_eig.warning == 1
                waitfor(warndlg('Value must be a positive integer.','Warning'))
            end
            h.show.n_eig.warning = 0;
            h.flowopt.eigs.maxit = floor(abs(h.flowopt.eigs.maxit));
            if h.flowopt.eigs.maxit == 0, h.flowopt.eigs.maxit = 1; end
            set(h.sy.et.maxit,'string',num2str(h.flowopt.eigs.maxit))
        end
end

guidata(hObject,h)

function sy_et_m(hObject, ~, ~)
h = guidata(hObject);

switch hObject
    case h.sy.et.m_start
        name = 'm_start';

    case h.sy.et.m_delta
        name = 'm_delta';

    case h.sy.et.m_end
        name = 'm_end';
end

if h.flowopt.ax == 1
    if floor(h.(name)) ~= h.(name) || h.(name) < 0
        if h.show.m.warning == 1
            waitfor(warndlg('Value must be 0 or a positive integer.','Warning'))
        end
        h.show.m.warning = 0;
        h.(name) = floor(abs(h.(name)));
        set(hObject,'string',num2str(h.(name)))
    end
else
    if strcmp(name,'m_delta')
        if h.(name) <= 0
            if h.show.m_delta.warning == 1
                waitfor(warndlg('Value must be positive.','Warning'))
            end
            h.show.m_delta.warning = 0;
            h.(name) = 1;
            set(hObject,'string',num2str(h.(name)))
        end
    else
        if h.(name) < 0
            if h.show.m.warning == 1
                waitfor(warndlg('Value must be greater or equal 0.','Warning'))
            end
            h.show.m.warning = 0;
            h.(name) = 0;
            set(hObject,'string',num2str(h.(name)))
        end
    end
end

guidata(hObject,h)

function sy_bg_parameter(hObject, ~, ~)
h = guidata(hObject);

selection = get(get(hObject,'SelectedObject'),'String');

switch selection
    case ' '
        h.dependent = 'delta_T';
        h.T_init = (h.T_d1l+h.T_d2l)/2;
    case '  '
        h.dependent = 'T_d1l';
        h.T_init = h.T_d1l;
        h.equal_T_d2 = {};
    case '   '
        h.dependent = 'T_d2l';
        h.T_init = h.T_d2l;
        h.equal_T_d1 = {};
end

guidata(hObject, h);

function sy_pb_most_d_mode(hObject, ~, ~)
h = guidata(hObject);

mk = 'km'; p = h.flowopt.ax+1; % pointer

% error message if there is no basic state solution
    flg = evalin('base','exist(''b1'',''var'') && isfield(b1,''w'') && residual<1e15');
    if flg == 0
        waitfor(errordlg('Basic state solution not found!','Error 404'))
        return
    end

% error message if invalid range of wavenumbers
    if isempty(h.m_start:h.m_delta:h.m_end)
        waitfor(errordlg('Invalid range of wavenumbers!','Error 007'))
        return
    end

% changing button options
    set(h.sy.st.run,'visible','on','string','running ...')
    set([h.tabs{:} h.sy.et.n_eig h.sy.et.n_eig_c h.sy.et.m_start h.sy.et.m_delta h.sy.et.m_end h.sy.rb.delta_T h.sy.rb.T_d1 h.sy.rb.T_d2 h.sy.pb.run(1) h.sy.pb.run(2) h.sy.pb.load h.sy.pb.save h.sy.pb.plot],'Enable','inactive')
    if h.flowopt.energy == 0
        set([h.sy.rb.delta_T h.sy.rb.T_d1 h.sy.rb.T_d2 h.sy.pb.run(2)],'Enable','off')
    end
    
% show axes for eigenvalues
    axes(h.sy.as.eigs), cla(h.sy.as.eigs,'reset');
    ax = gca; ax.Box = 'on'; xlabel('$s=\Re(\gamma)$ [1/s]','interpreter','latex'); ylabel('$\omega=\Im(\gamma)$ [1/s]','interpreter','latex'); hold on
    ylim([-1 1]); x_lim = xlim; h.hline = plot(x_lim,[0 0],'k-.','HandleVisibility','off');
    if hObject == h.sy.pb.run(2), ce_display_dependent(h,0), end
    axes(h.sy.as.eigs)
    % if this is the second run, plot previous gamma
    if isfield(h,'old')
        h.hline.Visible = 'off';
        i = 0;
        for m = h.m_start:h.m_delta:h.m_end
            i = i+1;
            if isfield(h.old,['gamma_m' strrep(num2str(m),'.','_')])
                plot(complex(h.old.("gamma_m" + strrep(num2str(m),'.','_'))),'o','DisplayName',[mk(p) ' = ' num2str(m)  ' (previous)' ],'color',h.colors(mod(i-1,size(h.colors,1))+1,:))
            end
            lgnd = legend('show'); legend('Location','northwest'); set(lgnd,'Interpreter','latex'); ylim('auto');
        end
        y_lim = ylim; y_lim = max(abs(y_lim))*[-1 1]; ylim(y_lim)
        x_lim = xlim; xlim('manual'); h.hline = plot(x_lim,[0 0],'k-.','HandleVisibility','off');
        if x_lim(1)*x_lim(2)<0
            plot([0 0],y_lim,'k--','HandleVisibility','off')
        end
    end
    
    updateValues_SFM; evalin('base','caseSetup_SFM_V3d2') % problem_specification
    
% run simulation
    evalin('base','gamma_m = []; m_s =  m_start:m_delta:m_end;')
    i = 0;
    for m = h.m_start:h.m_delta:h.m_end
        i = i+1;
        set(h.sy.st.m,'visible','on','string',[mk(p) ' = ' num2str(m) ' ...'])
        pause(1)
                
        assignin('base','m',m);
        evalin('base','flowopt.stability=1; flowopt.m=m;')
        updateValues_SFM; evalin('base','caseSetup_SFM_V3d2') % problem_specification
        evalin('base','grid_generation; block_assembly_V3d2; opts.p=flowopt.eigs.krylov; opts.maxit=flowopt.eigs.maxit; opts.tol=flowopt.eigs.tol;');
        evalin('base','[eigenvectors,gamma] = eigs(A_stability,B_stability,flowopt.eigs.n,''sm'',opts);')
        evalin('base','gamma=diag(gamma,0);gamma=unique([gamma; conj(gamma)]);[~,index_gamma]=sort(real(gamma));gamma=gamma(index_gamma,1);')
        h.gamma_plot = evalin('base','real(gamma)*-1+imag(gamma)*1i');
        
        if h.flowopt.eigs.n_cayley > 0
            evalin('base','r=5; reig=gamma(imag(gamma)<=min(flowopt.tolerance.newton,flowopt.tolerance.residuals)); user=1.2; alpha1=(real(gamma(1,1))*(-user-1)+2*real(reig(r+1,1)))/(-user+1); alpha2=2*real(reig(r+1,1))-alpha1;')
            evalin('base','[eigenvectors_cayley, gamma_cayley,flag]= eigs(A_stability-alpha2*B_stability,A_stability-alpha1*B_stability,flowopt.eigs.n_cayley,''lm'',opts);')
            evalin('base','gamma_cayley=diag(gamma_cayley,0); gamma_cayley(abs(abs(gamma_cayley)-1)<eps)=[]; if flag==1, gamma_cayley(gamma_cayley==0)=[]; end; gamma_cayley=(alpha2-gamma_cayley*alpha1)./(1-gamma_cayley); gamma_cayley=unique([gamma_cayley; conj(gamma_cayley)]); [~,index_gamma_cayley]=sort(real(gamma_cayley)); gamma_cayley=gamma_cayley(index_gamma_cayley,1);')
            h.gamma_plot = [h.gamma_plot; evalin('base','real(gamma_cayley)*-1+imag(gamma_cayley)*1i')];
            evalin('base','if real(gamma_cayley(1,1))<real(gamma(1,1)), eigenvector = eigenvectors_cayley(:,index_gamma_cayley(1,1)); gamma = [gamma_cayley; gamma(flowopt.eigs.n_cayley+1:end)];, else eigenvector = eigenvectors(:,index_gamma(1,1)); end')
        else
            evalin('base','eigenvector = eigenvectors(:,index_gamma(1,1));')
        end
        evalin('base','eval([''old.eigenvector_m_'' strrep(num2str(m),''.'',''_'') ''=eigenvector;'']); old.gamma=gamma(1,1);')

        % plot eigenvalues
        h.gamma = evalin('base','gamma');
        h.("gamma_m" + strrep(num2str(m),'.','_')) = h.gamma;
        h.old.("gamma_m" + strrep(num2str(m),'.','_')) = h.gamma_plot;
        
        h.hline.Visible = 'off'; xlim('auto'); ylim('auto'); hold on
        h.gamma_plot = complex(h.gamma_plot);
        plt = plot(h.gamma_plot,'o','DisplayName',[mk(p) ' = ' num2str(m)],'color',h.colors(mod(i-1,size(h.colors,1))+1,:)); set(plt,'MarkerFaceColor',get(plt,'Color'));
        xlabel('$s=\Re(\gamma)$ [1/s]','interpreter','latex'); ylabel('$\omega=\Im(\gamma)$ [1/s]','interpreter','latex'); lgnd = legend('show'); legend('Location','northwest'); set(lgnd,'Interpreter','latex')
        y_lim = ylim; y_lim = max(abs(y_lim))*[-1 1]; ylim(y_lim)
        x_lim = xlim; xlim('manual'); h.hline = plot(x_lim,[0 0],'k-.','HandleVisibility','off');
        if x_lim(1)*x_lim(2)<0
            if isfield(h,'vline') && isvalid(h.vline), h.vline.Visible = 'off'; end
            h.vline = plot([0 0],y_lim,'k--','HandleVisibility','off');
        end
        
        evalin('base','eval([''eigenvector_'' strrep(num2str(m),''.'',''_'') '' = eigenvector;''])')
        evalin('base','eval([''gamma_'' strrep(num2str(m),''.'',''_'') '' = gamma;''])')
        evalin('base','gamma_m = [gamma_m; min(real(gamma))];')
    end
    
    
    h.hline.Visible = 'off'; xlim('auto'); ylim('auto');
    y_lim = ylim; y_lim = max(abs(y_lim))*[-1 1]; ylim(y_lim)
    x_lim = xlim; xlim('manual'); h.hline = plot(x_lim,[0 0],'k-.','HandleVisibility','off');
    
    evalin('base','[~, idx] = min(real(gamma_m)); m = m_s(idx); gamma = eval([''gamma_'' strrep(num2str(m),''.'',''_'')]); eigenvector = eval([''eigenvector_'' strrep(num2str(m),''.'',''_'')]);')
    evalin('base','disp([''gamma = '', num2str(min(real(gamma)))])')
    
    set(h.sy.st.m,'visible','off')
    h.gamma = evalin('base','gamma'); ce_display_dependent(h);

    % show this message box only if 'most dangerous mode' was pressed
    m_c = evalin('base','m');
    if hObject == h.sy.pb.run(1)
        helpdlg(['Most dangerous wavenumber ' mk(p) ' = ' num2str(m_c) '!'])
        evalin('base','gamma = gamma(1);')
        updateValues_SFM; evalin('base','caseSetup_SFM_V3d2') % problem_specification
        evalin('base','flowopt.stability=0; grid_generation; block_assembly_V3d2')
        evalin('base','iteration = 0; block_assembly_V3d2; perturbation_amplitudes')
        
        set([h.tabs{:} h.sy.et.n_eig h.sy.et.n_eig_c h.sy.et.m_start h.sy.et.m_delta h.sy.et.m_end h.sy.rb.delta_T h.sy.rb.T_d1 h.sy.rb.T_d2 h.sy.pb.run(1) h.sy.pb.run(2) h.sy.pb.load h.sy.pb.save h.sy.pb.plot],'Enable','on')
        set(h.sy.st.run,'visible','off')

        if h.flowopt.energy == 0
            set([h.sy.rb.delta_T h.sy.rb.T_d1 h.sy.rb.T_d2 h.sy.pb.run(2)],'Enable','off')
        end
    end

guidata(hObject,h)

function sy_pb_critical_mode(hObject, ~, ~)
h = guidata(hObject);

% error message if there is no basic state solution
    flg = evalin('base','exist(''b1'',''var'') && isfield(b1,''w'') && residual<1e15');
    if flg == 0
        waitfor(errordlg('Basic state solution not found!','Error 404'))
        return
    elseif isempty(h.m_start:h.m_delta:h.m_end) % error message if invalid range of wavenumbers
        waitfor(errordlg('Invalid range of wavenumbers!','Error 007'))
        return
    else
        set(h.sy.pb.run(2),'Visible','off')
        set(h.sy.pb.stop,'Visible','on','userdata',0)
    end

% find out which of the temperatures are equal because they will remain
% equal during search for critical mode
    T_d1l = str2double(h.b1.bc.rhs.T.r{1}(5:end));
    T_d2l = str2double(h.b1.bc.rhs.T.r{2}(5:end));
    if h.r_c>0 || h.flowopt.ax==0
        T_rc = str2double(h.b1.bc.rhs.T.z{1}(5:end));
    else
        T_rc = NaN;
    end
    
    % initialize dependent parameter
    switch h.dependent
        case 'delta_T'
            h.T_init = abs(T_d2l-T_d1l);
            h.delta_T = h.T_init;
            h.T_0 = (T_d1l+T_d2l)/2;
        case 'T_d1l'
            h.T_init = T_d1l;
        case 'T_d2l'
            h.T_init = T_d2l;
    end

    T_list = {'T_d1l','T_d2l','T_rc'};
    T_s = [T_d1l, T_d2l, T_rc];

    h.equal_T_d1 = {}; h.idx_T_d1 = [];
    h.equal_T_d2 = {}; h.idx_T_d2 = [];
    for i = 1:length(T_s)
        if T_d1l == T_s(i) && ~strcmp(h.dependent,'T_d2l')
            h.equal_T_d1 = [h.equal_T_d1, T_list{i}];
            h.idx_T_d1 = [h.idx_T_d1 i];
        end
        if T_d2l == T_s(i) && ~strcmp(h.dependent,'T_d1l')
            h.equal_T_d2 = [h.equal_T_d2, T_list{i}];
            h.idx_T_d2 = [h.idx_T_d2 i];
        end
    end

% if gamma_m exists, one can skip the first eigenvalue computation
    gamma_flg = evalin('base','exist(''gamma_m'',''var'')');
    if gamma_flg == 0
        guidata(hObject,h); sy_pb_most_d_mode(hObject); h = guidata(hObject);
    else
        set(h.sy.st.run,'visible','on','string','running ...')
    end
    
    h.gamma = evalin('base','gamma');
    h.gamma_m = evalin('base','gamma_m');
    tol = h.flowopt.tolerance.growth;
    
% compute gamma
    if min(real(h.gamma_m))<0
        while min(real(h.gamma_m))<0 && get(h.sy.pb.stop,'userdata')==0
            evalin('base','gammaminus_m = gamma_m;')
            h.alpha_minus = h.(h.dependent);
            if (strcmp(h.dependent,'T_d1l') && h.T_d1l<h.T_d2l) || (strcmp(h.dependent,'T_d2l') && h.T_d1l>h.T_d2l)
                h.(h.dependent) = h.(h.dependent)*1.05;
            else
                h.(h.dependent) = h.(h.dependent)/1.05;
            end
            h.alpha = h.(h.dependent);
            disp([h.dependent ' = ' num2str(h.(h.dependent))])
            h = sy_update_parameter(h);
            
            % basic_state
            ce_basic_state(h)
            if get(h.sy.pb.stop,'userdata')
                return
            end
            
            % most dangerous mode
            guidata(hObject,h); sy_pb_most_d_mode(hObject); h = guidata(hObject);
            h.gamma = evalin('base','gamma');
            h.gamma_m = evalin('base','gamma_m');
        end

        % computation stopped by user
        if get(h.sy.pb.stop,'userdata')
            set(h.sy.st.run,'String','running ...','visible','off')
            set(h.sy.pb.run(2),'Visible','on')
            set(h.sy.pb.stop,'Visible','off')
            helpdlg('Simulation stopped!','Notice')
            set([h.tabs{:} h.sy.et.n_eig h.sy.et.n_eig_c h.sy.et.m_start h.sy.et.m_delta h.sy.et.m_end h.sy.rb.delta_T h.sy.rb.T_d1 h.sy.rb.T_d2 h.sy.pb.run(1) h.sy.pb.run(2) h.sy.pb.load h.sy.pb.save h.sy.pb.plot],'Enable','on')

            % update all edit text boxes according to new temperatures
            et_T_equal = [h.bc.et.axialLB.d1.temperature h.bc.et.axialLB.d2.temperature h.bc.et.radialLB.rc.temperature];
            for i = 1:length(T_s)
                if T_d1l == T_s(i) && ~strcmp(h.dependent,'T_d2l')
                    et_T_equal(i).String = num2str(h.T_d1l);
                end
                if T_d2l == T_s(i) && ~strcmp(h.dependent,'T_d1l')
                    et_T_equal(i).String = num2str(h.T_d2l);
                end
            end
            return
        end
            
        h.alpha_plus_m = h.alpha;
        h.gammaplus_m = h.gamma_m;
        h.gammaminus_m = evalin('base','gammaminus_m');
        [h.gamma_m,h.m_n] = sort(real(h.gammaminus_m));
        h.gamma_m = h.gamma_m(real(h.gamma_m)<0);
        h.gammaplus_m = h.gammaplus_m(h.m_n);
        h.gammaplus_m = h.gammaplus_m(real(h.gamma_m)<0);
        
    else
        while min(real(h.gamma_m))>0 && get(h.sy.pb.stop,'userdata')==0
            evalin('base','gammaplus_m = gamma_m;')
            h.alpha_plus_m = h.(h.dependent);
            if (strcmp(h.dependent,'T_d1l') && h.T_d1l<h.T_d2l) || (strcmp(h.dependent,'T_d2l') && h.T_d1l>h.T_d2l)
                h.(h.dependent) = h.(h.dependent)/1.05;
            else
                h.(h.dependent) = h.(h.dependent)*1.05;
            end
            h.alpha = h.(h.dependent);
            disp([h.dependent ' = ' num2str(h.(h.dependent))])
            h = sy_update_parameter(h);
            
            % basic_state
            ce_basic_state(h)
            if get(h.sy.pb.stop,'userdata')
                return
            end
            
            % most dangerous mode
            guidata(hObject,h); sy_pb_most_d_mode(hObject); h = guidata(hObject);
            h.gamma = evalin('base','gamma');
            h.gamma_m = evalin('base','gamma_m');
            
        end

        % computation stopped by user
        if get(h.sy.pb.stop,'userdata')
            set(h.sy.st.run,'String','running ...','visible','off')
            set(h.sy.pb.run(2),'Visible','on')
            set(h.sy.pb.stop,'Visible','off')
            helpdlg('Simulation stopped!','Notice')
            set([h.tabs{:} h.sy.et.n_eig h.sy.et.n_eig_c h.sy.et.m_start h.sy.et.m_delta h.sy.et.m_end h.sy.rb.delta_T h.sy.rb.T_d1 h.sy.rb.T_d2 h.sy.pb.run(1) h.sy.pb.run(2) h.sy.pb.load h.sy.pb.save h.sy.pb.plot],'Enable','on')

            % update all edit text boxes according to new temperatures
            et_T_equal = [h.bc.et.axialLB.d1.temperature h.bc.et.axialLB.d2.temperature h.bc.et.radialLB.rc.temperature];
            for i = 1:length(T_s)
                if T_d1l == T_s(i) && ~strcmp(h.dependent,'T_d2l')
                    et_T_equal(i).String = num2str(h.T_d1l);
                end
                if T_d2l == T_s(i) && ~strcmp(h.dependent,'T_d1l')
                    et_T_equal(i).String = num2str(h.T_d2l);
                end
            end
            return
        end
            
        h.alpha_minus = h.alpha;
        h.gammaplus_m = evalin('base','gammaplus_m');
        [h.gamma_m,h.m_n] = sort(real(h.gamma_m));
        h.gamma_m = h.gamma_m(real(h.gamma_m)<0);
        h.gammaplus_m = h.gammaplus_m(h.m_n);
        h.gammaplus_m = h.gammaplus_m(real(h.gamma_m)<0);
    end
    
    % Gottlieb10
    h.gamma = evalin('base','gamma');
    h.alpha_plus = h.alpha_plus_m;
    h.m = h.m_n(1)+h.m_start-h.m_delta;
    disp(['m = ', num2str(h.m)])   
    h.gammaminus = h.gamma_m(1);
    h.gammaplus = h.gammaplus_m(1);
    h.gammaplus_m = h.gammaplus_m(2:end);
    h.alpha_n = h.alpha_minus;
    h.gamman = h.gammaminus;
    while (abs(cos(h.alpha_plus)-cos(h.alpha_minus))>=tol || abs(sin(h.alpha_plus)-sin(h.alpha_minus))>=tol) && min(abs(real(h.gamma)))>=tol && get(h.sy.pb.stop,'userdata')==0
        h.D = (h.alpha_minus-h.alpha_plus)/2;
        h.alpha_n = (h.alpha_minus+h.alpha_plus)/2;
        h.alpha = h.alpha_n;
        h.(h.dependent) = h.alpha;
        disp([h.dependent ' = ' num2str(h.(h.dependent))])
        h = sy_update_parameter(h);
        
        % basic_state
        ce_basic_state(h)
        if get(h.sy.pb.stop,'userdata')
            return
        end
            
        % most dangerous mode
        guidata(hObject,h); sy_pb_most_d_mode(hObject); h = guidata(hObject);
        h.gamma = evalin('base','gamma');
        h.gamma_m = evalin('base','gamma_m');
            
        if (abs(cos(h.alpha_plus)-cos(h.alpha_minus))<=tol && abs(sin(h.alpha_plus)-sin(h.alpha_minus))<=tol) || min(abs(real(h.gamma)))<=tol, break, end
        h.gamman = min(real(h.gamma));
        h.a = (h.gammaminus+h.gammaplus-2*h.gamman)/(2*h.D^2);
        h.b = (h.gammaminus-h.gammaplus)/(2*h.D);
        h.alpha = h.alpha_n-2*h.gamman/(h.b*(1+sqrt(1-4*h.a*h.gamman/h.b^2)));
        h.(h.dependent) = h.alpha;
        disp([h.dependent ' = ' num2str(h.(h.dependent))])
        h = sy_update_parameter(h);
        
        % basic_state
        ce_basic_state(h)
        if get(h.sy.pb.stop,'userdata')
            return
        end
            
        % most dangerous mode
        guidata(hObject,h); sy_pb_most_d_mode(hObject); h = guidata(hObject);
        h.gamma = evalin('base','gamma');
        h.gamma_m = evalin('base','gamma_m');
        
        if min(real(h.gamma))<0
            h.alpha_minus = h.alpha;
            h.gammaminus = min(real(h.gamma));
            if h.gamman>0
                h.alpha_plus = h.alpha_n;
                h.gammaplus = h.gamman;
            end
        else
            h.alpha_plus = h.alpha;
            h.gammaplus = min(real(h.gamma));
            if h.gamman<0
                h.alpha_minus = h.alpha_n;
                h.gammaminus = h.gamman;
            end
        end
    end

    % update all edit text boxes according to new temperatures
    et_T_equal = [h.bc.et.axialLB.d1.temperature h.bc.et.axialLB.d2.temperature h.bc.et.radialLB.rc.temperature];
    for i = 1:length(T_s)
        if T_d1l == T_s(i) && ~strcmp(h.dependent,'T_d2l')
            et_T_equal(i).String = num2str(h.T_d1l);
        end
        if T_d2l == T_s(i) && ~strcmp(h.dependent,'T_d1l')
            et_T_equal(i).String = num2str(h.T_d2l);
        end
    end

    % computation stopped by user
    if get(h.sy.pb.stop,'userdata')
        set(h.sy.st.run,'String','running ...','visible','off')
        set(h.sy.pb.run(2),'Visible','on')
        set(h.sy.pb.stop,'Visible','off')
        helpdlg('Simulation stopped!','Notice')
        set([h.tabs{:} h.sy.et.n_eig h.sy.et.n_eig_c h.sy.et.m_start h.sy.et.m_delta h.sy.et.m_end h.sy.rb.delta_T h.sy.rb.T_d1 h.sy.rb.T_d2 h.sy.pb.run(1) h.sy.pb.run(2) h.sy.pb.load h.sy.pb.save h.sy.pb.plot],'Enable','on')

        return
    end
    
    mk = 'km'; p = h.flowopt.ax+1;
    m_c = evalin('base','m');
    warndlg(['Critical wavenumber ' mk(p) '_c = ' num2str(m_c) '!'],'Information',h.opts)
    evalin('base','gamma = gamma(1);')
    updateValues_SFM; evalin('base','caseSetup_SFM_V3d2') % problem_specification
    evalin('base','flowopt.stability=0; grid_generation; block_assembly_V3d2')
    evalin('base','iteration = 0; block_assembly_V3d2; perturbation_amplitudes')
    
set([h.tabs{:} h.sy.et.n_eig h.sy.et.n_eig_c h.sy.et.m_start h.sy.et.m_delta h.sy.et.m_end h.sy.rb.delta_T h.sy.rb.T_d1 h.sy.rb.T_d2 h.sy.pb.run(1) h.sy.pb.run(2) h.sy.pb.load h.sy.pb.save h.sy.pb.plot],'Enable','on')
set([h.sy.st.run h.sy.pb.stop],'visible','off')
set(h.sy.pb.run(2),'Visible','on')

guidata(hObject,h)

function sy_pb_stop(hObject, ~, ~)
h = guidata(hObject);

set(h.sy.pb.stop,'userdata',1)
set(h.sy.st.run,'String','stopping ...')

guidata(hObject,h)