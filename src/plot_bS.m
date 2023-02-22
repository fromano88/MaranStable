function plot_bS

f_plot = figure('units','pixels','menubar','none','name','Post Processing - Basic State',...
    'numbertitle','off','ToolBar','figure','resize','off','Visible','Off');

if ispc % windows
    f_size = [1264 650];
else    % linux
    f_size = [1296 760];
end
f_plot.Position(3:4) = f_size;
movegui(f_plot,'center')

initialPlot_basicState

% create figure
% windows or linux
    % w=1 --> windows
    % w=0 --> linux
    w = ispc;

% SFM or TFM
    h.blocks = evalin('base','blocks');  
    
% save writings
    h.latex = {'interpreter','Latex','FontSize',12};
    R = {'style','rad','units','pix','parent'};
    C = {'style','check','units','pix','parent'};
    E = {'style','edit','units','pix','parent'};
    A = {'unit','pix','visible','off','parent'};
    S = {'style','text','units','pix','HorizontalAlignment','left','parent'};

% main panels: ct = contourplot, lp = lineplot
    h.pl.ct_plot = uipanel('unit','pix','pos',[0 0 988-28*w 760-110*w],'BorderType','none','BackgroundColor',[1 1 1]);
    h.pl.lp_plot = uipanel('unit','pix','pos',[0 0 988-28*w 760-110*w],'BorderType','none','BackgroundColor',[1 1 1],'visible','off');
    h.pl.ct      = uipanel('unit','pix','pos',[988-28*w 33 308-4*w 570-77*w],'BorderType','none');
    h.pl.lp      = uipanel('unit','pix','pos',[988-28*w 33 308-4*w 570-77*w],'BorderType','none','visible','off');
    
% subpanels
    h.ct.pl.domain(1) = uipanel(h.pl.ct,'units','pix','pos',[8 484-67*w 296-4*w  86-10*w],'title','Domain');
    h.ct.pl.quantity  = uipanel(h.pl.ct,'units','pix','pos',[8 370-42*w 296-4*w 110-25*w],'title','Field Quantity');
    h.ct.pl.vector    = uipanel(h.pl.ct,'units','pix','pos',[8 191-30*w 296-4*w 175-12*w],'title','Velocity Vectors');
    h.ct.pl.lines(1)  = uipanel(h.pl.ct,'units','pix','pos',[8   8      296-4*w 179-30*w],'title','Isolines');
    % --------------------------------------------------------------------
    h.lp.pl.line_plot = uipanel(h.pl.lp_plot,'units','pix','pos',[0 0 870-20*w 760-110*w],'BorderType','none','BackgroundColor',[1 1 1]);
    h.lp.pl.outline   = uipanel(h.pl.lp_plot,'units','pix','pos',[870-25*w 8 110 200],'BorderType','none','BackgroundColor',[1 1 1]);
    h.lp.pl.quantity  = uipanel(h.pl.lp,'units','pix','pos',[8 422-49*w 296-4*w  94-21*w],'title','Field Quantity');
    h.lp.pl.domain    = uipanel(h.pl.lp,'units','pix','pos',[8 260-17*w 296-4*w  50-9*w],'title','Domain');
    h.lp.pl.pos(1)    = uipanel(h.pl.lp,'units','pix','pos',[8 126-11*w 296-4*w  130-6*w],'title','Position');
    h.lp.pl.interval  = uipanel(h.pl.lp,'units','pix','pos',[8   8      296-4*w 114-11*w],'title','Plot Interval');
    
% subsubpanels
    h.ct.pl.domain(2) = uipanel(h.ct.pl.domain(1),'units','pix','pos',[170 -5 131 96]);
    h.ct.pl.lines(2)  = uipanel(h.ct.pl.lines(1),'units','pix','pos',[-5 -5 306 88-11*w]);
    h.lp.pl.pos(2)    = uipanel(h.lp.pl.pos(1),'units','pix','pos',[-5 -5 306 89]);
    
% button groups
    h.bg.dimension = uibuttongroup('units','pix','pos',[996-30*w 663-85*w 296-4*w 34-4*w],'SelectionChangedFcn',{@bg_dimension,h});
    h.bg.unity     = uibuttongroup('units','pix','pos',[996-30*w 607-77*w 296-4*w 52-8*w],'title','Units','SelectionChangedFcn',{@bg_unity,h});
    % --------------------------------------------------------------------
    h.lp.bg.number = uibuttongroup(h.pl.lp,'units','pix','pos',[8 519-66*w 296-4*w 50-10*w],'title','Number of Plots','SelectionChangedFcn',{@lp_bg_number,h});
    h.lp.bg.axis   = uibuttongroup(h.pl.lp,'units','pix','pos',[8 368-36*w 296-4*w 50-10*w],'title','Plot versus','SelectionChangedFcn',{@lp_bg_axis,h});
    h.lp.bg.coord  = uibuttongroup(h.pl.lp,'units','pix','pos',[8 314-26*w 296-4*w 50-10*w],'title','Plotting Line','SelectionChangedFcn',{@lp_bg_coord,h});
    
% radio buttons
    h.rb.dimless = uicontrol(R{:},h.bg.dimension,'pos',[12+14*w 5-2*w 116 22],'string','Dimensionless');
    h.rb.dim = uicontrol(R{:},h.bg.dimension,'pos',[177 5-2*w 101 22],'string','Dimensional','value',1);
    h.rb.m = uicontrol(R{:},h.bg.unity,'pos',[26 10-1*w 52 20],'string','m');
    h.rb.dm = uicontrol(R{:},h.bg.unity,'pos',[89 10-1*w 52 20],'string','dm');
    h.rb.cm = uicontrol(R{:},h.bg.unity,'pos',[152 10-1*w 52 20],'string','cm');
    h.rb.mm = uicontrol(R{:},h.bg.unity,'pos',[213 10-1*w 52 20],'string','mm','value',1);
    % --------------------------------------------------------------------
    h.lp.rb.dual = uicontrol(R{:},h.lp.bg.number,'pos',[27 10-4*w 78 19],'string','Dual Plot','value',1);
    h.lp.rb.single = uicontrol(R{:},h.lp.bg.number,'pos',[174 10-4*w 90 19],'string','Single Plot');
    h.lp.rb.zAxis = uicontrol(R{:},h.lp.bg.axis,'pos',[27 10-4*w 78 19],'string','z-Axis','value',1);
    h.lp.rb.rAxis = uicontrol(R{:},h.lp.bg.axis,'pos',[174 10-4*w 80 19],'string','r-Axis');
    h.lp.rb.curved = uicontrol(R{:},h.lp.bg.coord,'pos',[27 10-4*w 120 19],'string','Coordinate Line','value',1);
    h.lp.rb.straight = uicontrol(R{:},h.lp.bg.coord,'pos',[174 10-4*w 80 19],'string','Cartesian');
    
% check boxes
    h.ct.cb.LB = uicontrol(C{:},h.ct.pl.domain(1),'pos',[20 43-5*w 130 20],'String','Liquid Phase','value',1,'call',{@ct_cb_domain,h});
    h.ct.cb.SG = uicontrol(C{:},h.ct.pl.domain(1),'pos',[20 14-5*w 130 20],'String','Gas Phase','value',1,'call',{@ct_cb_domain,h});
    h.ct.cb.equal = uicontrol(C{:},h.ct.pl.domain(2),'pos',[20 32-4*w 130 20],'String','Equal Axis','value',1,'call',{@ct_cb_equal,h});
    h.ct.cb.vector = uicontrol(C{:},h.ct.pl.vector,'pos',[20 128-5*w 50 20],'String','Off','call',{@ct_cb_vector,h});
    h.ct.cb.streamlines = uicontrol(C{:},h.ct.pl.lines(1),'pos',[20 132-22*w 100 20],'String','Streamlines','value',1,'call',{@ct_cb_isolines,h});
    h.ct.cb.temperatureIsolines = uicontrol(C{:},h.ct.pl.lines(2),'pos',[20 51-5*w 160 20],'String','Temperature','call',{@ct_cb_isolines,h});
    if evalin('base','flowopt.energy') == 0
        set(h.ct.cb.temperatureIsolines,'enable','off')
    end
    % --------------------------------------------------------------------
    h.lp.cb.LB = uicontrol(C{:},h.lp.pl.domain,'pos',[20 10-4*w 130 19],'String','Liquid Phase','value',1,'call',{@lp_cb_domain,h});
    h.lp.cb.SG = uicontrol(C{:},h.lp.pl.domain,'pos',[150 10-4*w 130 19],'String','Gas Phase','call',{@lp_cb_domain,h});
    
% edit texts
    h.ct.et.vector(1) = uicontrol(E{:},h.ct.pl.vector,'pos',[195-20*w 91-5*w 60 23],'string',num2str(h.vp.numR),'enable','off','call',{@ct_et_arrows,h});
    h.ct.et.vector(2) = uicontrol(E{:},h.ct.pl.vector,'pos',[195-20*w 55-5*w 60 23],'string',num2str(h.vp.numZ),'enable','off','call',{@ct_et_arrows,h});
    h.ct.et.streamlines = uicontrol(E{:},h.ct.pl.lines(1),'pos',[140 130-22*w 60 23],'string',num2str(h.numStreamlines),'call',{@ct_et_isolines,h});
    h.ct.et.temperatureIsolines = uicontrol(E{:},h.ct.pl.lines(2),'pos',[140 49-5*w 60 23],'string',num2str(h.numTemperatureIsolines),'enable','off','call',{@ct_et_isolines,h});
    % --------------------------------------------------------------------
    h.lp.et.pos = uicontrol(E{:},h.lp.pl.pos(2),'pos',[162 35 60 23],'string',num2str(h.lp.sum),'call',{@lp_et_sr_pos,h});
    h.lp.et.interval(1) = uicontrol(E{:},h.lp.pl.interval,'pos',[42 38-2*w 60 23],'string',num2str(h.lp.zoom(1)),'call',{@lp_et_sr_interval,h});
    h.lp.et.interval(2) = uicontrol(E{:},h.lp.pl.interval,'pos',[42 9-2*w 60 23],'string',num2str(h.lp.zoom(2)),'call',{@lp_et_sr_interval,h});

% popup menus
    h.pm.headline = uicontrol('style','pop','pos',[1120-40*w 714-99*w 110 24],'string',{'Contour Plot','Line Plot'},'call',{@pm_headline,h});
    % --------------------------------------------------------------------
    if evalin('base','flowopt.energy') == 1
        h.ct.pm.quantity = uicontrol(h.ct.pl.quantity,'style','pop','pos',[14 55-15*w 150 24],'string',{'Temperature','u - Velocity','w - Velocity','Velocity Magnitude','Pressure','None'},'call',{@ct_pm_quantity,h});
    else
        h.ct.pm.quantity = uicontrol(h.ct.pl.quantity,'style','pop','pos',[14 55-15*w 150 24],'string',{'u - Velocity','w - Velocity','Velocity Magnitude','Pressure','None'},'call',{@ct_pm_quantity,h});
    end
    h.ct.pm.colormap = uicontrol(h.ct.pl.quantity,'style','pop','pos',[94 14-5*w 115 24],'string',{'Cool2Warm','Cold2Warm','Rainbow','Parula','Thermal','Inferno','Black&White'},'call',{@ct_pm_colormap,h});
    h.ct.pm.color(1) = uicontrol(h.ct.pl.vector,'style','pop','pos',[70 14-5*w 84 24],'string',{'Black','Blue','Green','Red','White'},'enable','off','call',{@ct_pm_color,h});
    h.ct.pm.color(2) = uicontrol(h.ct.pl.lines(1),'style','pop','pos',[117 97-22*w 84 24],'string',{'Black','Blue','Green','Red','White'},'enable','on','call',{@ct_pm_color,h});
    h.ct.pm.color(3) = uicontrol(h.ct.pl.lines(2),'style','pop','pos',[117 14-4*w 84 24],'string',{'Black','Blue','Green','Red','White'},'enable','off','call',{@ct_pm_color,h});
    % --------------------------------------------------------------------
    quantities = ["Temperature","u - Velocity","w - Velocity","Velocity Magnitude","Pressure","Surface Deformation","Normal Heat Flux"];
    if h.dim == 0
        quantities = [quantities, "Biot"];
    end
    if evalin('base','flowopt.energy') == 0
        quantities = quantities(2:end-1);
    end
    h.lp.pm.quantity(1) = uicontrol(h.lp.pl.quantity,'style','pop','pos',[120-14*w 46-16*w 150 24],'string',cellstr(quantities),'call',{@lp_pm_quantity,h});
    h.lp.pm.quantity(2) = uicontrol(h.lp.pl.quantity,'style','pop','pos',[120-14*w 13-10*w 150 24],'string',cellstr(quantities),'value',2,'call',{@lp_pm_quantity,h});
    
% pushbuttons
    h.pb.export = uicontrol('units','pix','visible','on','string','Export',...
        'pos',[1214-32*w 8 76 25],'ForegroundColor',[0.8500 0.3250 0.0980],...
        'call',{@pb_export,h});
    
% slider
    h.lp.sr.pos = uicontrol(h.lp.pl.pos(2),'style','slider','pos',[40 10+2*w 224 20-3*w],'max',2,'value',1,'SliderStep', [1/200, 10/200],'call',{@lp_et_sr_pos,h});
    h.lp.sr.interval(1) = uicontrol(h.lp.pl.interval,'style','slider','pos',[120 40-1*w 154 20-3*w],'call',{@lp_et_sr_interval,h});
    h.lp.sr.interval(2) = uicontrol(h.lp.pl.interval,'style','slider','pos',[120 11 154 20-3*w],'value',1,'call',{@lp_et_sr_interval,h});
    
% axes
    h.ct.as.plot = axes(A{:},h.pl.ct_plot,'OuterPosition',h.pl.ct_plot.Position);
    h.lp.as.plot = axes(A{:},h.lp.pl.line_plot,'OuterPosition',h.lp.pl.line_plot.Position);
    h.lp.as.plot_copy = copy(h.lp.as.plot); % deep copy is needed because 'tiledlayout' deletes the current axis every time 'plot_basic_state_line' is called
    h.lp.as.outline = axes(A{:},h.lp.pl.outline,'Position',[5 5 h.lp.pl.outline.Position(3:4)-6]);
    h.lp.as.pos(1) = axes(A{:},h.lp.pl.pos(1),'pos',[20 91-2*w 60 23]);
    h.lp.as.pos(2) = axes(A{:},h.lp.pl.pos(2),'pos',[25 61 60 23]);
    h.lp.as.pos(3) = axes(A{:},h.lp.pl.pos(2),'pos',[162 35 60 23]);
    if evalin('base','flowopt.ax') == 1
        h.ct.as.vec(1) = axes(A{:},h.ct.pl.vector,'pos',h.ct.et.vector(1).Position); text(-0.63-0.3*w,0.47,'$r$:',h.latex{:});
        h.ct.as.vec(2) = axes(A{:},h.ct.pl.vector,'pos',h.ct.et.vector(2).Position); text(-0.63-0.3*w,0.47,'$z$:',h.latex{:});
    else
        h.ct.as.vec(1) = axes(A{:},h.ct.pl.vector,'pos',h.ct.et.vector(1).Position); text(-0.63-0.3*w,0.47,'$x$:',h.latex{:});
        h.ct.as.vec(2) = axes(A{:},h.ct.pl.vector,'pos',h.ct.et.vector(2).Position); text(-0.63-0.3*w,0.47,'$y$:',h.latex{:});
    end
    h.lp.as.interval = axes(A{:},h.lp.pl.interval,'pos',[11 68-5*w 60 23]);
    h.lp.as.alpha = axes(A{:},h.lp.pl.interval,'pos',h.lp.et.interval(1).Position);
    h.lp.as.beta  = axes(A{:},h.lp.pl.interval,'pos',h.lp.et.interval(2).Position);
    text_axes = [h.lp.as.outline,h.lp.as.pos,h.ct.as.vec,h.lp.as.interval,h.lp.as.alpha,h.lp.as.beta];
    for i = text_axes
        i.Toolbar.Visible = 'off';
    end
    axes(h.lp.as.pos(3)); cla(h.lp.as.pos(3));
    if length(h.blocks) == 2
        text(-0.05,0.47,'$\sum \varepsilon_i= $ ',h.latex{:},'HorizontalAlignment','right')
    else
        set(h.lp.as.pos(3),'pos',[135 50 60 23])
        text(-0.05,0.47,'$\varepsilon= $ ',h.latex{:},'HorizontalAlignment','right')
    end

% static texts
    uicontrol(S{:},f_plot,'pos',[1050-20*w 716-98*w 60-15*w 17],'string','Plot Type');
    uicontrol(S{:},h.ct.pl.quantity,'pos',[14 17-5*w 65 17],'string','Colormap');
    uicontrol(S{:},h.ct.pl.vector,'pos',[20 93-6*w 135-37*w 17],'string','Number of Arrows in');
    uicontrol(S{:},h.ct.pl.vector,'pos',[20 57-6*w 135-37*w 17],'string','Number of Arrows in');
    uicontrol(S{:},h.ct.pl.vector,'pos',[20 17-5*w 40 17],'string','Color');
    uicontrol(S{:},h.ct.pl.lines(1),'pos',[67 100-22*w 40 17],'string','Color');
    uicontrol(S{:},h.ct.pl.lines(2),'pos',[67 17-4*w 40 17],'string','Color');
    % --------------------------------------------------------------------
    h.lp.st.quantity(1) = uicontrol(S{:},h.lp.pl.quantity,'pos',[15 49-15*w 100-10*w 17],'string','Bottom (blue)');
    h.lp.st.quantity(2) = uicontrol(S{:},h.lp.pl.quantity,'pos',[15 16-9*w 90 17],'string','Top (red)');
    h.lp.st.plotting = uicontrol(S{:},h.pl.lp_plot,'pos',[5 5 100 17],'string','plotting ...','BackgroundColor',[1 1 1],'visible','off');
    
% difference between SFM and TFM
    if length(h.blocks) == 1
        set([h.ct.cb.LB h.lp.cb.LB],'enable','inactive')
        set([h.ct.cb.SG h.lp.cb.SG],'value',0,'enable','off')
    end

% update text on axes
    update_text_axes(h)
    
    if length(h.blocks) == 1
        set(h.lp.et.pos,'position',[135 50 60 23])
        set(h.lp.sr.pos,'max',1,'SliderStep',[1/100, 10/100],'pos',[40 18 224 20])
    end
    axes(h.lp.as.alpha);
    text(-0.04,0.47,'$\alpha= $ ',h.latex{:},'HorizontalAlignment','right')
    axes(h.lp.as.beta);
    text(-0.04,0.47,'$\beta= $ ',h.latex{:},'HorizontalAlignment','right')
    
% plot default setup
    axes(h.ct.as.plot); cla(h.ct.as.plot,'reset')
    evalin('base','plot_basic_state_contour')
    %
    h = reset_axes(h);
    evalin('base','plot_basic_state_line')
    %
    axes(h.lp.as.outline); cla(h.lp.as.outline,'reset')
    evalin('base','draw_contours')
    
% change radio button according to flowopt.ax
    if evalin('base','flowopt.ax') == 0
        set(h.lp.rb.rAxis,'string','x-Axis')
        set(h.lp.rb.zAxis,'string','y-Axis')
    else
        set(h.lp.rb.rAxis,'string','r-Axis')
        set(h.lp.rb.zAxis,'string','z-Axis')
    end
    
% showing errors or not
    h.show.arrow.warning    = 1;
    h.show.lines.warning    = 1;
    h.show.interval.error   = 1;
    h.show.interval.warning = 1;
    
% Windows style
    h.opts = struct('WindowStyle','non-modal','Interpreter','tex');

% correct figure size
    pause(0.1)
    f_plot.Position(3:4) = f_size;
    set(f_plot,'Visible','On')
    
guidata(f_plot,h)

function h = reset_axes(h,dlt)
if nargin == 1
    dlt = 0;
end
if dlt == 1 && evalin('base','exist(''ax'',''var'')')
    evalin('base', 'delete(ax)')
    evalin('base', 'clear ax')
end
if ~isgraphics(h.lp.as.plot)
    h.lp.as.plot = copy(h.lp.as.plot_copy); set(h.lp.as.plot,'parent',h.lp.pl.line_plot);
end
axes(h.lp.as.plot); cla(h.lp.as.plot,'reset');

function update_text_axes(h)
p = evalin('base','flowopt.ax + 1;'); % pointer
rx = 'xr'; zy = 'yz'; % coordinate: radial, axial

axes(h.lp.as.pos(1)); cla(h.lp.as.pos(1));
if length(h.blocks) == 2
    if strcmp(h.lp.axis,'z')
        if evalin('base','r_c') == 0
            if strcmp(h.lp.coord,'straight')
                text(0,0.4,['$' rx(p) '=\varepsilon_1 ' rx(p) '_i + \varepsilon_2 (' rx(p) '_o - ' rx(p) '_i)$'],h.latex{:})
            else
                text(0,0.4,['$' rx(p) '=\varepsilon_1 h(' zy(p) ') + \varepsilon_2 (' rx(p) '_o - h(' zy(p) '))$'],h.latex{:})
            end
        else
            if strcmp(h.lp.coord,'straight')
                text(-0.1,0.4,['$' rx(p) '=' rx(p) '_c + \varepsilon_1 (' rx(p) '_i-' rx(p) '_c) + \varepsilon_2 (' rx(p) '_o - ' rx(p) '_i)$'],h.latex{:})
            else
                text(-0.1,0.4,['$' rx(p) '=' rx(p) '_c + \varepsilon_1 (h(' zy(p) ')-' rx(p) '_c) + \varepsilon_2 (' rx(p) '_o - h(' zy(p) '))$'],h.latex{:})
            end
        end
    else
        text(0,0.4,['$' zy(p) '=\varepsilon_1 d_1 + \varepsilon_2 d + \varepsilon_3 d_2$'],h.latex{:})
    end
else
    if strcmp(h.lp.axis,'z')
        if evalin('base','r_c') == 0
            if strcmp(h.lp.coord,'straight')
                text(0,0.4,['$' rx(p) '=\varepsilon ' rx(p) '_i$'],h.latex{:})
            else
                text(0,0.4,['$' rx(p) '=\varepsilon h(' zy(p) ')$'],h.latex{:})
            end
        else
            if strcmp(h.lp.coord,'straight')
                text(-0.1,0.4,['$' rx(p) '=' rx(p) '_c + \varepsilon (' rx(p) '_i-' rx(p) '_c)$'],h.latex{:})
            else
                text(-0.1,0.4,['$' rx(p) '=' rx(p) '_c + \varepsilon (h(' zy(p) ')-' rx(p) '_c)$'],h.latex{:})
            end
        end
    else
        text(0,0.4,['$' zy(p) '=\varepsilon d$'],h.latex{:})
    end
end

axes(h.lp.as.pos(2)); cla(h.lp.as.pos(2));
if length(h.blocks) == 2
    if strcmp(h.lp.axis,'r')
        text(0,0.47,   ['$\varepsilon_1= $ ' num2str(1-h.lp.e_z(1))],  h.latex{:})
        text(1.64,0.47,['$\varepsilon_2= $ ' num2str(0.5-h.lp.e_z(2))],h.latex{:})
        text(3.28,0.47,['$\varepsilon_3= $ ' num2str(-h.lp.e_z(3))],    h.latex{:})
    else
        text(0,0.47,   ['$\varepsilon_1= $ ' num2str(h.lp.e_r(1))],h.latex{:})
        text(1.64,0.47,['$\varepsilon_2= $ ' num2str(h.lp.e_r(2))],h.latex{:})
    end
end

% interval
axes(h.lp.as.interval); cla(h.lp.as.interval);
if strcmp(h.lp.axis,'z')
    text(0,0.5,['$[' zy(p) '_\mathrm{min}\! +\! \alpha (' zy(p) '_\mathrm{max}\!-\! ' zy(p) '_\mathrm{min});\; ' zy(p) '_\mathrm{min}\! +\! \beta (' zy(p) '_\mathrm{max}\!-\!' zy(p) '_\mathrm{min}) ]$'],'interpreter','Latex','FontSize',11)
else
    text(0,0.5,['$[' rx(p) '_\mathrm{min}\! +\! \alpha (' rx(p) '_\mathrm{max}\!-\! ' rx(p) '_\mathrm{min});\; ' rx(p) '_\mathrm{min}\! +\! \beta (' rx(p) '_\mathrm{max}\!-\!' rx(p) '_\mathrm{min}) ]$'],'interpreter','Latex','FontSize',11)
end

function pb_export(hObject, ~, ~)
h = guidata(hObject);

str = get(h.pm.headline,'string'); val = get(h.pm.headline,'value');
quantities = get(h.lp.pm.quantity(1),'string');

% remove surface quantities
    if strcmp(h.lp.axis,'z') && (h.lp.sum ~=1 || strcmp(h.lp.coord,'straight'))
        quantities = quantities(1:end-(3-h.dim));
    end

switch str{val}
    case 'Contour Plot'
        export_vtk('bS')
        
    case 'Line Plot'
        export_dat(quantities)
end

function pm_headline(hObject, ~, ~)
h = guidata(hObject);

str=get(hObject, 'String'); val=get(hObject, 'Value');

switch str{val}
    case 'Contour Plot'
        set([h.pl.ct h.pl.ct_plot],'visible','on');
        set([h.pl.lp h.pl.lp_plot],'visible','off');
        set(h.pb.export,'ForegroundColor',[0.8500 0.3250 0.0980])
        
    case 'Line Plot'
        set([h.pl.ct h.pl.ct_plot],'visible','off');
        set([h.pl.lp h.pl.lp_plot],'visible','on');
        set(h.pb.export,'ForegroundColor',[0 0.4470 0.7410])

end

guidata(hObject,h)

function bg_dimension(hObject, ~, ~)
h = guidata(hObject);

selection = get(get(hObject,'SelectedObject'),'String');
quantities = get(h.lp.pm.quantity(1),'string');

switch selection
    case 'Dimensionless'
        waitfor(warndlg({'Scalings:','     - length: d','     - temperature: \DeltaT','     - velocity: \gamma\DeltaT / \mu','     - pressure: \gamma\DeltaT / d','','d ... length of liquid bridge','\gamma ... negative linear Taylor coefficient of surface tension','\mu ... dynamic viscosity of liquid','\DeltaT ... temperature difference T_{hot}-T_{cold}'},'Note',h.opts))
        set([h.rb.m h.rb.dm h.rb.cm h.rb.mm],'Enable','off','Value',0)
        h.dim = 0;
        if strcmp(h.lp.axis,'z')
            quantities = [quantities; 'Biot'];
        end
 
    case 'Dimensional'
        set([h.rb.m h.rb.dm h.rb.cm h.rb.mm],'Enable','on','Value',0)
        h.dim = 1;
        if strcmp(h.unity,'m')
            set(h.rb.m,'value',1)
        elseif strcmp(h.unity,'dm')
            set(h.rb.dm,'value',1)
        elseif strcmp(h.unity,'cm')
            set(h.rb.cm,'value',1)
        elseif strcmp(h.unity,'mm')
            set(h.rb.mm,'value',1)
        end
        if strcmp(h.lp.axis,'z')
            quantities = quantities(1:end-1);
        end
end

% update text on axes
    update_text_axes(h)
    
% update popupmenu of lp
    set(h.lp.pm.quantity(:),'string',quantities)

% update plots
    assignin('base','tmp',h.dim);
    evalin('base','dim=tmp;');
    % ------------------------------------------------------------
    axes(h.ct.as.plot); cla(h.ct.as.plot,'reset')
    evalin('base','plot_basic_state_contour')

    h = reset_axes(h);
    evalin('base','plot_basic_state_line')
    
guidata(hObject,h)

function bg_unity(hObject, ~, ~)
h = guidata(hObject);

selection = get(get(h.bg.unity,'SelectedObject'),'String');

switch selection
    case 'm'
        h.unity = 'm';
    case 'dm'
        h.unity = 'dm';
    case 'cm'
        h.unity = 'cm';
    case 'mm'
        h.unity = 'mm';
end

% update value
    assignin('base','tmp',h.unity);
    evalin('base','unity=tmp;');

% update plots
    axes(h.ct.as.plot); cla(h.ct.as.plot,'reset')
    evalin('base','plot_basic_state_contour')

    h = reset_axes(h);
    evalin('base','plot_basic_state_line')

guidata(hObject,h)

function ct_cb_domain(hObject, ~, ~)
h = guidata(hObject);

LB = get(h.ct.cb.LB,'value'); SG = get(h.ct.cb.SG,'value');

if LB + SG == 0
    set(hObject,'value',1)
    return
elseif LB + SG == 2
    h.ct.domain = [1 1];
elseif LB == 1
    h.ct.domain = [1 0];
else
    h.ct.domain = [0 1];
end

% update value
    assignin('base','tmp',h.ct.domain);
    evalin('base','ct.domain=tmp;');

% update plot
    axes(h.ct.as.plot); cla(h.ct.as.plot,'reset')
    evalin('base','plot_basic_state_contour')

% update number of shown arrows in terms of new domain
    if get(h.ct.cb.vector,'value')
        if h.ct.domain(1) == 0 % gas
            h.vp.n_r = evalin('base','round((b2.geom.r(end) - b2.geom.r(1))/(max(geom.r)-min(geom.r))*vp.numR)');
            set(h.ct.et.vector(1),'string',num2str(h.vp.n_r))
            set(h.ct.et.vector(2),'string',num2str(h.vp.numZ))
        elseif h.ct.domain(2) == 0 % liquid
            h.vp.n_r = evalin('base','round((b1.geom.r(end) - b1.geom.r(1))/(max(geom.r)-min(geom.r))*vp.numR)');
            h.vp.n_z = evalin('base','round((b1.geom.z(end) - b1.geom.z(1))/(max(geom.z)-min(geom.z))*vp.numZ)');
            set(h.ct.et.vector(1),'string',num2str(h.vp.n_r))
            if h.vp.n_z < 1
                set(h.ct.et.vector(2),'string',num2str(0))
            else
                set(h.ct.et.vector(2),'string',num2str(h.vp.n_z-1))
            end
        else
            set(h.ct.et.vector(1),'string',num2str(h.vp.numR))
            set(h.ct.et.vector(2),'string',num2str(h.vp.numZ))
        end
    end

% update text on axes
    update_text_axes(h)
    
guidata(hObject,h)

function ct_cb_equal(hObject, ~, ~)
h = guidata(hObject);

h.ct.equal = get(hObject,'value');

% update value
    assignin('base','tmp',h.ct.equal);
    evalin('base','ct.equal=tmp;');

% update plot
    axes(h.ct.as.plot); cla(h.ct.as.plot,'reset')
    evalin('base','plot_basic_state_contour')

guidata(hObject,h)

function ct_pm_quantity(hObject, ~, ~)
h = guidata(hObject);

str=get(hObject, 'String'); val=get(hObject, 'Value');

h.ct.draw = 1;
set(h.ct.pm.colormap,'enable','on')

switch str{val}
    case 'Temperature'
        h.ct.quantity = 'T';
        
    case 'w - Velocity'
        h.ct.quantity = 'w';
        
    case 'u - Velocity'
        h.ct.quantity = 'u';
        
    case 'Velocity Magnitude'
        h.ct.quantity = 'norm_u';
        
    case 'Pressure'
        h.ct.quantity = 'p';
        
    case 'None'
        h.ct.draw = 0;
        set(h.ct.pm.colormap,'enable','off')
end

% update values
    assignin('base','tmp',h.ct.quantity);
    evalin('base','ct.quantity=tmp;');
    % --------------------------------------------------------- %
    assignin('base','tmp',h.ct.draw);
    evalin('base','ct.draw=tmp;');

% update plot
    axes(h.ct.as.plot); cla(h.ct.as.plot,'reset')
    evalin('base','plot_basic_state_contour')

guidata(hObject,h)

function ct_pm_colormap(hObject, ~, ~)
h = guidata(hObject);

str=get(hObject,'String'); val=get(hObject,'Value');

switch str{val}
    case 'Cool2Warm'
        h.ct.colormap = evalin('base','myColormap(''Cool2Warm'');');
        
    case 'Cold2Warm'
        h.ct.colormap = evalin('base','myColormap(''BlueWhiteRed'');');
        
    case 'Black&White'
        h.ct.colormap = 'bone(256)';
        
    case 'Rainbow'
        h.ct.colormap = evalin('base','myColormap(''Rainbow'');');
        
    case 'Parula'
        h.ct.colormap = 'parula(256)';
        
    case 'Thermal'
        h.ct.colormap = evalin('base','myColormap(''Thermal'');');
        
    case 'Inferno'
        h.ct.colormap = evalin('base','myColormap(''Thermal2'');');

end

% update values
    assignin('base','tmp',h.ct.colormap);
    evalin('base','ct.colormap=tmp;');

% update plot
    axes(h.ct.as.plot); cla(h.ct.as.plot,'reset')
    evalin('base','plot_basic_state_contour')

guidata(hObject,h)

function ct_cb_vector(hObject, ~, ~)
h = guidata(hObject);

h.vp.draw = get(hObject,'Value');
if h.vp.draw == 1
    set(hObject,'string','On')
    set([h.ct.et.vector(1) h.ct.et.vector(2) h.ct.pm.color(1)],'Enable','on')
else
    set(hObject,'string','Off')
    set([h.ct.et.vector(1) h.ct.et.vector(2) h.ct.pm.color(1)],'Enable','off')
end

% update values
    assignin('base','tmp',h.vp.draw);
    evalin('base','vp.draw=tmp;');
    
% update plot
    axes(h.ct.as.plot); cla(h.ct.as.plot,'reset')
    evalin('base','plot_basic_state_contour')
    
% update number of shown arrows in terms of new domain
    if get(h.ct.cb.vector,'value')
        if h.ct.domain(1) == 0 % gas
            h.vp.n_r = evalin('base','round((b2.geom.r(end) - b2.geom.r(1))/(max(geom.r)-min(geom.r))*vp.numR)');
            set(h.ct.et.vector(1),'string',num2str(h.vp.n_r))
            set(h.ct.et.vector(2),'string',num2str(h.vp.numZ))
        elseif h.ct.domain(2) == 0 % liquid
            h.vp.n_r = evalin('base','round((b1.geom.r(end) - b1.geom.r(1))/(max(geom.r)-min(geom.r))*vp.numR)');
            h.vp.n_z = evalin('base','round((b1.geom.z(end) - b1.geom.z(1))/(max(geom.z)-min(geom.z))*vp.numZ)');
            set(h.ct.et.vector(1),'string',num2str(h.vp.n_r))
            if h.vp.n_z < 1
                set(h.ct.et.vector(2),'string',num2str(0))
            else
                set(h.ct.et.vector(2),'string',num2str(h.vp.n_z-1))
            end
        else
            set(h.ct.et.vector(1),'string',num2str(h.vp.numR))
            set(h.ct.et.vector(2),'string',num2str(h.vp.numZ))
        end
    end

guidata(hObject,h)

function ct_et_arrows(hObject, ~, ~)
h = guidata(hObject);

lowerLimit=-inf; upperLimit=inf; default=h.vp.numR;

[numArrow,~] = getValue(lowerLimit,upperLimit,default,hObject,[],h.show.arrow.warning);
    
if floor(numArrow) ~= numArrow || numArrow <= 0
    if h.show.arrow.warning == 1
        waitfor(warndlg('Value must be a positive integer.','Warning'))
    end
    h.show.arrow.warning = 0;
    numArrow = floor(norm(numArrow));
    set(hObject,'string',num2str(numArrow))
end

LB = get(h.ct.cb.LB,'value'); SG = get(h.ct.cb.SG,'value');

switch hObject
    case h.ct.et.vector(1)
        h.vp.numR = numArrow;
        
        % compute number of arrows with respect to whole geometry
        if LB == 1 && SG == 0
            h.vp.numR = evalin('base','round((max(geom.r)-min(geom.r))/(b1.geom.r(end) - b1.geom.r(1)))')*h.vp.numR;
        elseif LB == 0 && SG == 1
            h.vp.numR = evalin('base','round((max(geom.r)-min(geom.r))/(b2.geom.r(end) - b2.geom.r(1)))')*h.vp.numR;
        end
        % update value
        assignin('base','tmp',h.vp.numR);
        evalin('base','vp.numR=tmp;');
        
    case h.ct.et.vector(2)
        h.vp.numZ = numArrow;
        
        % compute number of arrows with respect to whole geometry
        if LB == 1 && SG == 0
            h.vp.numZ = evalin('base','round((max(geom.z)-min(geom.z))/(b1.geom.z(end) - b1.geom.z(1)))')*(h.vp.numZ+1);
        end
        % update value
        assignin('base','tmp',h.vp.numZ-1);
        evalin('base','vp.numZ=tmp;');
end
    
% update plot
    axes(h.ct.as.plot); cla(h.ct.as.plot,'reset')
    evalin('base','plot_basic_state_contour')

guidata(hObject,h)

function ct_pm_color(hObject, ~, ~)
h = guidata(hObject);

str=get(hObject,'String'); val=get(hObject,'Value');

switch str{val}
    case 'Black'
        clr = '''k''';
    case 'Blue'
        clr = 'mycolor.blue';
    case 'Green'
        clr = 'mycolor.green';
    case 'Red'
        clr = 'mycolor.red';
    case 'White'
        clr = '''w''';
end

% update value
    assignin('base','tmp',clr);
    switch hObject
        case h.ct.pm.color(1)
            evalin('base','vp.color=eval(tmp);');
        case h.ct.pm.color(2)
            evalin('base','ct.color.sf=eval(tmp);');
        case h.ct.pm.color(3)
            evalin('base','ct.color.T=eval(tmp);');
    end

% update plot
    axes(h.ct.as.plot); cla(h.ct.as.plot,'reset')
    evalin('base','plot_basic_state_contour')

guidata(hObject,h)

function ct_cb_isolines(hObject, ~, ~)
h = guidata(hObject);

switch hObject
    case h.ct.cb.streamlines
        h.streamlines = get(hObject,'Value');

        if get(hObject,'Value')
            set([h.ct.et.streamlines h.ct.pm.color(2)],'enable','on')
        else
            set([h.ct.et.streamlines h.ct.pm.color(2)],'enable','off')
        end
        % update value
        assignin('base','tmp',h.streamlines);
        evalin('base','ct.streamlines=tmp;');
    
    case h.ct.cb.temperatureIsolines
        h.temperatureIsolines = get(hObject,'Value');

        if get(hObject,'Value')
            set([h.ct.et.temperatureIsolines h.ct.pm.color(3)],'enable','on')
        else
            set([h.ct.et.temperatureIsolines h.ct.pm.color(3)],'enable','off')
        end
        % update value
        assignin('base','tmp',h.temperatureIsolines);
        evalin('base','ct.isolines=tmp;');        
end
     
% update plot
    axes(h.ct.as.plot); cla(h.ct.as.plot,'reset')
    evalin('base','plot_basic_state_contour')

guidata(hObject,h)

function ct_et_isolines(hObject, ~, ~)
h = guidata(hObject);

lowerLimit=-inf; upperLimit=inf; default=h.vp.numR;

[numLines,~] = getValue(lowerLimit,upperLimit,default,hObject,[],h.show.lines.warning);
    
if floor(numLines) ~= numLines || numLines <= 0
    if h.show.lines.warning == 1
        waitfor(warndlg('Value must be a positive integer.','Warning'))
    end
    h.show.lines.warning = 0;
    numLines = floor(norm(numLines));
    set(hObject,'string',num2str(numLines))
end

switch hObject
    case h.ct.et.streamlines
        h.numStreamlines = numLines;
        % update value
        assignin('base','tmp',h.numStreamlines);
        evalin('base','ct.numStreamlines=tmp;');
        
    case h.ct.et.temperatureIsolines
        h.numTemperatureIsolines = numLines;
        % update value
        assignin('base','tmp',h.numTemperatureIsolines);
        evalin('base','ct.numIsolines=tmp;');
end

% update plot
    axes(h.ct.as.plot); cla(h.ct.as.plot,'reset')
    evalin('base','plot_basic_state_contour')

guidata(hObject,h)






% ------------------ Lineplot -----------------------
function lp_bg_number(hObject, ~, ~)
h = guidata(hObject);

selection = get(get(hObject,'SelectedObject'),'String');
switch selection
    case 'Dual Plot'
        h.lp.n = 2;
        set(h.lp.pm.quantity(2),'Enable','on')
    case 'Single Plot'
        h.lp.n = 1;
        set(h.lp.pm.quantity(2),'Enable','off')
        
        if ismember(h.lp.quantity(2), {'h', 'n_nabla_T', 'Bi'})
            h.lp.quantity(2) = "w";
            str = get(h.lp.pm.quantity(2),'string');
            set(h.lp.pm.quantity(2),'value',find(contains(str,'u - Velocity')))
        end
        if ~ismember(h.lp.quantity(1), {'h', 'n_nabla_T', 'Bi'})
            set([h.lp.rb.rAxis h.lp.rb.zAxis h.lp.rb.curved h.lp.rb.straight h.lp.cb.LB h.lp.cb.SG h.lp.sr.pos h.lp.et.pos],'enable','on')
            if length(h.blocks) == 1
                set(h.lp.cb.LB,'enable','inactive')
                set(h.lp.cb.SG,'enable','off')
            end
        end
end       

% update value
    assignin('base','tmp',h.lp.n);
    evalin('base','lp.n=tmp;');
    % ----------------------------------------------------- %
    assignin('base','tmp',h.lp.quantity(2));
    evalin('base','lp.quantity(2)=tmp;');
    
% plot
    h = reset_axes(h,h.lp.axis=='r');
    evalin('base','plot_basic_state_line')

guidata(hObject,h)

function lp_pm_quantity(hObject, ~, ~)
h = guidata(hObject);

str=get(hObject,'String'); val=get(hObject,'Value');

switch str{val}
    case 'Temperature'
        quantity = "T";
    case 'w - Velocity'
        quantity = "w";
    case 'u - Velocity'
        quantity = "u";
    case 'Velocity Magnitude'
        quantity = "norm_u";
    case 'Pressure'
        quantity = "p";
    case 'Surface Deformation'
        quantity = "h";
    case 'Normal Heat Flux'
        quantity = "n_nabla_T";
    case 'Biot'
        quantity = "Bi";
end

assignin('base','tmp',quantity);
switch hObject
    case h.lp.pm.quantity(1)
        h.lp.quantity(1) = quantity;
        evalin('base','lp.quantity(1)=tmp;');
    case h.lp.pm.quantity(2)
        h.lp.quantity(2) = quantity;
        evalin('base','lp.quantity(2)=tmp;');
end

if max(ismember(h.lp.quantity(:), {'h', 'n_nabla_T', 'Bi'})) == 1
    h.lp.domain = [1 0];
    h.lp.e_r    = [1 0];
    h.lp.sum    = sum(h.lp.e_r);
    h.lp.coord  = 'curved';
    set([h.lp.cb.LB h.lp.rb.curved h.lp.sr.pos],'value',1)
    set(h.lp.cb.SG,'value',0)
    set(h.lp.et.pos,'string','1')
    set([h.lp.rb.rAxis h.lp.rb.zAxis h.lp.rb.curved h.lp.rb.straight h.lp.cb.LB h.lp.cb.SG h.lp.sr.pos h.lp.et.pos],'enable','off')
    
    % update values
    assignin('base','tmp',h.lp.domain);
    evalin('base','lp.domain=tmp;');
    % ----------------------------------------------------- %
    assignin('base','tmp',h.lp.coord);
    evalin('base','lp.coord=tmp;');
    % ----------------------------------------------------- %
    assignin('base','tmp',h.lp.e_r);
    evalin('base','lp.e_r=tmp;');
    
    % update text on axes
    update_text_axes(h)
    
else
    set([h.lp.rb.rAxis h.lp.rb.zAxis h.lp.rb.curved h.lp.rb.straight h.lp.cb.LB h.lp.cb.SG h.lp.sr.pos h.lp.et.pos],'enable','on')
    
end

if length(h.blocks) == 1
    set(h.lp.cb.LB,'enable','inactive')
    set(h.lp.cb.SG,'enable','off')
end

% plot
    h = reset_axes(h);
    evalin('base','plot_basic_state_line')

guidata(hObject,h)

function lp_bg_axis(hObject, ~, ~)
h = guidata(hObject);

selection  = get(get(hObject,'SelectedObject'),'String');
quantities = get(h.lp.pm.quantity(1),'string');

switch selection
    case {'r-Axis','x-Axis'}
        h.lp.axis = 'r';
        set(h.lp.st.quantity(1),'String','Left (blue)')
        set(h.lp.st.quantity(2),'String','Right (red)')
        % =================================================================
        h.lp.coord = 'straight';
        set(h.lp.rb.straight,'value',1)
        set(h.lp.rb.curved,'enable','off')
        % =================================================================
        if length(h.blocks) == 2
            set(h.lp.sr.pos,'min',-1.5,'max',1.5,'SliderStep',[1/300, 10/300])
        else
            set(h.lp.sr.pos,'min',-0.5,'max',0.5,'SliderStep',[1/100, 10/100])
        end
        h.lp.e_z = [1 0.5 0];
        if length(h.blocks) == 2
            h.lp.domain = [1 1];
            set([h.lp.cb.LB h.lp.cb.SG],'value',1)
        end
        % =================================================================
        h.lp.fsum = 1.5-sum(h.lp.e_z); 
        set(h.lp.et.pos,'string',num2str(h.lp.fsum))
        set(h.lp.sr.pos,'value',h.lp.fsum)
        % =================================================================
        if evalin('base','flowopt.energy') == 1
            quantities = quantities(1:5);
        else
            quantities = quantities(1:4);
        end
        set([h.lp.pm.quantity(1) h.lp.pm.quantity(2)],'string',quantities);
        
    case {'z-Axis','y-Axis'}
        h.lp.axis = 'z';
        set(h.lp.st.quantity(1),'String','Bottom (blue)')
        set(h.lp.st.quantity(2),'String','Top (red)')
        % =================================================================
        h.lp.coord = 'curved';
        set(h.lp.rb.curved,'value',1,'enable','on')
        % =================================================================
        set(h.lp.sr.pos,'min',0,'max',length(h.blocks),'SliderStep', [1/(length(h.blocks)*100) , 10/(length(h.blocks)*100)])
        h.lp.e_r = [1 0];
        h.lp.domain = [1 0];
        set(h.lp.cb.LB,'value',1); set(h.lp.cb.SG,'value',0)
        % =================================================================
        h.lp.sum = sum(h.lp.e_r);
        set(h.lp.et.pos,'string',num2str(h.lp.sum))
        set(h.lp.sr.pos,'value',h.lp.sum)
        % =================================================================
        quantities = [quantities; 'Surface Deformation'];
        if evalin('base','flowopt.energy') == 1
            quantities = [quantities; 'Normal Heat Flux'];
            if h.dim == 0
                quantities = [quantities; 'Biot'];
            end
        end
        set([h.lp.pm.quantity(1) h.lp.pm.quantity(2)],'string',quantities);
end

% update text on axes
    update_text_axes(h)

% update value
    assignin('base','tmp',h.lp.axis);
    evalin('base','lp.axis=tmp;');
    assignin('base','tmp',h.lp.domain);
    evalin('base','lp.domain=tmp;');
    assignin('base','tmp',h.lp.coord);
    evalin('base','lp.coord=tmp;');
    if strcmp(h.lp.axis,'r')
        assignin('base','tmp',h.lp.e_z);
        evalin('base','lp.e_z=tmp;');
    else
        assignin('base','tmp',h.lp.e_r);
        evalin('base','lp.e_r=tmp;');
    end    

% plot
    h = reset_axes(h,1);
    evalin('base','plot_basic_state_line')
    axes(h.lp.as.outline); cla(h.lp.as.outline,'reset')
    evalin('base','draw_contours')

guidata(hObject,h)

function lp_bg_coord(hObject, ~, ~)
h = guidata(hObject);
set(h.lp.st.plotting,'Visible','on')

selection  = get(get(hObject,'SelectedObject'),'String');

switch selection
    case 'Coordinate Line'
        h.lp.coord = 'curved';

        if h.lp.e_r(2) > 0
            h.lp.domain = [0 1];
        else
            h.lp.domain = [1 0];
        end

    case 'Cartesian'
        h.lp.coord = 'straight';

        % determine position of h_fs
        h_fs = evalin('base','b1.R.v(end,:)');
        [r_c, r_i] = evalin('base','deal(r_c,r_i)');
        if length(h.blocks) == 2
            r_o = evalin('base','r_o');
        else
            r_o = r_i;
        end
        pos = r_c + h.lp.e_r(1)*(r_i-r_c) + h.lp.e_r(2)*(r_o-r_i);

        if prod(pos<=h_fs) == 1 % only LB
            h.lp.domain = [1 0];
        elseif prod(pos>h_fs) == 1 % only SG
            h.lp.domain = [0 1];
        else % LB or both
            if max(abs(h_fs-r_i)/r_i) < 1e-12
                h.lp.domain = [1 0];
            else 
                h.lp.domain = [1 1];
            end
        end
        
end

set(h.lp.cb.LB,'value',h.lp.domain(1))
set(h.lp.cb.SG,'value',h.lp.domain(2))

% update value
    assignin('base','tmp',h.lp.domain);
    evalin('base','lp.domain=tmp;');
    assignin('base','tmp',h.lp.coord);
    evalin('base','lp.coord=tmp;');

% update text on axes
    update_text_axes(h)

% plot
    h = reset_axes(h,1);
    evalin('base','plot_basic_state_line')
    axes(h.lp.as.outline); cla(h.lp.as.outline,'reset')
    evalin('base','draw_contours')

set(h.lp.st.plotting,'Visible','off')
guidata(hObject,h)

function lp_cb_domain(hObject, ~, ~)
h = guidata(hObject);
set(h.lp.st.plotting,'Visible','on')

% determine position of h_fs
    h_fs = evalin('base','b1.R.v(end,:)');
    [r_c, r_i] = evalin('base','deal(r_c,r_i)');
    if length(h.blocks) == 2
        r_o = evalin('base','r_o');
    else
        r_o = r_i;
    end
    pos = r_c + h.lp.e_r(1)*(r_i-r_c) + h.lp.e_r(2)*(r_o-r_i);

LB = get(h.lp.cb.LB,'value'); SG = get(h.lp.cb.SG,'value');

if LB + SG == 0
    set(hObject,'value',1)
    return
elseif LB + SG == 2
    h.lp.domain = [1 1];
    if strcmp(h.lp.axis,'z')
        switch hObject
            case h.lp.cb.LB
                if strcmp(h.lp.coord,'curved') || prod(pos>h_fs) == 1 % only SG
                    h.lp.domain = [1 0];
                    set(h.lp.cb.SG,'value',0)
                    h.lp.e_r = [0.5 0];
                end
                
            case h.lp.cb.SG
                if strcmp(h.lp.coord,'curved') || prod(pos<h_fs) == 1 % only LB
                    h.lp.domain = [0 1];
                    set(h.lp.cb.LB,'value',0)
                    h.lp.e_r = [1 0.5];
                end
        end
        h.lp.sum = sum(h.lp.e_r);
        set(h.lp.et.pos,'string',num2str(h.lp.sum))
        set(h.lp.sr.pos,'value',h.lp.sum)
    else
        if h.lp.e_z(2) == 0 || h.lp.e_z(2) == 1
            h.lp.e_z = [1 0.5 0];
            h.lp.fsum = 1.5-sum(h.lp.e_z);
            set(h.lp.et.pos,'string',num2str(h.lp.fsum))
            set(h.lp.sr.pos,'value',h.lp.fsum)
        end
    end
    
elseif LB == 1
    h.lp.domain = [1 0];
else
    h.lp.domain = [0 1];
end

% update value
    assignin('base','tmp',h.lp.domain);
    evalin('base','lp.domain=tmp;');
    if strcmp(h.lp.axis,'z')
        assignin('base','tmp',h.lp.e_r);
        evalin('base','lp.e_r=tmp;');
    else
        assignin('base','tmp',h.lp.e_z);
        evalin('base','lp.e_z=tmp;');
    end
    
% update text on axes
    update_text_axes(h)
    
% plot
    h = reset_axes(h);
    evalin('base','plot_basic_state_line')
    axes(h.lp.as.outline); cla(h.lp.as.outline,'reset')
    evalin('base','draw_contours')

set(h.lp.st.plotting,'Visible','off')
guidata(hObject,h)

function lp_et_sr_pos(hObject, ~, ~)
h = guidata(hObject);
set(h.lp.st.plotting,'Visible','on')

old_e_z = h.lp.e_z;

% determine position of h_fs
    h_fs = evalin('base','b1.R.v(end,:)');
    [r_c, r_i] = evalin('base','deal(r_c,r_i)');
    if length(h.blocks) == 2
        r_o = evalin('base','r_o');
    else
        r_o = r_i;
    end
    old_pos = r_c + h.lp.e_r(1)*(r_i-r_c) + h.lp.e_r(2)*(r_o-r_i);

% get input value
switch hObject
    case h.lp.et.pos
        if strcmp(h.lp.axis,'r')
            lowerLimit=0.5-length(h.blocks); upperLimit=-0.5+length(h.blocks); default=0;
        else
            lowerLimit=0; upperLimit=length(h.blocks); default=1;
        end
        h.lp.fsum = getValue(lowerLimit, upperLimit, default, hObject);
        h.lp.fsum = floor(h.lp.fsum*1000)/1000; % cut number off after 4 comma digits
        
    case h.lp.sr.pos
        h.lp.fsum = get(hObject,'Value');
        h.lp.fsum = round(h.lp.fsum,3); % resolution step = 0.001

end

% convert fake sum into real sum
if strcmp(h.lp.axis,'r')
    h.lp.sum = 1.5-h.lp.fsum;
else
    h.lp.sum = h.lp.fsum;
end

% update edit text box and slider and checkboxes
    set(h.lp.et.pos,'string',num2str(h.lp.fsum))
    set(h.lp.sr.pos,'Value',h.lp.fsum)
    if strcmp(h.lp.axis,'r')
        h.lp.e_z = [min([1,h.lp.sum]) max([0,h.lp.sum-1]) max([0,h.lp.sum-2])];
        if h.lp.e_z(2)>1, h.lp.e_z(2)=1; end
        if h.lp.sum < 1
            h.lp.domain = [0 1];
        elseif h.lp.sum <= 2
            if sum(old_e_z)<1 || sum(old_e_z)>2
                h.lp.domain = [1 1];
            end
        else
            h.lp.domain = [0 1];
        end
    else
        h.lp.e_r = [min([1,h.lp.sum]) max([0,h.lp.sum-1])];
        pos = r_c + h.lp.e_r(1)*(r_i-r_c) + h.lp.e_r(2)*(r_o-r_i);
        if strcmp(h.lp.coord,'curved') || prod(pos<=h_fs) == 1 || prod(pos>h_fs) == 1
            if h.lp.sum <= 1
                h.lp.domain = [1 0];
            else
                h.lp.domain = [0 1];
            end
        else % straight and LB+SG
            if (prod(old_pos<=h_fs) == 1 || prod(old_pos>h_fs) == 1) && max(abs(h_fs-r_i)/r_i) > 1e-12
                h.lp.domain = [1 1];
            end
        end
    end
    set(h.lp.cb.LB,'value',h.lp.domain(1)); set(h.lp.cb.SG,'value',h.lp.domain(2))

% update value
    assignin('base','tmp',h.lp.axis);
    evalin('base','lp.axis=tmp;');
    % -----------------------------------------------------------
    assignin('base','tmp',h.lp.domain);
    evalin('base','lp.domain=tmp;');
    if strcmp(h.lp.axis,'z')
        assignin('base','tmp',h.lp.e_r);
        evalin('base','lp.e_r=tmp;');
    else
        assignin('base','tmp',h.lp.e_z);
        evalin('base','lp.e_z=tmp;');
    end
    
% update text on axes
    update_text_axes(h)
    
% plot
    h = reset_axes(h);
    evalin('base','plot_basic_state_line')
    axes(h.lp.as.outline); cla(h.lp.as.outline,'reset')
    evalin('base','draw_contours')

set(h.lp.st.plotting,'Visible','off')
guidata(hObject,h)

function lp_et_sr_interval(hObject, ~, ~)
h = guidata(hObject);
set(h.lp.st.plotting,'Visible','on')

errormsg = 'Note: 0 \leq \alpha < \beta \leq 1 !';

switch hObject
    case h.lp.et.interval(1)
        % get input value
        lowerLimit=0; upperLimit=h.lp.zoom(2); default=0; showOption = h.show.interval;
        [h.lp.zoom(1), showMessage] = getValue(lowerLimit, upperLimit, default, hObject, errormsg, showOption, h.opts);
        
        % change value such that new value is slightly smaller
        if h.lp.zoom(1) == h.lp.zoom(2)
            h.lp.zoom(1) = h.lp.zoom(2) - 0.001;
        end
        
    case h.lp.sr.interval(1)
        h.lp.zoom(1) = get(hObject,'Value');
        h.lp.zoom(1) = round(h.lp.zoom(1),2);
        if h.lp.zoom(1) >= h.lp.zoom(2)
            h.lp.zoom(1) = h.lp.zoom(2)-0.001;
            h.lp.zoom(1) = floor(h.lp.zoom(1)*100)/100;
        end
        
    case h.lp.et.interval(2)
        % get input value
        lowerLimit=h.lp.zoom(1); upperLimit=1; default=1; showOption = h.show.interval;
        [h.lp.zoom(1), showMessage] = getValue(lowerLimit, upperLimit, default, hObject, errormsg, showOption, h.opts);
        
        % change value such that new value is slightly greater
        if h.lp.zoom(2) == h.lp.zoom(1)
            h.lp.zoom(2) = h.lp.zoom(1) + 0.001;
        end
        
    case h.lp.sr.interval(2)
        h.lp.zoom(2) = get(hObject,'Value');
        h.lp.zoom(2) = round(h.lp.zoom(2),2);
        if h.lp.zoom(2) <= h.lp.zoom(1)
            h.lp.zoom(2) = h.lp.zoom(1)+0.001;
            h.lp.zoom(2) = ceil(h.lp.zoom(2)*100)/100;
        end
        
end

if hObject == h.lp.et.interval(1) || hObject == h.lp.et.interval(2)
    % update show options
    if isempty(showMessage.error) == 0
        h.show.interval.error = showMessage.error;
    end
    if isempty(showMessage.warning) == 0
        h.show.interval.warning = showMessage.warning;
    end
end

% cut number off after 4 comma digits
    h.lp.zoom(1) = floor(h.lp.zoom(1)*1000)/1000;
    h.lp.zoom(2) = ceil(h.lp.zoom(2)*1000)/1000;

% update edit text box and slider
    set(h.lp.et.interval(1),'string',num2str(h.lp.zoom(1)))
    set(h.lp.sr.interval(1),'Value',h.lp.zoom(1))
    set(h.lp.et.interval(2),'string',num2str(h.lp.zoom(2)));
    set(h.lp.sr.interval(2),'Value',h.lp.zoom(2))
   
% update value and plot
    h = reset_axes(h);
    assignin('base','tmp',h.lp.zoom);
    evalin('base','lp.zoom=tmp;');
    evalin('base','plot_basic_state_line')

set(h.lp.st.plotting,'Visible','off')
guidata(hObject,h)