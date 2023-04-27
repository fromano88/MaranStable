function plot_pF

f_plot = figure('units','pixels','menubar','none','name','Post Processing - Perturbation Flow',...
    'numbertitle','off','ToolBar','figure','resize','off','Visible','Off');

if ispc % windows
    f_size = [1264 650];
else    % linux
    f_size = [1296 760];
end
f_plot.Position(3:4) = f_size;
movegui(f_plot,'center')

initialPlot_perturbation

% create figure
% windows or linux
    % w=1 --> windows
    % w=0 --> linux
    w = ispc;

% save writings
    R = {'style','rad','units','pix','parent'};
    C = {'style','check','units','pix','parent'};
    E = {'style','edit','units','pix','parent'};
    A = {'unit','pix','visible','off','parent'};
    S = {'style','text','units','pix','HorizontalAlignment','left','parent'};

% main panels: ct = contourplot
    h.pl.plot = uipanel('unit','pix','pos',[0 0 988-28*w 760-110*w],'BorderType',...
        'none','BackgroundColor',[1 1 1]);
    
% subpanels
    h.ct.pl.position(1) = uipanel('units','pix','pos',[995-27*w 607-89*w 296-5*w 52-4*w],'title','Position');
    h.ct.pl.position(2) = uipanel('units','pix','pos',h.ct.pl.position(1).Position,'title','Position','Visible','off');
    h.ct.pl.position(3) = uipanel('units','pix','pos',h.ct.pl.position(1).Position,'title','Position','Visible','off');
    h.ct.pl.domain(1)   = uipanel('units','pix','pos',[995-27*w 517-79*w 296-5*w 86-10*w],'title','Domain');
    h.ct.pl.quantity    = uipanel('units','pix','pos',[995-27*w 403-47*w 296-5*w 110-32*w],'title','Contour Plot');
    h.ct.pl.vector      = uipanel('units','pix','pos',[995-27*w 238-32*w 296-5*w 161-15*w],'title','Velocity Vectors');
    h.ct.pl.lines(1)    = uipanel('units','pix','pos',[995-27*w 41-2*w 296-5*w 193-30*w],'title','Isolines - Basic State');
    
% subsubpanels
    h.ct.pl.domain(2) = uipanel(h.ct.pl.domain(1),'units','pix','pos',[170 -5 131 96]);
    h.ct.pl.lines(2) = uipanel(h.ct.pl.lines(1),'units','pix','pos',[0 -1 310 97.5-11*w],'title','Isolines - Perturbation Flow');
    
% button groups
    h.bg.cut =       uibuttongroup('units','pix','pos',[995-27*w 703-97*w 296-5*w 50-7*w],'title','2D space','SelectionChangedFcn',{@bg_cut,h});
    h.bg.dimension = uibuttongroup('units','pix','pos',[995-27*w 663-93*w 296-5*w 34-4*w],'SelectionChangedFcn',{@bg_dimension,h});
    
% radio buttons
    h.rb.dimless = uicontrol(R{:},h.bg.dimension,'pos',[12 5-1*w 116 22],'string','Dimensionless','value',1);
    h.rb.dim = uicontrol(R{:},h.bg.dimension,'pos',[177 5-1*w 101 22],'string','Dimensional');
    h.rb.r_z = uicontrol('style','rad','unit','pix','pos',[12 8-2*w 15 20],'string',' ','value',1,'parent',h.bg.cut);
    h.rb.r_phi = uicontrol('style','rad','unit','pix','pos',[12+100 8-2*w 15 20],'string','  ','parent',h.bg.cut);
    h.rb.z_phi = uicontrol('style','rad','unit','pix','pos',[12+200 8-2*w 15 20],'string','   ','parent',h.bg.cut);
    
% check boxes
    h.ct.cb.LB = uicontrol(C{:},h.ct.pl.domain(1),'pos',[20 43-5*w 130 20],'String','Liquid Phase','value',1,'call',{@ct_cb_domain,h});
    h.ct.cb.SG = uicontrol(C{:},h.ct.pl.domain(1),'pos',[20 12-3*w 130 20],'String','Gas Phase','value',1,'call',{@ct_cb_domain,h});
    h.ct.cb.equal = uicontrol(C{:},h.ct.pl.domain(2),'pos',[20 31-4*w 130 20],'String','Equal Axis','value',1,'call',{@ct_cb_equal,h});
    h.ct.cb.vector = uicontrol(C{:},h.ct.pl.vector,'pos',[20 118-10*w 50 20],'String','Off','value',0,'call',{@ct_cb_vector,h});
    h.ct.cb.iso_sf_bS = uicontrol(C{:},h.ct.pl.lines(1),'pos',[8 141-22*w 100 20],'String','Streamlines','value',h.sf_bS,'call',{@ct_cb_isolines,h});
    h.ct.cb.iso_temp_bS = uicontrol(C{:},h.ct.pl.lines(1),'pos',[8 107-19*w 160 20],'String','Temperature','value',h.temp_bS,'call',{@ct_cb_isolines,h});
    h.ct.cb.iso_prod_pF = uicontrol(C{:},h.ct.pl.lines(2),'pos',[8 48-5*w 160 20],'String','Th. Production','value',h.prod_pF,'call',{@ct_cb_isolines,h});
    h.ct.cb.iso_temp_pF = uicontrol(C{:},h.ct.pl.lines(2),'pos',[8 12-2*w 160 20],'String','Temperature','value',h.temp_pF,'call',{@ct_cb_isolines,h});
    if evalin('base','flowopt.energy') == 0
        set([h.ct.cb.iso_temp_bS h.ct.cb.iso_temp_pF h.ct.cb.iso_prod_pF],'enable','off')
    end
    
% edit texts
    h.ct.et.position(1) = uicontrol(E{:},h.ct.pl.position(1),'pos',[60 8 80 23],'string',num2str(h.ct.phi),'call',{@ct_et_position,h});
    h.ct.et.position(2) = uicontrol(E{:},h.ct.pl.position(2),'pos',[60 8 80 23],'string',num2str(h.ct.z),'call',{@ct_et_position,h});
    h.ct.et.position(3) = uicontrol(E{:},h.ct.pl.position(3),'pos',[60 8 80 23],'string',num2str(h.ct.r),'call',{@ct_et_position,h});
    h.ct.et.vector(1) = uicontrol(E{:},h.ct.pl.vector,'pos',[195-20*w 86-5*w 60 23],'string',num2str(h.vp.numR),'enable','off','call',{@ct_et_arrows,h});
    h.ct.et.vector(2) = uicontrol(E{:},h.ct.pl.vector,'pos',[195-20*w 55-5*w 60 23],'string',num2str(h.vp.numZ),'enable','off','call',{@ct_et_arrows,h});
    h.ct.et.iso_sf_bS = uicontrol(E{:},h.ct.pl.lines(1),'pos',[122 141-22*w 60 23],'string',num2str(h.iso_sf_bS),'enable','off','call',{@ct_et_isolines,h});
    h.ct.et.iso_temp_bS = uicontrol(E{:},h.ct.pl.lines(1),'pos',[122 105-19*w 60 23],'string',num2str(h.iso_temp_bS),'enable','off','call',{@ct_et_isolines,h});
    h.ct.et.iso_prod_pF = uicontrol(E{:},h.ct.pl.lines(2),'pos',[122 46-5*w 60 23],'string',num2str(h.iso_prod_pF),'enable','off','call',{@ct_et_isolines,h});
    h.ct.et.iso_temp_pF = uicontrol(E{:},h.ct.pl.lines(2),'pos',[122 10-2*w 60 23],'string',num2str(h.iso_temp_pF),'enable','off','call',{@ct_et_isolines,h});

% popup menus
    if evalin('base','flowopt.energy') == 1
        h.ct.pm.quantity = uicontrol(h.ct.pl.quantity,'style','pop','pos',[14 55-20*w 150 24],'string',{'Temperature','Thermal Production','u - Velocity','v - Velocity','w - Velocity','Pressure','None'},'call',{@ct_pm_quantity,h});
    else
        h.ct.pm.quantity = uicontrol(h.ct.pl.quantity,'style','pop','pos',[14 55-20*w 150 24],'string',{'u - Velocity','v - Velocity','w - Velocity','Pressure','None'},'call',{@ct_pm_quantity,h});
    end
    h.ct.pm.colormap = uicontrol(h.ct.pl.quantity,'style','pop','pos',[94 14-10*w 115 24],'string',{'Cold2Warm','Cool2Warm','Rainbow','Parula','Thermal','Inferno','Black&White'},'call',{@ct_pm_colormap,h});
    h.ct.pm.color(1) = uicontrol(h.ct.pl.vector,'style','pop','pos',[70 15-5*w 84 24],'string',{'Black','Blue','Green','Red','White'},'enable','off','call',{@ct_pm_color,h});
    h.ct.pm.color(2) = uicontrol(h.ct.pl.lines(1),'style','pop','pos',[196 141-24*w 84 23],'string',{'Black','Blue','Green','Red','White'},'enable','off','call',{@ct_pm_color,h});
    h.ct.pm.color(3) = uicontrol(h.ct.pl.lines(1),'style','pop','pos',[196 105-21*w 84 23],'string',{'Black','Blue','Green','Red','White'},'enable','off','call',{@ct_pm_color,h});
    h.ct.pm.color(4) = uicontrol(h.ct.pl.lines(2),'style','pop','pos',[196 46-7*w 84 23],'string',{'Black','Blue','Green','Red','White'},'enable','off','call',{@ct_pm_color,h});
    h.ct.pm.color(5) = uicontrol(h.ct.pl.lines(2),'style','pop','pos',[196 10-4*w 84 23],'string',{'Black','Blue','Green','Red','White'},'enable','off','call',{@ct_pm_color,h});

% pushbuttons
    h.pb.export = uicontrol('units','pix','visible','on','string','Export',...
        'pos',[1214-32*w 8 76 25],'ForegroundColor',[0.8500 0.3250 0.0980],...
        'call',{@pb_export,h});
        
% axes
    h.latex_l = {'interpreter','Latex','FontSize',12,'HorizontalAlignment','right'};
    h.latex_r = {'interpreter','Latex','FontSize',12,'HorizontalAlignment','left'};
    h.as.plot = axes(A{:},h.pl.plot,'OuterPosition',[0 0 988 760-110*w]);
    if evalin('base','flowopt.ax') == 1
        h.ct.as.position(1) = axes(A{:},h.ct.pl.position(1),'pos',[60 8 80 23]); text(-0.09,0.47,'$\varphi=$',h.latex_l{:}); text(1.05,0.47,'$\cdot\, \pi$',h.latex_r{:});
        h.ct.as.position(2) = axes(A{:},h.ct.pl.position(2),'pos',[60 8 80 23]); text(-0.09,0.47,'$z=$',h.latex_l{:});
        h.ct.as.position(3) = axes(A{:},h.ct.pl.position(3),'pos',[60 8 80 23]); text(-0.09,0.47,'$r=$',h.latex_l{:}); text(1.034,0.47,'$\cdot h(z)$',h.latex_r{:});
        h.ct.as.r_z = axes(A{:},h.bg.cut,'pos',[32 17-2*w 85 20]); text(0,0,'$r$ - $z$',h.latex_r{:})
        h.ct.as.r_phi = axes(A{:},h.bg.cut,'pos',[32+100 17-2*w 85 20]); text(0,0,'$r$ - $\varphi$',h.latex_r{:})
        h.ct.as.z_phi = axes(A{:},h.bg.cut,'pos',[32+100*2 17-2*w 85 20]); text(0,0,'$z$ - $\varphi$',h.latex_r{:})
        h.ct.as.vec(1) = axes(A{:},h.ct.pl.vector,'pos',h.ct.et.vector(1).Position); text(-0.63-0.3*w,0.47,'$r$:',h.latex_r{:});
        h.ct.as.vec(2) = axes(A{:},h.ct.pl.vector,'pos',h.ct.et.vector(2).Position); text(-0.63-0.3*w,0.47,'$z$:',h.latex_r{:});
    else
        h.ct.as.position(1) = axes(A{:},h.ct.pl.position(1),'pos',[60 8 80 23]); text(-0.09,0.47,'$z=$',h.latex_l{:});
        h.ct.as.position(2) = axes(A{:},h.ct.pl.position(2),'pos',[60 8 80 23]); text(-0.09,0.47,'$y=$',h.latex_l{:});
        h.ct.as.position(3) = axes(A{:},h.ct.pl.position(3),'pos',[60 8 80 23]); text(-0.09,0.47,'$x=$',h.latex_l{:}); text(1.034,0.47,'$\cdot h(y)$',h.latex_r{:});
        h.ct.as.r_z = axes(A{:},h.bg.cut,'pos',[32 17-2*w 85 20]); text(0,0,'$x$ - $y$',h.latex_r{:})
        h.ct.as.r_phi = axes(A{:},h.bg.cut,'pos',[32+100 17-2*w 85 20]); text(0,0,'$x$ - $z$',h.latex_r{:})
        h.ct.as.z_phi = axes(A{:},h.bg.cut,'pos',[32+100*2 17-2*w 85 20]); text(0,0,'$y$ - $z$',h.latex_r{:})
        h.ct.as.vec(1) = axes(A{:},h.ct.pl.vector,'pos',h.ct.et.vector(1).Position); text(-0.63-0.3*w,0.47,'$x$:',h.latex_r{:});
        h.ct.as.vec(2) = axes(A{:},h.ct.pl.vector,'pos',h.ct.et.vector(2).Position); text(-0.63-0.3*w,0.47,'$y$:',h.latex_r{:});
    end
    text_axes = [h.ct.as.position h.ct.as.r_z h.ct.as.r_phi h.ct.as.z_phi h.ct.as.vec];
    for i = text_axes
        i.Toolbar.Visible = 'off';
    end

% static texts
    uicontrol(S{:},h.ct.pl.quantity,'pos',[14 17-9*w 65 17],'string','Colormap');
    uicontrol(S{:},h.ct.pl.vector,'pos',[20 88-6*w 135-37*w 17],'string','Number of Arrows in');
    h.ct.st.vec = uicontrol(S{:},h.ct.pl.vector,'pos',[20 57-6*w 135-37*w 17],'string','Number of Arrows in');
    uicontrol(S{:},h.ct.pl.vector,'pos',[20 18-5*w 40 17],'string','Color');
    
% SFM or TFM
    h.blocks = evalin('base','blocks');
    h.l_lb = evalin('base','l_lb');
    if length(h.blocks) == 1
        h.l_d1 = 0;
        h.l_d2 = 0;
    else
        h.l_d1 = evalin('base','l_d1');
        h.l_d2 = evalin('base','l_d2');
    end
    
% only difference if SFM or TFM
    if length(h.blocks) == 1
        set(h.ct.cb.LB,'enable','inactive')
        set(h.ct.cb.SG,'value',0,'enable','off')
    end
    
% plot default setup
    h.cut = 'r-z';
    axes(h.as.plot);
    evalin('base','perturbations_r_z')
    set(h.as.plot,'visible','on')

% position values of the plot to fix the plot at one position
    h.position = get(h.as.plot,'Position');  
    
% showing errors or not
    h.show.position.warning = 1;
    h.show.arrow.warning = 1;
    h.show.lines.warning = 1;
    
% Windows style
    h.opts = struct('WindowStyle','non-modal','Interpreter','tex');

% correct figure size
    pause(0.1)
    f_plot.Position(3:4) = f_size;
    set(f_plot,'Visible','On')
    
guidata(f_plot,h)

function pb_export(~, ~, ~)
    export_vtk('pF')

function bg_cut(hObject, ~, ~)
h = guidata(hObject);

h.cut = get(get(hObject,'SelectedObject'),'String');
switch h.cut
    case ' '
        h.cut = 'r-z';
    case '  '
        h.cut = 'r-phi';
    case '   '
        h.cut = 'z-phi';
end

% domain for TFM
if length(h.blocks) == 2
    switch h.cut
        case 'r-z'
            set([h.ct.cb.LB h.ct.cb.SG],'enable','on');
    
        case 'r-phi'
            set(h.ct.cb.SG,'enable','on');
            if h.ct.z_shift < h.l_d1 || h.ct.z_shift > h.l_d1+h.l_lb
                h.ct.domain = [0 1];
                set(h.ct.cb.LB,'enable','off');
            end
                
        case 'z-phi'
            set(h.ct.cb.LB,'enable','on')
            set(h.ct.cb.SG,'enable','off')
            h.ct.domain = [1 0];
    end
    set(h.ct.cb.LB,'value',h.ct.domain(1)); set(h.ct.cb.SG,'value',h.ct.domain(2))
end

% update value
    assignin('base','tmp',h.ct.domain);
    evalin('base','ct.domain=tmp;');

% update plot
    w = ispc;
    axes(h.as.plot); cla(h.as.plot,'reset');
    switch h.cut
        case 'r-z'
            evalin('base','perturbations_r_z')
            
            set([h.ct.pl.position(1) h.ct.st.vec h.ct.et.vector(2)],'Visible','on')
            set(h.ct.pl.position(2:3),'Visible','off')
            set([h.ct.cb.iso_sf_bS h.ct.cb.iso_temp_bS],'enable','on')
            if h.sf_bS == 1
                set([h.ct.et.iso_sf_bS h.ct.pm.color(2)],'enable','on')
                set(h.ct.cb.iso_sf_bS,'value',1)
            end
            if h.temp_bS == 1
                set([h.ct.et.iso_temp_bS h.ct.pm.color(3)],'enable','on')
                set(h.ct.cb.iso_temp_bS,'value',1)
            end

            axes(h.ct.as.vec(1)); cla(h.ct.as.vec(1));
            if evalin('base','flowopt.ax') == 1
                text(-0.63-0.3*w,0.47,'$r$:',h.latex_r{:});
            else
                text(-0.63-0.3*w,0.47,'$x$:',h.latex_r{:});
            end
            % ----------------------------------------------
            axes(h.ct.as.vec(2)); cla(h.ct.as.vec(2));
            if evalin('base','flowopt.ax') == 1
                text(-0.63-0.3*w,0.47,'$z$:',h.latex_r{:});
            else
                text(-0.63-0.3*w,0.47,'$y$:',h.latex_r{:});
            end

        case 'r-phi'
            evalin('base','perturbations_r_phi')
            
            set(h.ct.pl.position(2),'Visible','on')
            set(h.ct.pl.position([1 3]),'Visible','off')
            set([h.ct.cb.iso_sf_bS h.ct.cb.iso_temp_bS],'value',0,'enable','off')
            set([h.ct.et.iso_sf_bS h.ct.et.iso_temp_bS h.ct.pm.color(2:3)],'enable','off')
            
            axes(h.ct.as.position(2)); cla(h.ct.as.position(2));
            if evalin('base','flowopt.ax') == 1
                text(-0.09,0.47,'$z=$',h.latex_l{:});
            else
                text(-0.09,0.47,'$y=$',h.latex_l{:});
            end
            if h.dim == 1
                text(1.05,0.47,'mm',h.latex_r{:});
            end
            % ----------------------------------------------
            axes(h.ct.as.vec(1)); cla(h.ct.as.vec(1));
            if evalin('base','flowopt.ax') == 1
                text(-0.63-0.3*w,0.47,'$r$:',h.latex_r{:});
            else
                text(-0.63-0.3*w,0.47,'$x$:',h.latex_r{:});
            end
            % ----------------------------------------------
            axes(h.ct.as.vec(2)); cla(h.ct.as.vec(2));
            if evalin('base','flowopt.ax') == 1
                set([h.ct.st.vec h.ct.et.vector(2) h.ct.as.vec(2)],'Visible','off')
            else
                text(-0.63-0.3*w,0.47,'$z$:',h.latex_r{:});
            end

        case 'z-phi'
            evalin('base','perturbations_z_phi')
            
            set([h.ct.pl.position(3) h.ct.st.vec h.ct.et.vector(2)],'Visible','on')
            set(h.ct.pl.position(1:2),'Visible','off')
            set([h.ct.cb.iso_sf_bS h.ct.cb.iso_temp_bS],'value',0,'enable','off')
            set([h.ct.et.iso_sf_bS h.ct.et.iso_temp_bS h.ct.pm.color(2:3)],'enable','off')

            axes(h.ct.as.vec(1)); cla(h.ct.as.vec(1));
            if evalin('base','flowopt.ax') == 1
                text(-0.63-0.3*w,0.47,'$\varphi$:',h.latex_r{:});
            else
                text(-0.63-0.3*w,0.47,'$z$:',h.latex_r{:});
            end
            % ----------------------------------------------
            axes(h.ct.as.vec(2)); cla(h.ct.as.vec(2));
            if evalin('base','flowopt.ax') == 1
                text(-0.63-0.3*w,0.47,'$z$:',h.latex_r{:});
            else
                text(-0.63-0.3*w,0.47,'$y$:',h.latex_r{:});
            end
    end
    set(h.as.plot,'Position',h.position)

    if evalin('base','flowopt.energy') == 0
        set([h.ct.cb.iso_temp_bS h.ct.cb.iso_temp_pF h.ct.cb.iso_prod_pF h.ct.pm.color(3:5)],'enable','off')
    end
    
    
% update number of shown arrows in terms of new domain
    if get(h.ct.cb.vector,'value')
        switch h.cut
            case 'r-z'
                if get(h.ct.cb.SG,'value')==1 && get(h.ct.cb.LB,'value')==0
                    h.vp.n_r = evalin('base','round((b2.geom.r(end) - b2.geom.r(1))/(max(geom.r)-min(geom.r))*vp.numR)');
                    set(h.ct.et.vector(2),'string',num2str(h.vp.numZ))
                    set(h.ct.et.vector(1),'string',num2str(h.vp.n_r))
                elseif get(h.ct.cb.LB,'value')==1 && get(h.ct.cb.SG,'value')==0
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
            case 'r-phi'
                set(h.ct.et.vector(1),'string',num2str(h.vp.numR))
                set(h.ct.et.vector(2),'string',num2str(h.vp.numPhi))
            case 'z-phi'
                set(h.ct.et.vector(1),'string',num2str(h.vp.numPhi))
                set(h.ct.et.vector(2),'string',num2str(h.vp.numZ))
        end
        
    end
    
guidata(hObject,h)

function bg_dimension(hObject, ~, ~)
h = guidata(hObject);

selection = get(get(hObject,'SelectedObject'),'String');

switch selection
    case 'Dimensionless'
        h.dim = 0; h.f = 1/h.l_lb;
        h.ct.z = h.ct.z/h.l_lb/1000;
        set(h.ct.et.position(2),'string',num2str(h.ct.z))
 
    case 'Dimensional'
        h.dim = 1; h.f = 1000;
        h.ct.z = h.ct.z*h.l_lb*1000;
        set(h.ct.et.position(2),'string',num2str(h.ct.z))
end

axes(h.ct.as.position(2)); cla(h.ct.as.position(2));
text(-0.09,0.47,'$z=$',h.latex_l{:});
if h.dim == 1
    text(1.05,0.47,'mm',h.latex_r{:});
end

% update plot
    assignin('base','tmp',h.dim);
    evalin('base','dim=tmp;');
    % ------------------------------------------------------------
    h.ct.z_shift = h.l_lb/2+h.l_d1-h.ct.z/h.f;
    assignin('base','tmp',h.ct.z_shift);
    evalin('base','ct.z=tmp;');
    % ------------------------------------------------------------
    axes(h.as.plot); cla(h.as.plot,'reset');
    switch h.cut
        case 'r-z'
            evalin('base','perturbations_r_z') 
        case 'r-phi'
            evalin('base','perturbations_r_phi')
        case 'z-phi'
            evalin('base','perturbations_z_phi')
    end
    set(h.as.plot,'Position',h.position)
    
guidata(hObject,h)

function ct_et_position(hObject, ~, ~)
h = guidata(hObject);

h.l_lb = evalin('base','l_lb');

if h.dim == 0
    h.f = 1/h.l_lb;
else
    h.f = 1000; %mm
end

switch hObject
    case h.ct.et.position(1)
        lowerLimit = -inf; upperLimit = inf; default = 0;
        [h.ct.phi,~] = getValue(lowerLimit,upperLimit,default,hObject);
        assignin('base','tmp',h.ct.phi);
        evalin('base','ct.phi=tmp;');
        
    case h.ct.et.position(2)
        lowerLimit = -(h.l_lb/2+h.l_d2)*h.f; upperLimit = (h.l_lb/2+h.l_d1)*h.f; default = h.z_T_max*h.f;
        [h.ct.z,~] = getValue(lowerLimit,upperLimit,default,hObject);
        h.ct.z_shift = h.l_d1+h.l_lb/2-h.ct.z/h.f;
        if length(h.blocks) == 2
            if h.ct.z_shift < h.l_d1 || h.ct.z_shift > h.l_d1+h.l_lb
                h.ct.domain = [0 1];
                set(h.ct.cb.LB,'enable','off');
            else
                if evalin('base','ct.z') < h.l_d1 || evalin('base','ct.z') > h.l_d1+h.l_lb
                    h.ct.domain = [1 1];
                end
                set(h.ct.cb.LB,'enable','on');
            end
            set(h.ct.cb.LB,'value',h.ct.domain(1)); set(h.ct.cb.SG,'value',h.ct.domain(2))
        end

        assignin('base','tmp',h.ct.z_shift);
        evalin('base','ct.z=tmp;');
        % --------------------------------------------------------- %
        assignin('base','tmp',h.ct.domain);
        evalin('base','ct.domain=tmp;');

    case h.ct.et.position(3)
        lowerLimit = 0; upperLimit = 1; default = 1;
        [h.ct.r,~] = getValue(lowerLimit,upperLimit,default,hObject);
        assignin('base','tmp',h.ct.r);
        evalin('base','ct.r=tmp;');
         
end

% update plot
    axes(h.as.plot); cla(h.as.plot,'reset');
    switch h.cut
        case 'r-z'
            evalin('base','perturbations_r_z') 
        case 'r-phi'
            evalin('base','perturbations_r_phi')
        case 'z-phi'
            evalin('base','perturbations_z_phi')
    end
    set(h.as.plot,'Position',h.position)

    
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
    axes(h.as.plot); cla(h.as.plot,'reset');
    switch h.cut
        case 'r-z'
            evalin('base','perturbations_r_z') 
        case 'r-phi'
            evalin('base','perturbations_r_phi')
        case 'z-phi'
            evalin('base','perturbations_z_phi')
    end
    set(h.as.plot,'Position',h.position)

% update number of shown arrows in terms of new domain
    if get(h.ct.cb.vector,'value')
        if strcmp(h.ct.domain,'SG')
            h.vp.n_r = evalin('base','round((b2.geom.r(end) - b2.geom.r(1))/(max(geom.r)-min(geom.r))*vp.numR)');
            set(h.ct.et.vector(1),'string',num2str(h.vp.n_r))
            set(h.ct.et.vector(2),'string',num2str(h.vp.numZ))
        elseif strcmp(h.ct.domain,'LB')
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

function ct_cb_equal(hObject, ~, ~)
h = guidata(hObject);

h.ct.equal = get(hObject,'value');

% update value
    assignin('base','tmp',h.ct.equal);
    evalin('base','ct.equal=tmp;');

% update plot
    axes(h.as.plot); cla(h.as.plot,'reset');
    switch h.cut
        case 'r-z'
            evalin('base','perturbations_r_z') 
        case 'r-phi'
            evalin('base','perturbations_r_phi')
        case 'z-phi'
            evalin('base','perturbations_z_phi')
    end
    set(h.as.plot,'Position',h.position)

guidata(hObject,h)

function ct_pm_quantity(hObject, ~, ~)
h = guidata(hObject);

str=get(hObject, 'String'); val=get(hObject, 'Value');

h.ct.draw = 1;
set(h.ct.pm.colormap,'enable','on')

switch str{val}
    case 'Temperature'
        h.ct.quantity = 'T';
    
    case 'Thermal Production'
        h.ct.quantity = 'th_E';

    case 'u - Velocity'
        h.ct.quantity = 'u';

    case 'v - Velocity'
        if evalin('base','flowopt.ax') == 1
            h.ct.quantity = 'v';
        else
            h.ct.quantity = 'w';
        end
        
    case 'w - Velocity'
        if evalin('base','flowopt.ax') == 1
            h.ct.quantity = 'w';
        else
            h.ct.quantity = 'v';
        end
        
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
    axes(h.as.plot); cla(h.as.plot,'reset');
    switch h.cut
        case 'r-z'
            evalin('base','perturbations_r_z') 
        case 'r-phi'
            evalin('base','perturbations_r_phi')
        case 'z-phi'
            evalin('base','perturbations_z_phi')
    end
    set(h.as.plot,'Position',h.position)

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
    axes(h.as.plot); cla(h.as.plot,'reset');
    switch h.cut
        case 'r-z'
            evalin('base','perturbations_r_z') 
        case 'r-phi'
            evalin('base','perturbations_r_phi')
        case 'z-phi'
            evalin('base','perturbations_z_phi')
    end
    set(h.as.plot,'Position',h.position)

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
    axes(h.as.plot); cla(h.as.plot,'reset');
    switch h.cut
        case 'r-z'
            evalin('base','perturbations_r_z') 
        case 'r-phi'
            evalin('base','perturbations_r_phi')
        case 'z-phi'
            evalin('base','perturbations_z_phi')
    end
    set(h.as.plot,'Position',h.position)
    
% update number of shown arrows in terms of new domain
    if get(h.ct.cb.vector,'value')
        if strcmp(h.ct.domain,'SG')
            h.vp.n_r = evalin('base','round((b2.geom.r(end) - b2.geom.r(1))/(max(geom.r)-min(geom.r))*vp.numR)');
            set(h.ct.et.vector(2),'string',num2str(h.vp.numZ))
            set(h.ct.et.vector(1),'string',num2str(h.vp.n_r))
        elseif strcmp(h.ct.domain,'LB')
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
    if strcmp(h.cut,'z-phi')
        set(h.ct.et.vector(1),'string',num2str(h.vp.numPhi))
    end

guidata(hObject,h)

function ct_et_arrows(hObject, ~, ~)
h = guidata(hObject);

lowerLimit=0; upperLimit=inf; default=h.vp.numR;

switch hObject
    case h.ct.et.vector(1)
        default=h.vp.numR;
end

[numArrow,~] = getValue(lowerLimit,upperLimit,default,hObject);
    
if floor(numArrow) ~= numArrow || numArrow <= 0
    numArrow = floor(norm(numArrow));
    set(hObject,'string',num2str(numArrow))
end

LB = get(h.ct.cb.LB,'value'); SG = get(h.ct.cb.SG,'value');

switch hObject
    case h.ct.et.vector(1)
        switch h.cut
            case {'r-phi','r-z'}
                h.vp.numR = numArrow;
                % compute number of arrows with respect to whole geometry
                if LB == 1 && SG == 0
                    h.vp.numR = evalin('base','round((max(geom.r)-min(geom.r))/(b1.geom.r(end) - b1.geom.r(1)))')*h.vp.numR;
                elseif LB == 0 && SG == 1
                    h.vp.numR = evalin('base','round((max(geom.r)-min(geom.r))/(b2.geom.r(end) - b2.geom.r(1)))')*h.vp.numR;
                end
                
            case 'z-phi'
                h.vp.numPhi = numArrow;
                
        end
        % update value
        assignin('base','tmp',h.vp.numR);
        evalin('base','vp.numR=tmp;');
        % ------------------------------------------------------
        assignin('base','tmp',h.vp.numPhi);
        evalin('base','vp.numPhi=tmp;');
        
    case h.ct.et.vector(2)
        switch h.cut
            case {'r-z','z-phi'}
                h.vp.numZ = numArrow;
                % compute number of arrows with respect to whole geometry
                if LB == 1 && SG == 0
                    h.vp.numZ = evalin('base','round((max(geom.z)-min(geom.z))/(b1.geom.z(end) - b1.geom.z(1)))')*(h.vp.numZ+1);
                end
                % update value
                assignin('base','tmp',h.vp.numZ-1);
                evalin('base','vp.numZ=tmp;');

            case 'r-phi'
                h.vp.numPhi = numArrow;
                % update value
                assignin('base','tmp',h.vp.numPhi);
                evalin('base','vp.numPhi=tmp;');
        end
end
    
% update plot
    axes(h.as.plot); cla(h.as.plot,'reset');
    switch h.cut
        case 'r-z'
            evalin('base','perturbations_r_z') 
        case 'r-phi'
            evalin('base','perturbations_r_phi')
        case 'z-phi'
            evalin('base','perturbations_z_phi')
    end
    set(h.as.plot,'Position',h.position)

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
            evalin('base','ct.color.sf_bS=eval(tmp);');
        case h.ct.pm.color(3)
            evalin('base','ct.color.T_bS=eval(tmp);');
        case h.ct.pm.color(4)
            evalin('base','ct.color.prod_pF=eval(tmp);');
        case h.ct.pm.color(5)
            evalin('base','ct.color.T_pF=eval(tmp);');
    end

% update plot
    axes(h.as.plot); cla(h.as.plot,'reset');
    switch h.cut
        case 'r-z'
            evalin('base','perturbations_r_z') 
        case 'r-phi'
            evalin('base','perturbations_r_phi')
        case 'z-phi'
            evalin('base','perturbations_z_phi')
    end
    set(h.as.plot,'Position',h.position)

guidata(hObject,h)

function ct_cb_isolines(hObject, ~, ~)
h = guidata(hObject);

switch hObject
    case h.ct.cb.iso_sf_bS
        h.sf_bS = get(hObject,'Value');

        if get(hObject,'Value')
            set([h.ct.et.iso_sf_bS h.ct.pm.color(2)],'enable','on')
        else
            set([h.ct.et.iso_sf_bS h.ct.pm.color(2)],'enable','off')
        end
        % update value
        assignin('base','tmp',h.sf_bS);
        evalin('base','ct.sf_bS=tmp;');
    
    case h.ct.cb.iso_temp_bS
        h.temp_bS = get(hObject,'Value');

        if get(hObject,'Value')
            set([h.ct.et.iso_temp_bS h.ct.pm.color(3)],'enable','on')
        else
            set([h.ct.et.iso_temp_bS h.ct.pm.color(3)],'enable','off')
        end
        % update value
        assignin('base','tmp',h.temp_bS);
        evalin('base','ct.temp_bS=tmp;');
        
    case h.ct.cb.iso_prod_pF
        h.prod_pF = get(hObject,'Value');

        if get(hObject,'Value')
            set([h.ct.et.iso_prod_pF h.ct.pm.color(4)],'enable','on')
        else
            set([h.ct.et.iso_prod_pF h.ct.pm.color(4)],'enable','off')
        end
        % update value
        assignin('base','tmp',h.prod_pF);
        evalin('base','ct.prod_pF=tmp;');
        
    case h.ct.cb.iso_temp_pF
        h.temp_pF = get(hObject,'Value');

        if get(hObject,'Value')
            set([h.ct.et.iso_temp_pF h.ct.pm.color(5)],'enable','on')
        else
            set([h.ct.et.iso_temp_pF h.ct.pm.color(5)],'enable','off')
        end
        % update value
        assignin('base','tmp',h.temp_pF);
        evalin('base','ct.temp_pF=tmp;');
end
     
% update plot
    axes(h.as.plot); cla(h.as.plot,'reset');
    switch h.cut
        case 'r-z'
            evalin('base','perturbations_r_z') 
        case 'r-phi'
            evalin('base','perturbations_r_phi')
        case 'z-phi'
            evalin('base','perturbations_z_phi')
    end
    set(h.as.plot,'Position',h.position)

guidata(hObject,h)

function ct_et_isolines(hObject, ~, ~)
h = guidata(hObject);

lowerLimit=-inf; upperLimit=inf; default=h.vp.numR;

[numLines,~] = getValue(lowerLimit,upperLimit,default,hObject);
    
if floor(numLines) ~= numLines || numLines <= 0
    numLines = floor(norm(numLines));
    set(hObject,'string',num2str(numLines))
end

switch hObject
    case h.ct.et.iso_sf_bS
        h.iso_sf_bS = numLines;
        % update value
        assignin('base','tmp',h.iso_sf_bS);
        evalin('base','ct.iso_sf_bS=tmp;');
        
    case h.ct.et.iso_temp_bS
        h.iso_temp_bS = numLines;
        % update value
        assignin('base','tmp',h.iso_temp_bS);
        evalin('base','ct.iso_temp_bS=tmp;');
        
    
    case h.ct.et.iso_prod_pF
        h.iso_prod_pF = numLines;
        % update value
        assignin('base','tmp',h.iso_prod_pF);
        evalin('base','ct.iso_prod_pF=tmp;');
        
    
    case h.ct.et.iso_temp_pF
        h.iso_temp_pF = numLines;
        % update value
        assignin('base','tmp',h.iso_temp_pF);
        evalin('base','ct.iso_temp_pF=tmp;');
end

% update plot
    axes(h.as.plot); cla(h.as.plot,'reset');
    switch h.cut
        case 'r-z'
            evalin('base','perturbations_r_z') 
        case 'r-phi'
            evalin('base','perturbations_r_phi')
        case 'z-phi'
            evalin('base','perturbations_z_phi')
    end
    set(h.as.plot,'Position',h.position)

guidata(hObject,h)