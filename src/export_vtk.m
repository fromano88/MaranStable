function export_vtk(flow)

% flow is either bS (basic state) or pF (perturbation flow)
    if strcmp(flow,'bS')
        s = 0;
        h.export.bS = 1; h.export.pF = 0;
    else
        s = 25;
        h.export.bS = 0; h.export.pF = 1;
    end

f_export = figure('units','pixels','menubar','none','numbertitle',...
    'off','resize','off','name','Convert to VTK','Visible','Off');

f_size = [300 256+s];
f_export.Position(3:4) = f_size;
movegui(f_export,'center')

% default values
    h.export.T = 1;
    h.export.v = 1;
    h.export.p = 0;
    h.domain   = [1 0];
    h.n        = 100;

% windows vs linux
    if ispc
        h.sep = '\'; % separator for paths
    else
        h.sep = '/';
    end
    
% cartesian vs polar
    ax = evalin('base','flowopt.ax');
    if ax == 1
        dir = 'azimuthal';
        txt = '$n_{\varphi} = $';
    else
        dir = 'spanwise';
        txt = '$n_z = $';
    end
    
% energy equation solved?
    energy = evalin('base','flowopt.energy');
    
% SFM or TFM
    n_blocks = evalin('base','length(blocks)');

h.pl.var = uipanel('unit','pix','pos',[9 148 128 100+s],'title','Field Quantity');
    h.cb.T = uicontrol('parent',h.pl.var,'style','check','pos',[10 58+s 120 20],'value',h.export.T,'String','Temperature');
    if strcmp(flow,'pF')
        h.cb.th_E = uicontrol('parent',h.pl.var,'style','check','pos',[10 58 120 20],'value',h.export.T,'String','Th. Production');
    end
    h.cb.v = uicontrol('parent',h.pl.var,'style','check','pos',[10 33 120 20],'value',h.export.v,'String','Velocity');
    h.cb.p = uicontrol('parent',h.pl.var,'style','check','pos',[10 8 120 20],'value',h.export.p,'String','Pressure');
    if energy == 0
        set(h.cb.T,'value',0,'enable','off')
        if strcmp(flow,'pF')
            set(h.cb.th_E,'value',0,'enable','off')
        end
    end
    
h.bg.domain = uibuttongroup('unit','pix','pos',[145 148 148 100+s],'title','Domain');
    h.rb.LB = uicontrol('style','rad','units','pix','parent',h.bg.domain,'pos',[10 48+3*s/4 110 22],'string','Liquid Phase','value',1);
    h.rb.SG = uicontrol('style','rad','units','pix','parent',h.bg.domain,'pos',[10 22+s/4 130 22],'string','Gas Phase');
    if n_blocks == 1
        set(h.rb.SG,'enable','off');
    end
    
h.pl.n = uipanel('unit','pix','pos',[9 90 285 50],'title',['Number of points in ' dir ' direction']);
    h.et.n = uicontrol('style','edit','units','pix',...
        'parent',h.pl.n,'string',num2str(h.n),'pos',[60 7 75 23],'call',{@et_n,h});
    axes('unit','pix','visible','off','parent',h.pl.n,'pos',h.et.n.Position);
    if ispc
        text(-0.09,0.5,txt,'interpreter','Latex','FontSize',12,'HorizontalAlignment','right')
    else
        text(-0.09,0.47,txt,'interpreter','Latex','FontSize',12,'HorizontalAlignment','right')
    end
    
h.et.path = uicontrol('style','edit','units','pix','pos',[8 40 285 23],'HorizontalAlignment','left');
    if evalin('base','exist(''file_path'',''var'')') == 1 && evalin('base','strcmp(file_path(end-3:end),''.vtk'')')
        set(h.et.path,'string',evalin('base','file_path'))
    else
        set(h.et.path,'string',[pwd h.sep 'newData.vtk'])
    end
    uicontrol('style','text','units','pix','HorizontalAlignment','left','pos',[8 64 80 17],'string','File name');
    
h.pb.export = uicontrol('units','pix','string','Export',...
    'pos',[216 8 76 25],'ForegroundColor',[0.8500 0.3250 0.0980],...
    'call',{@pb_create_vtk,h});

% correct figure size
    pause(0.1)
    f_export.Position(3:4) = f_size;
    set(f_export,'Visible','On')

guidata(f_export,h)


    function et_n(hObject, ~, ~)
        h = guidata(hObject);
        
        h.n = getValue(1,inf,h.n,hObject);

        if floor(h.n) ~= h.n || h.n <= 0
            waitfor(warndlg('Value must be a positive integer.','Warning'))
            h.n = floor(abs(h.n));
            set(h.et.n,'string',num2str(h.n))
        end

        guidata(hObject,h)
    end


    function pb_create_vtk(hObject, ~, ~)
        h = guidata(hObject);
        
        % get exporting variables
            assignin('base','Nphi',h.n);

            h.export.T = get(h.cb.T,'value');
            h.export.v = get(h.cb.v,'value');
            h.export.p = get(h.cb.p,'value');
            if strcmp(flow,'pF')
                h.export.th_E = get(h.cb.th_E,'value');
            else
                h.export.th_E = 0;
            end
            
            if h.export.T+h.export.v+h.export.p+h.export.th_E == 0
                waitfor(warndlg('Choose at least one field quantity to export.','Warning'))
                return
            else
                assignin('base','tmp',h.export);
                evalin('base','export = tmp;');
            end
            
        % get domain
            str = get(get(h.bg.domain,'SelectedObject'),'String');
            switch str
                case 'Liquid Phase'
                    h.domain = [1 0];
                case 'Gas Phase'
                    h.domain = [0 1];
            end
            assignin('base','domain',h.domain);
            
        % get file name
            f_path = get(h.et.path,'string');
            idx = strfind(f_path,h.sep);
            path = f_path(1:idx(end));
            
            if ~strcmp(f_path(end-3:end),'.vtk')
                waitfor(errordlg('Invalid file name extension!','Error'))
                return
            elseif ~isfolder(path(1:idx(end)))
                waitfor(errordlg('Invalid Path! Directory does not exist','Error'))
                return
            elseif isfile(f_path)
                selection = questdlg({['A file named "' f_path(idx(end)+1:end) '" already exists!'],'Do you want to replace it?'},'','Cancel','Replace','Replace');
                switch selection
                    case 'Cancel'
                        return
                    case 'Replace'
                        assignin('base','file_path',f_path);
                end
            else
                assignin('base','file_path',f_path);
            end
            
        % convert to vtk
            evalin('base','convert_to_vtk')
            delete(f_export)

    end

end