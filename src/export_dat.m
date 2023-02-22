function export_dat(quantities)

% number of possible quantities
    h.n_q = length(quantities);
    
% panel's height
    h_box = 50+25*(h.n_q-1);

f_export = figure('units','pixels','menubar','none','numbertitle',...
    'off','resize','off','name','Convert to DAT','Visible','Off');

f_size = [300 98+h_box];
f_export.Position(3:4) = f_size;
movegui(f_export,'center')

% windows vs linux
    if ispc
        h.sep = '\'; % separator for paths
    else
        h.sep = '/';
    end

h.pl.var = uipanel('unit','pix','pos',[9 90 285 h_box],'title','Field Quantity');
    for i = 1:h.n_q
        h.cb(i) = uicontrol('parent',h.pl.var,'style','check','pos',...
            [10 10+25*(h.n_q-i) 200 20],'value',1,'String',quantities{i});
    end

h.et.path = uicontrol('style','edit','units','pix','pos',[8 40 285 23],'HorizontalAlignment','left');
    if evalin('base','exist(''file_path'',''var'')') == 1 && evalin('base','strcmp(file_path(end-3:end),''.dat'')')
        set(h.et.path,'string',evalin('base','file_path'))
    else
        set(h.et.path,'string',[pwd h.sep 'newData.dat'])
    end
    uicontrol('style','text','units','pix','HorizontalAlignment','left','pos',[8 64 80 17],'string','File name');
    
h.pb.export = uicontrol('units','pix','string','Export',...
    'pos',[216 8 76 25],'ForegroundColor',[0 0.4470 0.7410],...
    'call',{@pb_create_dat,h});

% correct figure size
    pause(0.1)
    f_export.Position(3:4) = f_size;
    set(f_export,'Visible','On')

guidata(f_export,h)


    function pb_create_dat(hObject, ~, ~)
        h = guidata(hObject);
        
        % get exporting variables n_e
            n_e = zeros(h.n_q,1);
            for j = 1:h.n_q
                n_e(j) = get(h.cb(j),'value')*j;
            end
            n_e = nonzeros(n_e);
            
            if isempty(n_e)
                waitfor(warndlg('Choose at least one field quantity to export.','Warning'))
                return
            else
                e_quantities = quantities(n_e);
                e_quantities = convertCharsToStrings(e_quantities);
                assignin('base','tmp',e_quantities);
                evalin('base','lp.e_quantity = tmp;');
            end
        
        % get file name
            f_path = get(h.et.path,'string');
            idx = strfind(f_path,h.sep);
            path = f_path(1:idx(end));
            
            if ~strcmp(f_path(end-3:end),'.dat')
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
            
        % convert to dat
            evalin('base','convert_to_dat')
            delete(f_export)

    end

end