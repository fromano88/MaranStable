function Temperature = inputTemperature_SFM(coordinate,r_c,ax)

if nargin == 0
    coordinate = 'r'; r_c = 1; ax = 1;
end
p  = ax+1; % pointer
rx = 'xr';
h.coordinate = coordinate;

Temperature = [];
f_temperature = figure('units','pixels','menubar','none','numbertitle',...
    'off','name','Function Editor','resize','off','Visible','Off');

f_size = [333 295];
f_temperature.Position(3:4) = f_size;
movegui(f_temperature,'center')

background = uipanel('parent',f_temperature,'unit','pix','pos',[6 4 324 288],'BackgroundColor',[1 1 1]);
Color = [1 0.7 0.7];

B = {'parent',background,'units','pix','pos'};

for j = 1:3
    for i = 1:3
    h.pb(i+(j-1)*3) = uicontrol(B{:},[(81*i-70) 48+(j-1)*37 56 27],'string',num2str(i+(j-1)*3));
    end
end

h.dot = uicontrol(B{:},[11 11 56 27],'string','.','enable','off','call',{@dot_call,h});
h.pb(10) = uicontrol(B{:},[92 11 56 27],'string','0');
h.var = uicontrol(B{:},[173 11 56 27],...
    'string',coordinate,'BackgroundColor',Color,'call',{@var_call,h});
h.Divide = uicontrol(B{:},[254 11 56 27],'string','/','enable','off');
h.Multiply = uicontrol(B{:},[254 48 56 27],'string','*','enable','off');
h.Minus = uicontrol(B{:},[254 85 56 27],'string','-');
h.Plus = uicontrol(B{:},[254 122 56 27],'string','+','enable','off');
h.Power = uicontrol(B{:},[11 159 56 27],'string','^','enable','off');
h.br.left = uicontrol(B{:},[92 159 56 27],'string','(');
h.br.right = uicontrol(B{:},[173 159 56 27],'string',')','enable','off','call',{@var_call,h});
h.Clear = uicontrol(B{:},[254 159 56 27],'string','Clear','call',{@clear_call,h});
h.lb = uicontrol(B{:},[11 196 56 27],'string','d','call',{@var_call,h});
h.rc = uicontrol(B{:},[92 196 56 27],'string',[rx(p) '_c'],'call',{@var_call,h});
h.ri = uicontrol(B{:},[173 196 56 27],'string',[rx(p) '_i'],'call',{@var_call,h});
h.Delete = uicontrol(B{:},[254 196 56 27],'string','DEL','call',{@del_call,h});
h.enter = uicontrol(B{:},[254 239 56 27],'string','Enter','BackgroundColor',Color,'enable','off','call',{@enter_call,h});
h.display.pl = uipanel(B{:},[11 239 218 27],'BackgroundColor',Color);
h.display.st = uicontrol('parent',h.display.pl,'units','pix','style',...
    'text','HorizontalAlignment','center','pos',[4 1 208 19],...
    'BackgroundColor',Color);

if r_c == 0
    set(h.rc,'enable','off')
end

h.dsp = [];
h.count.r = 0; h.count.l = 0;
h.del = 0;

set([h.pb(1) h.pb(2) h.pb(3) h.pb(4) h.pb(5) h.pb(6) h.pb(7) h.pb(8) h.pb(9) h.pb(10)],'call',{@numbers,h})
set([h.Divide h.Multiply h.Minus h.Plus h.Power h.br.left],'call',{@operators,h})
set(f_temperature,'CloseRequestFcn',{@fig_delet,h})

% correct figure size
    pause(0.1)
    f_temperature.Position(3:4) = f_size;
    set(f_temperature,'Visible','On')

uiwait(f_temperature)

    function numbers(varargin)
        if h.del == 0
            dsp = get(varargin{1},'string'); h.dsp = strcat(h.dsp,dsp);
            if length(h.dsp) > 24
                set(h.display.st,'String',h.dsp(end-24:end)) % cut string if it is too long
            else
                set(h.display.st,'string',h.dsp)
            end
        end

        set(h.pb(:),'enable','on'); set([h.Divide h.Multiply h.Minus h.Plus h.dot h.Power],'enable','on')
        set([h.br.left h.var h.lb h.rc h.ri],'Enable','off')

        % right bracket can only be pressed if the number of opened brackets is 
        % higher or rather not equal to number of closed brackets    
        if count(h.dsp,'(') == count(h.dsp,')')
            set(h.br.right,'Enable','off')
        else
            set(h.br.right,'Enable','on')
        end

        if varargin{1} == h.pb(10) && (length(h.dsp) == 1 || strcmp(h.dsp(end-1),'/'))
            set(h.pb(:),'enable','off'); set([h.Divide h.Multiply h.Minus h.Plus h.br.right],'enable','off')
        end

        if count(h.dsp,h.coordinate) > count(h.dsp,[h.coordinate '_']) && count(h.dsp,'(') == count(h.dsp,')')
            set(h.enter,'Enable','on')
        else
            set(h.enter,'Enable','off')
        end
    end

    function operators(varargin)
        if h.del == 0
            dsp = get(varargin{1},'string'); h.dsp = strcat(h.dsp,dsp);
            if length(h.dsp) > 24
                set(h.display.st,'String',h.dsp(end-24:end)) % cut string if it is too long
            else
                set(h.display.st,'string',h.dsp)
            end
        end

        set(h.pb(:),'enable','on'); set([h.br.left h.var h.lb h.rc h.ri],'Enable','on')
        set([h.Divide h.Multiply h.Minus h.Plus h.Power h.dot h.br.right h.enter],'enable','off')

        if varargin{1} == h.Power || varargin{1} == h.br.left
            set(h.Minus,'enable','on')
        end

        h.count.r = 0; h.count.l = 0;
    end

    function dot_call(varargin)
        if h.del == 0
            dsp = get(varargin{1},'string'); h.dsp = strcat(h.dsp,dsp);
            if length(h.dsp) > 24
                set(h.display.st,'String',h.dsp(end-24:end)) % cut string if it is too long
            else
                set(h.display.st,'string',h.dsp)
            end
        end

        set(h.pb(:),'enable','on');
        set([h.Divide h.Multiply h.Minus h.Plus h.br.left h.var h.lb h.rc h.ri h.Power h.dot h.br.right h.enter],'enable','off')
    end

    function var_call(varargin)
        if h.del == 0
            dsp = get(varargin{1},'string'); h.dsp = strcat(h.dsp,dsp);
            if length(h.dsp) > 24
                set(h.display.st,'String',h.dsp(end-24:end)) % cut string if it is too long
            else
                set(h.display.st,'string',h.dsp)
            end
        end

        set([h.Divide h.Multiply h.Minus h.Plus h.Power],'enable','on');
        set(h.pb(:),'enable','off'); set([h.br.left h.var h.lb h.rc h.ri h.dot],'enable','off')

        if count(h.dsp,'(') == count(h.dsp,')')
            set(h.br.right,'Enable','off')
        else
            set(h.br.right,'Enable','on')
        end
        if count(h.dsp,h.coordinate) > count(h.dsp,[h.coordinate '_']) && count(h.dsp,'(') == count(h.dsp,')')
            set(h.enter,'Enable','on')
        else
            set(h.enter,'Enable','off')
        end
    end

    function clear_call(varargin)
        h.dsp = []; set(h.display.st,'string',h.dsp)

        set(h.pb(:),'enable','on'); set([h.br.left h.var h.lb h.rc h.ri h.Minus],'enable','on')
        set([h.Divide h.Multiply h.Minus h.Plus h.Power h.dot h.br.right h.enter],'enable','off')

        h.count.r = 0; h.count.l = 0;
    end

    function del_call(varargin)
        if ~isempty(h.dsp)
            if max(strcmp(h.dsp(end),{'c','i'}))
                if length(h.dsp) > 3
                    h.dsp = h.dsp(1:end-3);
                else
                    h.dsp = '';
                end
            else
                h.dsp = h.dsp(1:end-1);
            end
        end
        if length(h.dsp) > 24
            set(h.display.st,'String',h.dsp(end-24:end)) % cut string if it is too long
        else
            set(h.display.st,'string',h.dsp)
        end

        % enable/disable buttons
        h.del = 1;
        guidata(varargin{1},h)
        if isempty(h.dsp)
            clear_call(varargin{1})
        elseif max(strcmp(h.dsp(end),{'0','1','2','3','4','5','6','7','8','9'}))
            numbers(varargin{1})
        elseif strcmp(h.dsp(end),'.')
            dot_call(varargin{1})
        elseif max(strcmp(h.dsp(end),{h.coordinate,')','c','i','b'}))
            var_call(varargin{1})
        elseif max(strcmp(h.dsp(end),{'+','-','*','/','^','('}))
            operators(varargin{1})
        end
        h.del = 0;
    end

    function enter_call(varargin)
        Temperature = h.dsp;
        delete(f_temperature)
    end

    function fig_delet(varargin)
        if isempty(get(h.display.st,'String'))
            selection = questdlg('Empty input! Close anyway?','','Yes','No','No');
        else
            selection = questdlg('Input not confirmed! Close anyway?','','Yes','No','No');
        end
        switch selection
            case 'Yes'
                delete(f_temperature)
            case 'No'
                return
        end
    end
end