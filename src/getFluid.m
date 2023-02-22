function selectedFluid = getFluid(fluid,h)

if nargin == 0
    fluid = 'liquid';
end

selectedFluid = [];
f_fluid = figure('units','pixels','menubar','none','numbertitle',...
    'off','resize','off','Visible','Off');

f_size = [160 330];
f_fluid.Position(3:4) = f_size;
movegui(f_fluid,'center')

h.lb = uicontrol('parent',f_fluid,'style','listbox','units','pix','pos',[6 4 150 322]);

if strcmp(fluid,'liquid')
    set(f_fluid,'name','Liquids')
    set(h.lb,'string',h.liquid_list)
else
    set(f_fluid,'name','Gases')
    set(h.lb,'string',h.gas_list)
end

set(h.lb,'call',{@fluid_list,h})

% correct figure size
    pause(0.1)
    f_fluid.Position(3:4) = f_size;
    set(f_fluid,'Visible','On')


uiwait(f_fluid)

    function fluid_list(varargin)
       if strcmp(get(gcf,'selectiontype'),'open')  
          str = get(h.lb,'string');
          val = get(h.lb,'value');
          selectedFluid = str{val};
          close(f_fluid)
       end
    end
end
