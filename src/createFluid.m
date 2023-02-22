function [name, prop, fluid] = createFluid(model,h)

if nargin == 0
    model = 'TFM';
end

% default values
    h.mu(1)     = 1e-3; h.mu(2)     = 0;
    h.rho(1)    = 1000; h.rho(2)    = 0;
    h.lambda(1) = 0.1;  h.lambda(2) = 0;
    h.cp(1)     = 1800; h.cp(2)     = 0;
    h.sigma(1)  = 0.02; h.sigma(2)  = 7e-5;

% default outputs
    name = 'MyFluid'; prop = []; fluid = 'liquid';

f_create = figure('units','pixels','menubar','none','numbertitle',...
    'off','resize','off','Name','Create Fluid','Visible','Off');

f_size = [480 475];
f_create.Position(3:4) = f_size;
movegui(f_create,'center')

% save writings
    h.latex_l = {'interpreter','Latex','FontSize',12,'HorizontalAlignment','right'};
    h.latex_r = {'interpreter','Latex','FontSize',12,'HorizontalAlignment','left'};
    E = {'style','edit','units','pix','parent'};
    A = {'unit','pix','visible','off','parent'};
    R = {'style','rad','units','pix','parent'};

% text position for axes
    h.tp_l = -0.05; h.tp_r = 1.1;
    if ispc
        h.tp_y = 0.5;
    else
        h.tp_y = 0.47;
    end

% panels
    ph = 70; % panel height
    ps = 41; % panel start (y coordinate)
    h.pl.name   = uipanel('unit','pix','pos',[6 ps+5*5+5*ph 232 55],'title','Name');
    h.pl.mu     = uipanel('unit','pix','pos',[6 ps+4*5+4*ph 470 ph],'title','Dynamic Viscosity');
    h.pl.rho    = uipanel('unit','pix','pos',[6 ps+3*5+3*ph 470 ph],'title','Density');
    h.pl.lambda = uipanel('unit','pix','pos',[6 ps+2*5+2*ph 470 ph],'title','Thermal Conductivity');
    h.pl.cp     = uipanel('unit','pix','pos',[6 ps+1*5+1*ph 470 ph],'title','Specific Heat Capacity');
    h.pl.sigma  = uipanel('unit','pix','pos',[6 ps          470 ph],'title','Surface Tension');

% button groups
    h.bg.fluid = uibuttongroup('units','pix','pos',[244 ps+5*5+5*ph 232 55],'title','Fluid','SelectionChangedFcn',{@bg_fluid,h});

% radio buttons
    h.rb.liquid = uicontrol(R{:},h.bg.fluid,'pos',[ 30 10 80 22],'string','Liquid','value',1);
    h.rb.gas    = uicontrol(R{:},h.bg.fluid,'pos',[140 10 80 22],'string','Gas');
    if strcmp(model,'SFM')
        set(h.rb.gas,'enable','off')
    end

% edit texts
    h.et.name = uicontrol(E{:},h.pl.name,'pos',[15 10 180 23],'string','MyFluid','HorizontalAlignment','left','call',{@et_name,h});
    h.et.mu(1)     = uicontrol(E{:},h.pl.mu,    'pos',[ 80 17 70 23],     'string',num2str(h.mu(1)),    'call',{@et_value,h});
    h.et.mu(2)     = uicontrol(E{:},h.pl.mu,    'pos',[205 17 70 23],     'string',num2str(h.mu(2)),    'call',{@et_value,h});
    h.et.rho(1)    = uicontrol(E{:},h.pl.rho,   'pos',h.et.mu(1).Position,'string',num2str(h.rho(1)),   'call',{@et_value,h});
    h.et.rho(2)    = uicontrol(E{:},h.pl.rho,   'pos',h.et.mu(2).Position,'string',num2str(h.rho(2)),   'call',{@et_value,h});
    h.et.lambda(1) = uicontrol(E{:},h.pl.lambda,'pos',h.et.mu(1).Position,'string',num2str(h.lambda(1)),'call',{@et_value,h});
    h.et.lambda(2) = uicontrol(E{:},h.pl.lambda,'pos',h.et.mu(2).Position,'string',num2str(h.lambda(2)),'call',{@et_value,h});
    h.et.cp(1)     = uicontrol(E{:},h.pl.cp,    'pos',h.et.mu(1).Position,'string',num2str(h.cp(1)),    'call',{@et_value,h});
    h.et.cp(2)     = uicontrol(E{:},h.pl.cp,    'pos',h.et.mu(2).Position,'string',num2str(h.cp(2)),    'call',{@et_value,h});
    h.et.sigma(1)  = uicontrol(E{:},h.pl.sigma, 'pos',h.et.mu(1).Position,'string',num2str(h.sigma(1)), 'call',{@et_value,h});
    h.et.sigma(2)  = uicontrol(E{:},h.pl.sigma, 'pos',[177 17 70 23],     'string',num2str(h.sigma(2)), 'call',{@et_value,h});

% pushbuttons
    h.pb.cancel = uicontrol('units','pix','string','Cancel','pos',[316 8 76 25],'call',{@pb_cancel,h});
    h.pb.create = uicontrol('units','pix','string','Create','pos',[398 8 76 25],...
        'ForegroundColor',[0.8500 0.3250 0.0980],'call',{@pb_create,h});

% axes
    h.as.mu(1)     = axes(A{:},h.pl.mu,'pos',h.et.mu(1).Position);
    h.as.mu(2)     = axes(A{:},h.pl.mu,'pos',h.et.mu(2).Position);
    h.as.mu(3)     = axes(A{:},h.pl.mu,'pos',[225 17 70 23]);
    h.as.rho(1)    = axes(A{:},h.pl.rho,'pos',h.et.rho(1).Position);
    h.as.rho(2)    = axes(A{:},h.pl.rho,'pos',h.et.rho(2).Position);
    h.as.rho(3)    = axes(A{:},h.pl.rho,'pos',[225 17 70 23]);
    h.as.lambda(1) = axes(A{:},h.pl.lambda,'pos',h.et.lambda(1).Position);
    h.as.lambda(2) = axes(A{:},h.pl.lambda,'pos',h.et.lambda(2).Position);
    h.as.lambda(3) = axes(A{:},h.pl.lambda,'pos',[225 17 70 23]);
    h.as.cp(1)     = axes(A{:},h.pl.cp,'pos',h.et.cp(1).Position);
    h.as.cp(2)     = axes(A{:},h.pl.cp,'pos',h.et.cp(2).Position);
    h.as.cp(3)     = axes(A{:},h.pl.cp,'pos',[225 17 70 23]);
    h.as.sigma(1)  = axes(A{:},h.pl.sigma,'pos',h.et.sigma(1).Position);
    h.as.sigma(2)  = axes(A{:},h.pl.sigma,'pos',h.et.sigma(2).Position);
    h.as.sigma(3)  = axes(A{:},h.pl.sigma,'pos',[225 17 70 23]);

    axes(h.as.mu(1)); cla(h.as.mu(1));
    text(h.tp_l,h.tp_y,'$\mu_\mathrm{l}(T)=$ ',h.latex_l{:})
    text(h.tp_r,h.tp_y,'$\cdot\; \Big(1\,-$ ',h.latex_r{:})
    %
    axes(h.as.mu(2)); cla(h.as.mu(2));
    text(h.tp_r,h.tp_y,'$\cdot\; (T\,-25^\circ \mathrm{C} )\Big)$ ',h.latex_r{:})
    %
    axes(h.as.mu(3)); cla(h.as.mu(3));
    text(3,h.tp_y,'$\left[\displaystyle\mathrm{Pa}\,\mathrm{s}\right]$ ',h.latex_r{:},'HorizontalAlignment','center')
    %
    %
    %
    axes(h.as.rho(1)); cla(h.as.rho(1));
    text(h.tp_l,h.tp_y,'$\rho_\mathrm{l}(T)=$ ',h.latex_l{:})
    text(h.tp_r,h.tp_y,'$\cdot\; \Big(1\,-$ ',h.latex_r{:})
    %
    axes(h.as.rho(2)); cla(h.as.rho(2));
    text(h.tp_r,h.tp_y,'$\cdot\; (T\,-25^\circ \mathrm{C} )\Big)$ ',h.latex_r{:})
    %
    axes(h.as.rho(3)); cla(h.as.rho(3));
    text(3,h.tp_y,'$\left[\displaystyle\frac{\mathrm{kg}}{\mathrm{m}^3}\right]$ ',h.latex_r{:},'HorizontalAlignment','center')
    %
    %
    %
    axes(h.as.lambda(1)); cla(h.as.lambda(1));
    text(h.tp_l,h.tp_y,'$\lambda_\mathrm{l}(T)=$ ',h.latex_l{:})
    text(h.tp_r,h.tp_y,'$\cdot\; \Big(1\,+$ ',h.latex_r{:})
    %
    axes(h.as.lambda(2)); cla(h.as.lambda(2));
    text(h.tp_r,h.tp_y,'$\cdot\; (T\,-25^\circ \mathrm{C} )\Big)$ ',h.latex_r{:})
    %
    axes(h.as.lambda(3)); cla(h.as.lambda(3));
    text(3,h.tp_y,'$\left[\displaystyle\frac{\mathrm{W}}{\mathrm{m K}}\right]$ ',h.latex_r{:},'HorizontalAlignment','center')
    %
    %
    %
    axes(h.as.cp(1)); cla(h.as.cp(1));
    text(h.tp_l,h.tp_y,'$c_{p\mathrm{l}}(T)=$ ',h.latex_l{:})
    text(h.tp_r,h.tp_y,'$\cdot\; \Big(1\,+$ ',h.latex_r{:})
    %
    axes(h.as.cp(2)); cla(h.as.cp(2));
    text(h.tp_r,h.tp_y,'$\cdot\; (T\,-25^\circ \mathrm{C} )\Big)$ ',h.latex_r{:})
    %
    axes(h.as.cp(3)); cla(h.as.cp(3));
    text(3,h.tp_y,'$\left[\displaystyle\frac{\mathrm{J}}{\mathrm{kg K}}\right]$ ',h.latex_r{:},'HorizontalAlignment','center')
    %
    axes(h.as.sigma(1)); cla(h.as.sigma(1));
    text(h.tp_l,h.tp_y,'$\sigma (T)=$ ',h.latex_l{:})
    text(h.tp_r,h.tp_y,'$-$ ',h.latex_r{:})
    %
    axes(h.as.sigma(2)); cla(h.as.sigma(2));
    text(h.tp_r,h.tp_y,'$\cdot\; (T\,-25^\circ \mathrm{C} )$ ',h.latex_r{:})
    %
    axes(h.as.sigma(3)); cla(h.as.sigma(3));
    text(3,h.tp_y,'$\left[\displaystyle\frac{\mathrm{N}}{\mathrm{m}}\right]$ ',h.latex_r{:},'HorizontalAlignment','center')

% function called if user closes the window or presses 'cancel'
    set(f_create,'CloseRequestFcn',{@pb_cancel,h})

% correct figure size
    pause(0.1)
    f_create.Position(3:4) = f_size;
    movegui(f_create,'center')
    set(f_create,'Visible','On')

uiwait(f_create)

    function et_name(varargin)
        name = get(varargin{1},'string');
        if ~isletter(name(1))
            waitfor(errordlg('Invalid name. First character must be a letter.','Error'))
            name = 'MyFluid';
            set(varargin{1},'string',name);
        elseif ismember(name,who) || strcmp(name,'b')
            waitfor(errordlg('Invalid name!','Error'))
            name = 'MyFluid';
            set(varargin{1},'string',name);
        end
    end

    function bg_fluid(varargin)
        fluid = lower(varargin{1}.SelectedObject.String);
        f = fluid(1);

        h.mu(2)     = 0;
        h.rho(2)    = 0;
        h.lambda(2) = 0;
        h.cp(2)     = 0;

        switch fluid
            case 'liquid'
                h.mu(1)     = 1e-3;
                h.rho(1)    = 1000;
                h.lambda(1) = 0.1;
                h.cp(1)     = 1800;

                set(h.et.sigma(:),'enable','on')
            case 'gas'
                h.mu(1)     = 1e-5;
                h.rho(1)    = 1;
                h.lambda(1) = 0.01;
                h.cp(1)     = 1000;

                set(h.et.sigma(:),'enable','off')
        end
        set(h.et.mu(1),'string',num2str(h.mu(1)))
        set(h.et.rho(1),'string',num2str(h.rho(1)))
        set(h.et.lambda(1),'string',num2str(h.lambda(1)))
        set(h.et.cp(1),'string',num2str(h.cp(1)))

        set(h.et.mu(2),'string',num2str(h.mu(2)))
        set(h.et.rho(2),'string',num2str(h.rho(2)))
        set(h.et.lambda(2),'string',num2str(h.lambda(2)))
        set(h.et.cp(2),'string',num2str(h.cp(2)))

        axes(h.as.rho(1)); cla(h.as.rho(1));
        text(h.tp_l,h.tp_y,['$\rho_\mathrm{' f '}(T)=$'],h.latex_l{:})
        text(h.tp_r,h.tp_y,'$\cdot\; \Big(1\,-$ ',h.latex_r{:})
        %
        axes(h.as.lambda(1)); cla(h.as.lambda(1));
        text(h.tp_l,h.tp_y,['$\lambda_\mathrm{' f '}(T)=$'],h.latex_l{:})
        text(h.tp_r,h.tp_y,'$\cdot\; \Big(1\,+$ ',h.latex_r{:})
        %
        axes(h.as.cp(1)); cla(h.as.cp(1));
        text(h.tp_l,h.tp_y,['$c_{p\mathrm{' f '}}(T)=$'],h.latex_l{:})
        text(h.tp_r,h.tp_y,'$\cdot\; \Big(1\,+$ ',h.latex_r{:})
    end

    function et_value(varargin)
        object = varargin{1};
        switch varargin{1}
            case h.et.mu(1)
                default = h.mu(1);
                [h.mu(1),~] = getValue(-inf,inf,default,object);
            case h.et.rho(1)
                default = h.rho(1);
                [h.rho(1),~] = getValue(-inf,inf,default,object);
            case h.et.lambda(1)
                default = h.lambda(1);
                [h.lambda(1),~] = getValue(-inf,inf,default,object);
            case h.et.cp(1)
                default = h.cp(1);
                [h.cp(1),~] = getValue(-inf,inf,default,object);
            case h.et.sigma(1)
                default = h.sigma(1);
                [h.sigma(1),~] = getValue(-inf,inf,default,object);

            case h.et.mu(2)
                default = h.mu(2);
                [h.mu(2),~] = getValue(-inf,inf,default,object);
            case h.et.rho(2)
                default = h.rho(2);
                [h.rho(2),~] = getValue(-inf,inf,default,object);
            case h.et.lambda(2)
                default = h.lambda(2);
                [h.lambda(2),~] = getValue(-inf,inf,default,object);
            case h.et.cp(2)
                default = h.cp(2);
                [h.cp(2),~] = getValue(-inf,inf,default,object);
            case h.et.sigma(2)
                default = h.sigma(2);
                [h.sigma(2),~] = getValue(-inf,inf,default,object);
        end
    end

    function pb_create(varargin)
        if ismember(name,h.liquid_list) || (strcmp(model,'TFM') && ismember(name,h.gas_list))
            waitfor(warndlg('Fluid cannot be created. Name already exists.','Attention'))
            return
        end
        prop.mu     = [num2str(h.mu(1)) '.*(1-' num2str(h.mu(2)) '*(theta-25))'];
        prop.rho    = [num2str(h.rho(1)) '.*(1-' num2str(h.rho(2)) '*(theta-25))'];
        prop.lambda = [num2str(h.lambda(1)) '.*(1+' num2str(h.lambda(2)) '*(theta-25))'];
        prop.cp     = [num2str(h.cp(1)) '.*(1+' num2str(h.cp(2)) '*(theta-25))'];
        
        prop.dmu     = num2str(-h.mu(1)*h.mu(2));
        prop.drho    = num2str(-h.rho(1)*h.rho(2));
        prop.dlambda = num2str(h.lambda(1)*h.lambda(2));
        prop.dcp     = num2str(h.cp(1)*h.cp(2));
        
        if strcmp(fluid,'liquid')
            prop.sigma  = [num2str(h.sigma(1)) '-' num2str(h.sigma(2)) '*(theta-25)'];
            prop.dsigma = num2str(-h.sigma(2));

        end

        delete(f_create)

    end

    function pb_cancel(varargin)
        name = [];
        delete(f_create)
    end
end