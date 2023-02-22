unity = 'mm'; dim = 1;

mesh.figure = figure('NumberTitle', 'off', 'Name', 'Mesh'); hold on;

% unity
if dim == 0
    f = 1/l_lb;
else
    switch unity
        case 'm'
            f = 1;
        case 'dm'
            f = 10;
        case 'cm'
            f = 100;
        case 'mm'
            f = 1000;
    end
end

% make figure visible after all properties are set
    set(mesh.figure,'visible','off')
    % change figure width and height
        mesh.figurePosition = mesh.figure.Position;

    % change figure position
        movegui(mesh.figure,'center')


axis('equal')

% grid continuity/thermal energy
    % Horizontal grid lines
    shift = b1.geom.z(1)+0.5*diff(b1.geom.z);
    plot(b1.R.v*f,(shift-b1.Z.v)*f,'color',[20 180 200]./255)
    if length(blocks) == 2
        plot(b2.R.v*f,(shift-b2.Z.v)*f,'color',[1 1 1]*100/255)
    end
    
    % Vertical grid lines
    plot(b1.R.v'*f,(shift-b1.Z.v')*f,'color',[20 180 200]./255)
    if length(blocks) == 2
        plot(b2.R.v'*f,(shift-b2.Z.v')*f,'color',[1 1 1]*100/255)
    end

if length(blocks) == 2
    axis([0 geom.r(end) shift-geom.z(end) shift-geom.z(1)]*f)
else
    axis([0 1.2*max(b1.R.v(:)) -l_lb/2 l_lb/2]*f)
end

% legends
p = flowopt.ax+1; rx = 'xr'; zy = 'yz';
if dim == 1
    xlabel(['$' rx(p) '$ [' unity ']'],'Interpreter','latex','FontSize',12)
    ylabel(['$' zy(p) '$ [' unity ']'],'Interpreter','latex','FontSize',12)
else
    xlabel(['$' rx(p) '$'],'Interpreter','latex','FontSize',12)
    ylabel(['$' zy(p) '$'],'Interpreter','latex','FontSize',12)
end

set(gcf,'color','white');
box on

set(mesh.figure,'visible','on')