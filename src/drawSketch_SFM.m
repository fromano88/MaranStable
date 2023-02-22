V_r = h.V_r;
g   = h.flowopt.g; ax = h.flowopt.ax;
surface = h.b1.bc.z(2,2);
r_c     = h.r_c;
delta_T = h.temperature.ri_d1-h.temperature.ri_d2;
energy  = h.flowopt.energy;

%--- change parameters in order to change the sketch ---%
w = 1.5;
r_i = 3.5; r_C = r_i*0.2;
l_d2 = 3.2; l_lb = 2.8; l_d1 = 3.4;
d = 0.1;
l.orgn = r_i/4; l.g = 1;
vPos.rc = 0.63*l_lb; vPos.ri = 0.4*l_d2;
n.lines = 15; n.colors = 100;
color.liquid = [181 255 255]./255; color.liquid_light = [218 255 255]./255;

if energy == 0
    color.up   = 0.3+[0 0 0];
    color.down = 0.3+[0 0 0];
    color.up_light   = 0.9+[0 0 0];
    color.down_light = 0.9+[0 0 0];
elseif delta_T >= 0
    color.up   = [0.8500 0.3250 0.0980];
    color.down = [0 0.4470 0.7410];
    color.up_light   = [255 145 98]/255;
    color.down_light = [50 174 255]/255;
else
    color.up   = [0 0.4470 0.7410];
    color.down = [0.8500 0.3250 0.0980];
    color.up_light   = [50 174 255]/255;
    color.down_light = [255 145 98]/255;
end
%-------------------------------------------------------%

% max and min values for plotting
if V_r > 1.6
    V_r = tanh(V_r+0.4)+1.6-tanh(2);
elseif V_r < 0.3
    V_r = tanh(V_r-2.3)+0.3+tanh(2);
end
if g > 30 || g < -30
    g = 30*tanh(g)/tanh(30);
end

% surface shape
if strcmp(surface,'r') % rigid
    a = 0; k = 0;
else
    if g == 0
        a = 5*r_i/l_lb^2-r_i*sqrt(30*V_r-5)/l_lb^2;
    else
        a = -0.025*g;
        if g > 0
            k = l_lb/2-5*r_i/(a*l_lb^2)-sqrt(5*r_i^2*(6*V_r-1)/(a^2*l_lb^4)-l_lb^2/28);
        elseif g < 0
            k = l_lb/2-5*r_i/(a*l_lb^2)+sqrt(5*r_i^2*(6*V_r-1)/(a^2*l_lb^4)-l_lb^2/28);
        end
    end
    
end
if g == 0
    f = @(x) a*(x.^2-l_lb.*x);
else
    f = @(x) a*(x.^3-(l_lb+k)*x.^2+k*l_lb*x);
end

l_ges = l_d2+l_lb+l_d1;

% background
delta_x = 0.3;
delta_y = 0.4;
rectangle('Position',[-delta_x -delta_y r_i+6.5*delta_x l_ges+2*delta_y],'FaceColor',[1 1 1],'LineStyle','none')
hold on

% liquid bridge
color.CTable = [color.liquid_light; color.liquid];
[sketch.x0,sketch.y0] = meshgrid(1:3,linspace(0,1,2));
[sketch.xn,sketch.yn] = meshgrid(1:3,linspace(0,1,n.colors));
cmap.liquid = interp2(sketch.x0,sketch.y0,color.CTable,sketch.xn,sketch.yn);
z_vec = fliplr(linspace(0,l_lb));
if r_c==0 && ax==1
    fill([0 r_i+f(z_vec) 0],l_d2+[0 fliplr(z_vec) l_lb],[1 zeros(size(z_vec)) 1])
    colormap(h.as.sketch,cmap.liquid)
else
    fill([0 r_i+f(z_vec) 0],l_d2+[0 fliplr(z_vec) l_lb],color.liquid_light)
end
plot(f(linspace(0,l_lb))+r_i,flip(linspace(0,l_lb))+l_d2,'lineWidth',w,'color',color.liquid)

% inclined wall lines
for i = 1:n.lines
    plot(r_i/n.lines*[i-1 i],l_ges+[0 r_i/n.lines],'lineWidth',w*0.25,'color','k')
    plot(r_i/n.lines*[i-1 i],[-r_i/n.lines 0],'lineWidth',w*0.25,'color','k')
end

% upper rod
if ax == 1
    color.CTable = [color.up_light; color.up];
    [sketch.x0,sketch.y0] = meshgrid(1:3,linspace(0,1,2));
    [sketch.xn,sketch.yn] = meshgrid(1:3,linspace(0,1,n.colors));
    cmap.up = interp2(sketch.x0,sketch.y0,color.CTable,sketch.xn,sketch.yn);
    rectangle('Position',[0 l_d2+l_lb r_i l_d1],'FaceColor',color.up_light,'LineStyle','none');
    for i = 1:n.colors
        rectangle('Position',[0 l_d2+l_lb (1-i/n.colors)*r_i l_d1],'FaceColor',cmap.up(i,:),'LineStyle','none')
    end
else
    rectangle('Position',[0 l_d2+l_lb r_i l_d1],'FaceColor',color.up_light,'LineStyle','none');
end

% lower rod
if ax == 1
    color.CTable = [color.down_light; color.down];
    [sketch.x0,sketch.y0] = meshgrid(1:3,linspace(0,1,2));
    [sketch.xn,sketch.yn] = meshgrid(1:3,linspace(0,1,n.colors));
    cmap.down = interp2(sketch.x0,sketch.y0,color.CTable,sketch.xn,sketch.yn);
    rectangle('Position',[0 0 r_i l_d2],'FaceColor',color.down_light,'LineStyle','none');
    for i = 1:n.colors
        rectangle('Position',[0 0 (1-i/n.colors)*r_i l_d2],'FaceColor',cmap.down(i,:),'LineStyle','none')
    end
else
    rectangle('Position',[0 0 r_i l_d2],'FaceColor',color.down_light,'LineStyle','none');
end

% dimensioning vectors for r_i & r_o
plot([0 r_i-3*d],vPos.ri*[1 1],'lineWidth',w,'color','k')
fill(r_i-3*d*[1 0 1],vPos.ri+d*[-1 0 1],'k','LineStyle','none')

% lines for the rods
plot([0 r_i],[0 0],'color',color.down,'lineWidth',w)
plot([0 r_i],l_d2*[1 1],'color',color.down,'lineWidth',w)
plot([0 r_i],l_d2+l_lb*[1 1],'color',color.up,'lineWidth',w)
plot([0 r_i],l_ges*[1 1],'color',color.up,'lineWidth',w)
plot(r_i*[1 1],[0 l_d2],'color',color.down,'lineWidth',w)
plot(r_i*[1 1],l_ges-[l_d1 0],'color',color.up,'lineWidth',w)

% central rod + dimensioning vector for r_c
if r_c ~= 0
    if ax == 1
        color.CTable = [color.up; cmap.up(round((1-r_C/r_i)*n.colors),:); 0.85*[1 1 1]; cmap.down(round((1-r_C/r_i)*n.colors),:); color.down];
    else
        color.CTable = [color.up_light; 0.85*[1 1 1]; color.down_light];
    end
    [sketch.x0,sketch.y0] = meshgrid(1:3,linspace(0,1,size(color.CTable,1)));
    [sketch.xn,sketch.yn] = meshgrid(1:3,linspace(0,1,n.colors));
    cmap.central = interp2(sketch.x0,sketch.y0,color.CTable,sketch.xn,sketch.yn);
    if ax == 1
        fill([0 r_C r_C r_C 0 0],l_d2+[-0.01*l_lb -0.01*l_lb l_lb/2 1.01*l_lb 1.01*l_lb l_lb/2],[1 0.75 0.5 0.25 0 0.5],'EdgeColor','none');
    else
        fill([0 r_C r_C r_C 0 0],l_d2+[-0.01*l_lb -0.01*l_lb l_lb/2 1.01*l_lb 1.01*l_lb l_lb/2],[1 1 0.5 0 0 0.5],'EdgeColor','none');
    end
    colormap(h.as.sketch,cmap.central)

    color.CTable = [color.up; 0.85*[1 1 1]; color.down];
    [sketch.x0,sketch.y0] = meshgrid(1:3,linspace(0,1,size(color.CTable,1)));
    [sketch.xn,sketch.yn] = meshgrid(1:3,linspace(0,1,n.colors));
    cmap.line = interp2(sketch.x0,sketch.y0,color.CTable,sketch.xn,sketch.yn);
    for i = 1:n.colors-1
        plot(r_C*[1 1],l_d2+z_vec(i:i+1),'Color',cmap.line(i,:),'LineWidth',w)
    end
    
    plot([0 r_C-3*d],l_d2+vPos.rc*[1 1],'color','k','lineWidth',w)
    fill(r_C-3*d*[1 0 1],l_d2+vPos.rc+d*[-1 0 1],'k','LineStyle','none')
end

% symmetry line
plot([0 0],[-r_i/n.lines l_ges+r_i/n.lines],'color','k','lineWidth',w*0.5,'lineStyle','-.')

% coordinate system
plot([0 l.orgn],l_d2+l_lb/2*[1 1],'lineWidth',1.2*w,'color','k')
plot([0 0],l_d2+l_lb/2+[0 l.orgn],'lineWidth',1.2*w,'color','k')
rectangle('Position',[-d/2 l_d2+l_lb/2-d/2 d d],'Curvature',[1 1],'FaceColor',[1 1 1],'lineWidth',1.2*w)
fill(l.orgn+[0 3*d 0],l_d2+l_lb/2+[-d 0 d],'w','lineWidth',1.2*w)
fill([-d 0 d],l_d2+l_lb/2+l.orgn+[0 3*d 0],'w','lineWidth',1.2*w)

% gravitiy vector
if g ~= 0
    plot(r_i+4*delta_x*[1 1],l_d2+l_lb+0.2*l_d1-[0 l.g],'lineWidth',1.2*w,'color','k')
    if g > 0
        fill(r_i+4*delta_x+d*[-1 0 1],l_d2+l_lb+0.2*l_d1-l.g-[0 3*d 0],'k','lineWidth',1.2*w)
    else
        fill(r_i+4*delta_x+d*[-1 0 1],l_d2+l_lb+0.2*l_d1+[0 3*d 0],'k','lineWidth',1.2*w)
    end
end

axis equal
set(gca,'visible','off')