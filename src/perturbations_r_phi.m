if isempty(get(gcf,'Name'))
    clf
    % scaling
        unity = 'm'; dim = 1;
    % perturbation contour plot
        ct.draw = 1; ct.equal = 1;
        ct.domain   = [1 0]; % [1 0]-->LB, [0 1]-->SG, [1 1]-->all
        ct.quantity = 'v';
        ct.iso = 'th_E';
    % perturbation vector plot
        vp.draw = 1; vp.color = [0 0 0];
        vp.numR = 20; vp.numPhi = 10;
    % isoline
        ct.prod_pF = 0; ct.temp_pF = 0;
    % position
        ct.z = z_T_max;
    % colors and colormaps:
        ct.color.prod_pF = [0 0 0]; ct.color.T_pF = [0 0 0];
    %   myColormaps: Cool2Warm, Cold2Warm, BlueWhiteRed, thermal, thermal2
        ct.colormap = myColormap('BlueWhiteRed');
end

hold on

ax = flowopt.ax;
if ax == 1
    vp.numr = vp.numR*200;
else
    vp.numr = vp.numR*300;
end

time = time_T_max;
z    = ct.z; %z = l_d1/2%+l_lb/2

n_phi = 100;
if m == 0
    phi_T_max = 0;
end
if ax==0 && m~=0
    dphi = 2*pi/m/(n_phi+1);
    phi  = (0:dphi:2*pi/m);
else
    dphi = 2*pi/(n_phi+1);
    phi = phi_T_max+(0:dphi:2*pi);
end

% % if isempty(get(gcf,'Name')), clf, end; hold on

% dimensional vs dimensionless
if dim == 0
    f = 1/l_lb;
else
    f = 1000;
end

% adjust domain if z coordinate is outside LB
    if z < b1.geom.z(1) || z > b1.geom.z(end)
        ct.domain = [0 1];
    end

% adjust domain for SFM
    if length(blocks) == 1
        ct.domain = [1 0];
    end

% transform ct.domain into 1, 2 or [1 2] if not already transformed by GUI
    if length(ct.domain)==2 && sum(ct.domain) < 3
        ct.domain = nonzeros(ct.domain.*[1 2])';
    end

r=[];
for i = 1:length(blocks)
    b = eval(blocks{i});
    if (z<=b1.geom.z(1) || z>=b1.geom.z(end)) && i==1
        r = [r; interp2(b.Z.v,b.R_cyl.v,b.R.v, b.Z.v(1,1),b.geom.r(1)+[0; cumsum(ones(round(vp.numr*(b.geom.r(end) - b.geom.r(1))),1)*(b.geom.r(end) - b.geom.r(1))/round(vp.numr*(b.geom.r(end) - b.geom.r(1))))],'spline')*f];
    else
        r = [r; interp2(b.Z.v,b.R_cyl.v,b.R.v, z,b.geom.r(1)+[0; cumsum(ones(round(vp.numr*(b.geom.r(end) - b.geom.r(1))),1)*(b.geom.r(end) - b.geom.r(1))/round(vp.numr*(b.geom.r(end) - b.geom.r(1))))],'spline')*f];
    end
    r = sort(r);
end

max_v = []; ct.max = [];
for i = ct.domain
    b = eval(blocks{i});
    rho = str2func(b.rho);
    cp  = str2func(b.cp);
    R = interp2(b.Z.v,b.R_cyl.v,b.R.v, z,b.R_cyl.v(:,1),'spline')*f;

    p_hat.z = interp2(b.Z.v,b.R_cyl.v,b.p_hat.v, z,b.R_cyl.v(:,1),'spline');
    u_hat.z = interp2(b.Z.v,b.R_cyl.v,b.u_hat.v./b.XI_r.v - b.w_hat.v.*b.JA.v.*b.XI_z.v, z,b.R_cyl.v(:,1),'spline');
    v_hat.z = interp2(b.Z.v,b.R_cyl.v,b.v_hat.v, z,b.R_cyl.v(:,1),'spline');
    w_hat.z = interp2(b.Z.v,b.R_cyl.v,b.w_hat.v./b.ETA_z.v, z,b.R_cyl.v(:,1),'spline');
    U_hat.z = interp2(b.Z.v,b.R_cyl.v,b.u_hat.v, z,b.R_cyl.v(:,1),'spline');
    W_hat.z = interp2(b.Z.v,b.R_cyl.v,b.w_hat.v, z,b.R_cyl.v(:,1),'spline');

    if flowopt.energy == 1
        T_hat.z = interp2(b.Z.v,b.R_cyl.v,b.T_hat.v, z,b.R_cyl.v(:,1),'spline');

        cpT_xi.u  = (cp(b.T.T(2:end,:)).*b.T.T(2:end,:)-cp(b.T.T(1:end-1,:)).*b.T.T(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm b.DXI.rm(:,end)];
        cpT_eta.w = (cp(b.T.T(:,2:end)).*b.T.T(:,2:end)-cp(b.T.T(:,1:end-1)).*b.T.T(:,1:end-1))./[b.DETA.zm(1,:); b.DETA.zm; b.DETA.zm(end,:)];
        cpT_xi.v  = interp2(b.ETA.u,b.XI.u,cpT_xi.u,b.ETA.v,b.XI.v,'spline');
        cpT_eta.v = interp2(b.ETA.w,b.XI.w,cpT_eta.w,b.ETA.v,b.XI.v,'spline');
        rho_cpT_xi  = interp2(b.Z.v,b.R_cyl.v,rho(b.T.v).*cpT_xi.v,  z,b.R_cyl.v(:,1),'spline');
        rho_cpT_eta = interp2(b.Z.v,b.R_cyl.v,rho(b.T.v).*cpT_eta.v, z,b.R_cyl.v(:,1),'spline');
        th_prod = rho_cpT_xi.*plus_cc(T_hat.z,m,phi,gamma,time).*plus_cc(U_hat.z,m,phi,gamma,time);
        th_prod = -th_prod - rho_cpT_eta.*plus_cc(T_hat.z,m,phi,gamma,time).*plus_cc(W_hat.z,m,phi,gamma,time);
    end
    
    if z >= b.geom.z(1) && z <= b.geom.z(end)    
        switch ct.quantity
            case 'T'
                ct.value = plus_cc(T_hat.z,m,phi,gamma,time);
                ct.max   = T_max;
                ct.label = '$\vartheta''$';
                ct.ticklabel = {'cold','hot'};

            case 'u'
                ct.value = plus_cc(u_hat.z,m,phi,gamma,time);
                ct.max   = max([ct.max max(abs(ct.value(:)))]);
                ct.label = '$u''$';
                ct.ticklabel = {'low','high'};

            case 'v'
                ct.value = plus_cc(v_hat.z,m,phi,gamma,time);
                ct.max   = max([ct.max max(abs(ct.value(:)))]);
                ct.label = '$v''$';
                ct.ticklabel = {'low','high'};

            case 'w'
                ct.value = plus_cc(w_hat.z,m,phi,gamma,time);
                ct.max   = max([ct.max max(abs(ct.value(:)))]);
                ct.label = '$w''$';
                ct.ticklabel = {'low','high'};

            case 'p'
                ct.value = plus_cc(p_hat.z,m,phi,gamma,time);
                ct.max   = max([ct.max max(abs(ct.value(:)))]);
                ct.label = '$p''$';
                ct.ticklabel = {'low','high'};

            case 'th_E'
                ct.value  = th_prod;
                ct.max    = production_max;
                ct.label  = 'thermal energy production';
                ct.ticklabel = {'low','high'};
        end
        
        u_z = plus_cc(u_hat.z,m,phi,gamma,time);
        v_z = plus_cc(v_hat.z,m,phi,gamma,time);
        w_z = plus_cc(w_hat.z,m,phi,gamma,time);
        
        if flowopt.ax == 1
            [PHI,R] = meshgrid(-phi_T_max+phi,R);
            [X,Y]   = pol2cart(PHI,R);
        else
            [X,Y] = meshgrid(R,phi*l_lb*f);
        end
        
        if ct.draw == 1
            if ax == 1
                [~,linie] = contourf(X,Y,ct.value,256);
            else
                [~,linie] = contourf(X,Y,ct.value',256);
            end
            set(linie,'Edgecolor','none');
            max_v = max([max_v,max(ct.value(:))]);
        end
        if ct.prod_pF == 1
            if ax == 1
                [~, c] = contour(X,Y,th_prod,'color',ct.color.prod_pF);
            else
                [~, c] = contour(X,Y,th_prod','color',ct.color.prod_pF);
            end
            c.LevelList = c.LevelList(c.LevelList>0);
            c.LevelList = linspace(c.LevelList(1),c.LevelList(end),ct.iso_prod_pF+1);
        end
        if ct.temp_pF == 1
            if ax == 1
                [~, c] = contour(X,Y,plus_cc(T_hat.z,m,phi,gamma,time),'color',ct.color.T_pF);
            else
                [~, c] = contour(X,Y,plus_cc(T_hat.z,m,phi,gamma,time)','color',ct.color.T_pF);
            end
            c.LevelList = linspace(c.LevelList(1),c.LevelList(end),ct.iso_temp_pF+1);
        end
            
        if flowopt.ax == 1    
            line(R(1)*cos(0:2*pi/n_phi/100:2*pi)+0, R(1)*sin(0:2*pi/n_phi/100:2*pi)+0, 'color', 'black');
            line(R(end)*cos(0:2*pi/n_phi/100:2*pi)+0, R(end)*sin(0:2*pi/n_phi/100:2*pi)+0, 'color', 'black');
        else
            plot(X(:,[1 end]),Y(:,[1 end]),'k');
            plot(X([1 end],:)',Y([1 end],:)','k');
        end
       
        if vp.draw == 1
            if ax == 1
                ux = u_z.*cos(PHI)-v_z.*sin(PHI);
                uy = u_z.*sin(PHI)+v_z.*cos(PHI);
                warning off
                UX = scatteredInterpolant(reshape(X,[],1),reshape(Y,[],1),reshape(ux,[],1));
                UY = scatteredInterpolant(reshape(X,[],1),reshape(Y,[],1),reshape(uy,[],1));
                warning on
                [x,y] = meshgrid([-flipud(r); r(2:end)].',[-flipud(r); r(2:end)]);
                ux = UX(x,y);
                uy = UY(x,y);
                quiver(x(sqrt(x.^2+y.^2)>=R(1) & sqrt(x.^2+y.^2)<=R(end)),y(sqrt(x.^2+y.^2)>=R(1) & sqrt(x.^2+y.^2)<=R(end)),ux(sqrt(x.^2+y.^2)>=R(1) & sqrt(x.^2+y.^2)<=R(end)),uy(sqrt(x.^2+y.^2)>=R(1) & sqrt(x.^2+y.^2)<=R(end)),1,'Color',vp.color);

            else
                PHI = linspace(phi(1),phi(end),vp.numPhi);
                [x,y] = meshgrid(r(r>=R(1) & r<=R(end)),PHI*l_lb*f);
                ux = interp2(X,Y,u_z',x,y,'spline');
                uy = interp2(X,Y,v_z',x,y,'spline');
                quiver(x,y,ux,uy,1,'Color',vp.color);

            end
        end
        
    end
end

% colormap options
if ct.draw == 1
    colormap(ct.colormap)
    if max_v ~= 0
        caxis(max_v*[-1 1])
        ct.cb = colorbar('YTick', 0.8*max_v*[-1 1],'YTickLabel',ct.ticklabel);
    else
        caxis([0 1])
        ct.cb = colorbar('YTick', [0.1 0.9],'YTickLabel',ct.ticklabel);
    end
    
    set(get(ct.cb,'ylabel'),'String',ct.label,'Interpreter','latex','FontSize', 12)
end

p = ax+1; zy = 'zy';
if dim == 1
    xlabel('$x$ [mm]','Interpreter','latex','FontSize',12)
    ylabel(['$' zy(p) '$ [mm]'],'Interpreter','latex','FontSize',12)
else
    xlabel('$x$','Interpreter','latex','FontSize',12)
    ylabel(['$' zy(p) '$'],'Interpreter','latex','FontSize',12)
end

if imag(gamma(1,1))~=0 && m~=0 && ax==1
    if exist('r_o','var') && ct.domain(end) == 2
        quiver(0.8*r_o*f,0.8*r_o*f,0.1*r_o*f,-0.1*r_o*f,'Color','black','MaxHeadSize',10)
    else
        quiver(0.8*r_i*f,0.8*r_i*f,0.1*r_i*f,-0.1*r_i*f,'Color','black','MaxHeadSize',10)
    end
end

% figure options
if ct.equal == 1
    axis('equal')
end
box on
if isempty(get(gcf,'Name'))
    set(gcf,'color','white');   
end

clear b
if isfield(ct,'value')
    ct = rmfield(ct,'value');
end
hold off