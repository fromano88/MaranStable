if isempty(get(gcf,'Name'))
    clf
    % scaling
        unity = 'm'; dim = 1;
    % contour plot
        ct.draw = 1; ct.equal = 1;
        ct.domain   = [1 1]; % [1 0]-->LB, [0 1]-->SG, [1 1]-->all
        ct.quantity = 'T';
    % vector plot
        vp.draw = 0; vp.color = [0 0 0];
        vp.numR = 50; vp.numZ = 50;
    % isolines
        ct.streamlines = 1; ct.isolines = 0; ct.color.sf = [0 0 0];
        ct.numStreamlines = 8; ct.numIsolines = 20; ct.color.T = [0 0 0];

    % colormaps:
    %   myColormaps: Cold2Warm, BlueWhiteRed, thermal, thermal2
        ct.colormap = myColormap('Cool2Warm');
end

hold on

% needed to non-dimensionialize
    rho     = str2func(b1.rho);
    mu      = str2func(b1.mu);
    dsigma  = str2func(b1.dsigma);
    lambda  = str2func(b1.lambda);
    T_0     = (T_d1l+T_d2l)/2;
    delta_T = abs(T_d1l-T_d2l);

% coordinate shift
    z_shift = b1.geom.z(1)+0.5*diff(b1.geom.z);

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

% scales
    switch ct.quantity
        case 'T'
            scale = delta_T;
        case {'w','u','norm_u'}
            scale = abs(dsigma(T_0))*delta_T/mu(T_0);
        case 'p'
            scale = abs(dsigma(T_0))*delta_T/l_lb;
    end

% adjust domain for SFM
    if length(blocks) == 1
        ct.domain = [1 0];
    end

% transform ct.domain into 1, 2 or [1 2] if not already transformed by GUI
    if length(ct.domain)==2 && sum(ct.domain) < 3
        ct.domain = nonzeros(ct.domain.*[1 2])';
    end

for i = ct.domain
    b = eval(blocks{i});
    switch ct.quantity
        case 'T'
            if dim == 1
                ct.value = b.T.v;
                ct.label = '$T\; [^{\circ}$C]';
            else
                ct.value = (b.T.v-T_0)/scale;
                ct.label = '$\vartheta$';
            end
            
        case 'w'
            if dim == 1
                ct.value = -1./b.ETA_z.v.*b.w.v*f;
                ct.label = strcat('$w$ [',unity,'/s]');
            else
                ct.value = -1./b.ETA_z.v.*b.w.v/scale;
                ct.label = '$w$';
            end
            
        case 'u'
            if dim == 1
                ct.value = (1./b.XI_r.v.*b.u.v-b.JA.v.*b.XI_z.v.*b.w.v)*f;
                ct.label = strcat('$u$ [',unity,'/s]');
            else
                ct.value = (1./b.XI_r.v.*b.u.v-b.JA.v.*b.XI_z.v.*b.w.v)/scale;
                ct.label = '$u$';
            end
            
        case 'norm_u'
            if dim == 1
                ct.value = sqrt((1./b.ETA_z.v.*b.w.v).^2+(1./b.XI_r.v.*b.u.v-b.JA.v.*b.XI_z.v.*b.w.v).^2)*f;
                ct.label = strcat('$| \vec{u} |$ [',unity,'/s]');
            else
                ct.value = sqrt((1./b.ETA_z.v.*b.w.v).^2+(1./b.XI_r.v.*b.u.v-b.JA.v.*b.XI_z.v.*b.w.v).^2)/scale;
                ct.label = '$| \vec{u} |$';
            end
            
        case 'p'
            if isfield(b1,'dp_lb') && i==1
                p = b.p.v+b1.dp_lb;
            else
                p = b.p.v;
            end
            if dim == 1
                ct.value = p;
                ct.label = '$p$ [Pa]';
            else
                ct.value = p/scale;
                ct.label = '$p$';
            end
    end

    % ------------------------- contour plot ------------------------------
    if exist('ct.min','var')
        ct.min = min(ct.min,min(ct.value(:)));
        ct.max = max(ct.max,max(ct.value(:)));
    else
        ct.min = min(ct.value(:));
        ct.max = max(ct.value(:));
    end
    if ct.draw == 1
        [~,linie] = contourf(b.R.v*f,(z_shift-b.Z.v)*f,ct.value,100);
        set(linie,'Edgecolor','none');
    end
    plot(b.R.v(:,[1 end]) *f,(z_shift-b.Z.v(:,[1 end])) *f,'k');
    plot(b.R.v([1 end],:)'*f,(z_shift-b.Z.v([1 end],:)')*f,'k')
    
    % stream lines & temperature isolines
    if ct.isolines == 1
        [~, c] = contour(b.R.v*f,(z_shift-b.Z.v)*f,(b.T.v-T_0)/delta_T,'color',ct.color.T);
        c.LevelList = linspace(c.LevelList(1),c.LevelList(end),ct.numIsolines+1);
    end
    if flowopt.ax == 1
        % stoke's streamfunction (Nienhueser02 p.62)
        sf = full(-cumsum(b.R.w(2:end-1,1).*b.JA.w(2:end-1,1).*rho(b.T.w(2:end-1,1)).*b.w.w(2:end-1,1).*b.DXI.c(:,1))*ones(1,b.J) + cumsum(b.R.p.*b.JA.p.*rho(b.T.p).*b.u.p.*b.DETA.c,2));
    else
        sf = full(-cumsum(b.JA.w(2:end-1,1).*rho(b.T.w(2:end-1,1)).*b.w.w(2:end-1,1).*b.DXI.c(:,1))*ones(1,b.J) + cumsum(b.JA.p.*rho(b.T.p).*b.u.p.*b.DETA.c,2));
    end
    if ct.streamlines == 1
        [~, c] = contour(b.R.p*f,(z_shift-b.Z.p)*f,sf,'color',ct.color.sf);
        if ~isempty(c.LevelList)
            c.LevelList = linspace(c.LevelList(1),c.LevelList(end),ct.numStreamlines+1);
        end
    end
    % trick: if nothing should be plotted, then draw white (invisible) streamlines
    if ct.draw == 0 && ct.streamlines == 0 && ct.isolines == 0
        contour(b.R.p*f,(z_shift-b.Z.p)*f,sf,ct.numStreamlines,'w');
    end
    
    % -------------------------- vector plot ------------------------------
    if vp.draw == 1
        vp.n_r = round((b.geom.r(end)-b.geom.r(1))/(max(geom.r)-min(geom.r))*vp.numR);
        vp.n_z = round((b.geom.z(end)-b.geom.z(1))/(max(geom.z)-min(geom.z))*vp.numZ);
        vp.delta_r = (b.geom.r(end)-b.geom.r(1))/vp.n_r;
        vp.delta_z = (b.geom.z(end)-b.geom.z(1))/vp.n_z;
        [vp.Z, vp.R] = meshgrid(b.geom.z(1) + [0 cumsum(ones(1,vp.n_z)*vp.delta_z)],b.geom.r(1) + [0 cumsum(ones(1,vp.n_r)*vp.delta_r)]);
        
        vp.u = interp2(b.Z.u,b.R_cyl.u,1./b.XI_r.u.*b.u.u-b.JA.u.*b.XI_z.u.*b.w.u,vp.Z,vp.R, 'spline');
        vp.w = -interp2(b.Z.w,b.R_cyl.w,1./b.ETA_z.w.*b.w.w, vp.Z,vp.R,'spline');
        vp.R1   = ones(vp.n_r+1,1)*interp1(b.Z.u(1,:),b.R1.u(1,:),vp.Z(1,:),'spline');
        vp.Rend = ones(vp.n_r+1,1)*interp1(b.Z.u(end,:),b.Rend.u(end,:),vp.Z(end,:),'spline');
        
        R    = @(R1,Rend,R_cyl) R1+(R_cyl-b.geom.r(1)).*(Rend-R1)/(b.geom.r(end)-b.geom.r(1));
        vp.R = R(vp.R1,vp.Rend,vp.R);

        quiver(vp.R*f,(z_shift-vp.Z)*f,vp.u,vp.w,1,'LineWidth',1.2,'Color',vp.color)
    end
end

% colormap options
if ct.draw == 1 && ct.min~=ct.max
    ct.cb = colorbar;
    set(get(ct.cb,'ylabel'),'String',ct.label,'Interpreter','latex','FontSize',12)
    caxis([ct.min,ct.max])
    caxis manual
    colormap(ct.colormap)
end

% label options
p = flowopt.ax+1; rx = 'xr'; zy = 'yz';
if dim == 1
    xlabel(['$' rx(p) '$ [' unity ']'],'Interpreter','latex','FontSize',12)
    ylabel(['$' zy(p) '$ [' unity ']'],'Interpreter','latex','FontSize',12)
else
    xlabel(['$' rx(p) '$'],'Interpreter','latex','FontSize',12)
    ylabel(['$' zy(p) '$'],'Interpreter','latex','FontSize',12)
end


% figure options
box on
if ct.equal == 1
    axis('equal')
end

if ct.domain(end) == 2
    ylim([-l_lb/2-l_d2 l_lb/2+l_d1]*f)
else
    ylim([-l_lb/2 l_lb/2]*f)
end

if isempty(get(gcf,'Name'))
    set(gcf,'color','white');   
end

clear b
hold off