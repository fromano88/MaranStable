if isempty(get(gcf,'Name'))
    clf
    % scaling
        dim = 0;
    % perturbation contour plot
        ct.draw = 1; ct.equal = 1;
        ct.domain   = [1 1]; % [1 0]-->LB, [0 1]-->SG, [1 1]-->all
        ct.quantity = 'v';
    % perturbation vector plot
        vp.draw = 0; vp.color = [0 0 0];
        vp.numR = 50; vp.numZ = 50;
    % isolines
        ct.temp_bS = 0; ct.sf_bS = 0; ct.prod_pF = 0; ct.temp_pF = 0;                 % flags for the isolines
        ct.iso_temp_bS = 8; ct.iso_sf_bS = 8; ct.iso_prod_pF = 8; ct.iso_temp_pF = 8; % number of isolines
    % position
        ct.phi = 0;
    % colors and colormaps:
        ct.color.T_bS = [0 0 0]; ct.color.sf_bS = [0 0 0]; ct.color.prod_pF = [0 0 0]; ct.color.T_pF = [0 0 0];
    %   myColormaps: Cool2Warm, Cold2Warm, BlueWhiteRed, thermal, thermal2
        ct.colormap = myColormap('BlueWhiteRed');
end

hold on

if m==0 || flowopt.ax==0
    phi = ct.phi;
else
    phi = phi_T_max+ct.phi*pi;
end
time = time_T_max;%0;%zeile/250*2*pi/imag(gamma(1,1));

if isempty(get(gcf,'Name')), clf, end; hold on

% needed to non-dimensionialize
    rho     = str2func(b1.rho);
    T_0     = (T_d1l+T_d2l)/2;
    delta_T = abs(T_d1l-T_d2l);

% coordinate shift
    z_shift = b1.geom.z(1)+0.5*diff(b1.geom.z);

% dimensional vs dimensionless
    if dim == 0
        f = 1/l_lb;
    else
        f = 1000; %mm
    end

% adjust domain for SFM
    if length(blocks) == 1
        ct.domain = [1 0];
    end

% transform ct.domain into 1, 2 or [1 2] if not already transformed by GUI
    if length(ct.domain)==2 && sum(ct.domain) < 3
        ct.domain = nonzeros(ct.domain.*[1 2])';
    end


if exist('ct','var') && isfield(ct,'value'), ct = rmfield(ct,'value'); end
ct.max = [];
for i = ct.domain
    b = eval(blocks{i});
    rho = str2func(b.rho);
    cp  = str2func(b.cp);
    
    if flowopt.energy == 1
        cpT_xi.u  = (cp(b.T.T(2:end,:)).*b.T.T(2:end,:)-cp(b.T.T(1:end-1,:)).*b.T.T(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm b.DXI.rm(:,end)];
        cpT_eta.w = (cp(b.T.T(:,2:end)).*b.T.T(:,2:end)-cp(b.T.T(:,1:end-1)).*b.T.T(:,1:end-1))./[b.DETA.zm(1,:); b.DETA.zm; b.DETA.zm(end,:)];
        cpT_xi.v  = interp2(b.ETA.u,b.XI.u,cpT_xi.u,b.ETA.v,b.XI.v,'spline');
        cpT_eta.v = interp2(b.ETA.w,b.XI.w,cpT_eta.w,b.ETA.v,b.XI.v,'spline');
        th_prod = rho(b.T.v).*cpT_xi.v.*plus_cc(b.T_hat.v,m,phi,gamma,time).*plus_cc(b.u_hat.v,m,phi,gamma,time);
        th_prod = -th_prod-rho(b.T.v).*cpT_eta.v.*plus_cc(b.T_hat.v,m,phi,gamma,time).*plus_cc(b.w_hat.v,m,phi,gamma,time);
    end
            
    switch ct.quantity
        case 'T'
            ct.value{i} = plus_cc(b.T_hat.v,m,phi,gamma,time);
            ct.max      = T_max;
            ct.label    = '$\vartheta''$';
            ct.ticklabel = {'cold','hot'};

        case 'u'
            ct.value{i} = 1./b.XI_r.v.*plus_cc(b.u_hat.v,m,phi,gamma,time)-b.JA.v.*b.XI_z.v.*plus_cc(b.w_hat.v,m,phi,gamma,time);
            ct.max      = max([ct.max max(abs(ct.value{i}(:)))]);
            ct.label    = '$u''$';
            ct.ticklabel = {'low','high'};

        case 'v'
            ct.value{i} = plus_cc(b.v_hat.v,m,phi,gamma,time);
            ct.max      = max([ct.max max(abs(ct.value{i}(:)))]);
            ct.label    = '$v''$';
            ct.ticklabel = {'low','high'};

        case 'w'
            ct.value{i} = 1./b.ETA_z.v.*plus_cc(b.w_hat.v,m,phi,gamma,time);
            ct.max      = max([ct.max max(abs(ct.value{i}(:)))]);
            ct.label    = '$w''$';
            ct.ticklabel = {'low','high'};

        case 'p'
            ct.value{i} = plus_cc(b.p_hat.v,m,phi,gamma,time);
            ct.max      = max([ct.max max(abs(ct.value{i}(:)))]);
            ct.label    = '$p''$';
            ct.ticklabel = {'low','high'};
            
        case 'th_E'
            ct.value{i} = th_prod;
            ct.max      = production_max;
            ct.label    = 'thermal energy production';
            ct.ticklabel = {'low','high'};
    end
    
    % ------------------------- contour plot ------------------------------    
    if ct.draw == 1
        [~,linie] = contourf(b.R.v*f,(z_shift-b.Z.v)*f,ct.value{i},250);
        set(linie,'Edgecolor','none');
    end
    plot(b.R.v(:,[1 end]) *f,(z_shift-b.Z.v(:,[1 end])) *f,'k');
    plot(b.R.v([1 end],:)'*f,(z_shift-b.Z.v([1 end],:)')*f,'k')
    
    % stream lines & temperature isolines
    if ct.temp_bS == 1
        [~, c] = contour(b.R.v*f,(z_shift-b.Z.v)*f,(b.T.v-T_0)/delta_T,'color',ct.color.T_bS);
        c.LevelList = linspace(c.LevelList(1),c.LevelList(end),ct.iso_temp_bS+1);
    end
    % stoke's streamfunction (Nienhueser02 p.62)
    if flowopt.ax == 1
        % stoke's streamfunction (Nienhueser02 p.62)
        sf = full(-cumsum(b.R.w(2:end-1,1).*b.JA.w(2:end-1,1).*rho(b.T.w(2:end-1,1)).*b.w.w(2:end-1,1).*b.DXI.c(:,1))*ones(1,b.J) + cumsum(b.R.p.*b.JA.p.*rho(b.T.p).*b.u.p.*b.DETA.c,2));
    else
        sf = full(-cumsum(b.JA.w(2:end-1,1).*rho(b.T.w(2:end-1,1)).*b.w.w(2:end-1,1).*b.DXI.c(:,1))*ones(1,b.J) + cumsum(b.JA.p.*rho(b.T.p).*b.u.p.*b.DETA.c,2));
    end
    if ct.sf_bS == 1
        [~, c] = contour(b.R.p*f,(z_shift-b.Z.p)*f,sf,'color',ct.color.sf_bS);
        if ~isempty(c.LevelList)
            c.LevelList = linspace(c.LevelList(1),c.LevelList(end),ct.iso_sf_bS+1);
        end
    end
    if ct.prod_pF == 1
        [~, c] = contour(b.R.v*f,(z_shift-b.Z.v)*f,th_prod,'color',ct.color.prod_pF);
        c.LevelList = c.LevelList(c.LevelList>=0); %take only positive values for production
        c.LevelList = linspace(c.LevelList(1),c.LevelList(end),ct.iso_prod_pF+1);
    end
    if ct.temp_pF == 1
        [~, c] = contour(b.R.v*f,(z_shift-b.Z.v)*f,plus_cc(b.T_hat.v,m,phi,gamma,time),'color',ct.color.T_pF);
        c.LevelList = linspace(c.LevelList(1),c.LevelList(end),ct.iso_temp_pF+1);
    end
    
    % -------------------------- vector plot ------------------------------
    if vp.draw == 1
        vp.n_r = round((b.geom.r(end)-b.geom.r(1))/(max(geom.r)-min(geom.r))*vp.numR);
        vp.n_z = round((b.geom.z(end)-b.geom.z(1))/(max(geom.z)-min(geom.z))*vp.numZ);
        vp.delta_r = (b.geom.r(end)-b.geom.r(1))/vp.n_r;
        vp.delta_z = (b.geom.z(end)-b.geom.z(1))/vp.n_z;
        [vp.Z, vp.R] = meshgrid(b.geom.z(1) + [0 cumsum(ones(1,vp.n_z)*vp.delta_z)],b.geom.r(1) + [0 cumsum(ones(1,vp.n_r)*vp.delta_r)]);
        
        vp.u =  interp2(b.Z.v,b.R_cyl.v,1./b.XI_r.v.*plus_cc(b.u_hat.v,m,phi,gamma,time)-b.JA.v.*b.XI_z.v.*plus_cc(b.w_hat.v,m,phi,gamma,time),vp.Z,vp.R,'spline');
        vp.w = -interp2(b.Z.v,b.R_cyl.v,1./b.ETA_z.v.*plus_cc(b.w_hat.v,m,phi,gamma,time),vp.Z,vp.R,'spline');
        vp.R1   = ones(vp.n_r+1,1)*interp1(b.Z.u(1,:),b.R1.u(1,:),vp.Z(1,:),'spline');
        vp.Rend = ones(vp.n_r+1,1)*interp1(b.Z.u(end,:),b.Rend.u(end,:),vp.Z(end,:),'spline');
        
        R    = @(R1,Rend,R_cyl) R1+(R_cyl-b.geom.r(1)).*(Rend-R1)/(b.geom.r(end)-b.geom.r(1));
        vp.R = R(vp.R1,vp.Rend,vp.R);

        q = quiver(vp.R*f,(z_shift-vp.Z)*f,vp.u,vp.w,1,'LineWidth',0.5,'Color',vp.color);
        q.MarkerSize  = 20;
        q.MaxHeadSize = 1;
    end
end

% colormap options
if ct.max ~= 0
    caxis([-ct.max ct.max])
end
caxis manual
colormap(ct.colormap)

% legend options
if ct.max ~= 0
    ct.cb = colorbar('YTick', 0.8*[-ct.max,ct.max],'YTickLabel',ct.ticklabel);
else
    ct.cb = colorbar('YTick', [0.1,0.9],'YTickLabel',ct.ticklabel);
end
set(get(ct.cb,'ylabel'),'String',ct.label,'Interpreter','latex','FontSize', 12)
p = flowopt.ax+1; rx = 'xr'; zy = 'yz';
if dim == 1
    xlabel(['$' rx(p) '$ [mm]'],'Interpreter','latex','FontSize',12)
    ylabel(['$' zy(p) '$ [mm]'],'Interpreter','latex','FontSize',12)
else
    xlabel(['$' rx(p) '$'],'Interpreter','latex','FontSize',12)
    ylabel(['$' zy(p) '$'],'Interpreter','latex','FontSize',12)
end

% figure options
if ct.equal == 1
    axis('equal')
end
box on
if isempty(get(gcf,'Name'))
    set(gcf,'color','white');   
end

h = b1.R.T(end,:); % free surface shape
if ct.domain(end) == 2
    ylim([-l_lb/2-l_d2 l_lb/2+l_d1]*f)
else
    ylim([-l_lb/2 l_lb/2]*f)
end

clear b
ct = rmfield(ct,'value');
hold off