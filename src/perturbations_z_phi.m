if isempty(get(gcf,'Name'))
    % scaling
        dim = 1;
    % perturbation contour plot
        ct.draw = 1; ct.equal = 1; ct.r = 0; % ct.r ... value between 0 and 1
    % perturbation vector plot
        vp.draw = 1; vp.color = [0 0 0]; ct.quantity = 'u';
        vp.numZ = 10; vp.numPhi = 60; %r25
    % isolines
        ct.temp_pF = 0; ct.prod_pF = 1;
        ct.iso_temp_pF = 8; ct.iso_prod_pF = 15;
    % colors and colormaps:
        ct.color.prod_pF = [0 0 0]; ct.color.T_pF = [0 0 0];
    %   myColormaps: Cool2Warm, Cold2Warm, BlueWhiteRed, thermal, thermal2
        ct.colormap = myColormap('BlueWhiteRed');
end

ax = flowopt.ax;
vp.numz = vp.numZ*1;

z_shift = b1.geom.z(1)+0.5*diff(b1.geom.z);

% dimensional vs dimensionless
if dim == 0
    f = 1/l_lb;
else
    f = 1000; %mm
end

time = time_T_max;%0;

n_phi = 100;

if ax==0
    if m~=0
        dphi = 2*pi/m/(n_phi+1);
        phi  = (0:dphi:2*pi/m);
    else
        dphi = 2*pi/(n_phi+1);
        phi  = (0:dphi:2*pi);
    end
else
    dphi = 2*pi/(n_phi+1);
    phi  = (0:dphi:2*pi)-pi;
end

if isempty(get(gcf,'Name')), clf, end; hold on

[PHI,Z]  = meshgrid(phi,(z_shift-b1.Z.v(1,:))*f);

if ct.r == 1
    idx = size(b1.R.v,1);
else
    idx = find(b1.R.v(:,1)>=ct.r*r_i,1);
end

if ax == 1
    r = b1.R.v(idx,1)*f;
else
    r = l_lb*f;
end
p =   b1.p_hat.v(idx,:).'*exp(-gamma(1,1)*time+1i*m*phi)+conj(b1.p_hat.v(idx,:).'*exp(-gamma(1,1)*time+1i*m*phi));
u =  (b1.u_hat.v(idx,:).'*exp(-gamma(1,1)*time+1i*m*phi)+conj(b1.u_hat.v(idx,:).'*exp(-gamma(1,1)*time+1i*m*phi))).*((1./b1.XI_r.v(idx,:).*sqrt(1+b1.R_z.v(idx,:).^2)).'*ones(size(phi))) - (b1.w_hat.v(idx,:).'*exp(-gamma(1,1)*time+1i*m*phi)+conj(b1.w_hat.v(idx,:).'*exp(-gamma(1,1)*time+1i*m*phi))).*((b1.JA.v(idx,:).*b1.XI_z.v(idx,:).*sqrt(1+b1.R_z.v(idx,:).^2)).'*ones(size(phi)));
v =   b1.v_hat.v(idx,:).'*exp(-gamma(1,1)*time+1i*m*phi)+conj(b1.v_hat.v(idx,:).'*exp(-gamma(1,1)*time+1i*m*phi));
w = -(b1.w_hat.v(idx,:).'*exp(-gamma(1,1)*time+1i*m*phi)+conj(b1.w_hat.v(idx,:).'*exp(-gamma(1,1)*time+1i*m*phi))).*((1./b1.ETA_z.v(idx,:).*sqrt(1+b1.R_z.v(idx,:).^2)).'*ones(size(phi)));

if flowopt.energy == 1
    T = plus_cc(b1.T_hat.v(idx,:).',m,phi,gamma,time);

    rho = str2func(b1.rho);
    cp  = str2func(b1.cp);
    cpT_xi.u  = (cp(b1.T.T(2:end,:)).*b1.T.T(2:end,:)-cp(b1.T.T(1:end-1,:)).*b1.T.T(1:end-1,:))./[b1.DXI.rm(:,1) b1.DXI.rm b1.DXI.rm(:,end)];
    cpT_eta.w = (cp(b1.T.T(:,2:end)).*b1.T.T(:,2:end)-cp(b1.T.T(:,1:end-1)).*b1.T.T(:,1:end-1))./[b1.DETA.zm(1,:); b1.DETA.zm; b1.DETA.zm(end,:)];
    cpT_xi.v  = interp2(b1.ETA.u,b1.XI.u,cpT_xi.u,b1.ETA.v,b1.XI.v,'spline');
    cpT_eta.v = interp2(b1.ETA.w,b1.XI.w,cpT_eta.w,b1.ETA.v,b1.XI.v,'spline');
    th_prod = rho(b1.T.v(idx,:))'.*cpT_xi.v(idx,:)'.*T.*plus_cc(b1.u_hat.v(idx,:).',m,phi,gamma,time);
    th_prod = -th_prod-rho(b1.T.v(idx,:))'.*cpT_eta.v(idx,:)'.*T.*plus_cc(b1.w_hat.v(idx,:).',m,phi,gamma,time);
end

if ct.draw == 1
    switch ct.quantity
        case 'T'
            [~,linie] = contourf((PHI)*r,Z,T,100);
            ct.label = '$\vartheta''$';
            ct.max = max(abs(T(:)));
            ct.ticklabel = {'cold','hot'};

        case 'u'
            [~,linie] = contourf((PHI)*r,Z,u,100);
            ct.label = '$u''$';
            ct.max = max(abs(u(:)));
            ct.ticklabel = {'low','high'};

        case 'v'
            [~,linie] = contourf((PHI)*r,Z,v,100);
            ct.label = '$v''$';
            ct.max = max(abs(v(:)));
            ct.ticklabel = {'low','high'};

        case 'w'
            [~,linie] = contourf((PHI)*r,Z,w,100);
            ct.label = '$w''$';
            ct.max = max(abs(w(:)));
            ct.ticklabel = {'low','high'};

        case 'p'
            [~,linie] = contourf((PHI)*r,Z,p,100);
            ct.label = '$p''$';
            ct.max = max(abs(p(:)));
            ct.ticklabel = {'low','high'};

        case 'th_E'
            [~,linie] = contourf((PHI)*r,Z,th_prod,100);
            ct.label = 'thermal energy production';
            ct.max = max(abs(th_prod(:)));
            ct.ticklabel = {'low','high'};
    end
    
    set(linie,'Edgecolor','none');
end
plot((PHI(:,[1 end]))*r,Z(:,[1 end]),'k')
plot((PHI([1 end],:)).'*r,Z([1 end],:).','k')

if ct.temp_pF == 1
    [~, c] = contour((PHI)*r,Z,T,'color',ct.color.T_pF);
    c.LevelList = linspace(c.LevelList(1),c.LevelList(end),ct.iso_temp_pF+1);
end

if ct.prod_pF == 1
    [~, c] = contour((PHI)*r,Z,th_prod,'color',ct.color.prod_pF);
    c.LevelList = c.LevelList(c.LevelList>=0); %take only positive values for production
    c.LevelList = linspace(c.LevelList(1),c.LevelList(end),ct.iso_prod_pF+1);
end

if vp.draw == 1
    vp.n_z = round((b1.geom.z(end) - b1.geom.z(1))/(max(geom.z)-min(geom.z))*vp.numz);
    vp.delta_z = (b1.geom.z(end) - b1.geom.z(1))/vp.n_z;
    [vp.PHI, vp.Z] = meshgrid(min(phi):2*pi/vp.numPhi:max(phi), (l_lb/2-[0 cumsum(ones(1,vp.n_z)*vp.delta_z)])*f);
    vp.w = interp2(PHI, Z, w, vp.PHI, vp.Z, 'spline');
    vp.v = interp2(PHI, Z, v, vp.PHI, vp.Z, 'spline');

    quiver((vp.PHI)*r,vp.Z,vp.v,vp.w,'Color',vp.color,'Linewidth',0.25,'MaxHeadSize',10);
    
end

if imag(gamma(1,1))~=0 && m~=0
    quiver(2*r_i*pi*f*0.02,1.2*l_lb/2*f,-r*pi/10,0,'Color',[1 1 1]*100/255,'MaxHeadSize',10)
end

if ct.draw == 1
    caxis(ct.max*[-1 1])
    caxis manual
    colormap(ct.colormap)
    
    ct.cb = colorbar('YTick', 0.8*ct.max*[-1 1],'YTickLabel',ct.ticklabel);
    set(get(ct.cb,'ylabel'),'String',ct.label,'Interpreter','latex','FontSize', 14)
end

yticks(l_lb/2*f*[-1 1])
if dim == 1
    unity = ' [mm]';
else
    unity = '';
end
if ax == 1
    xlim(r*pi*[-1 1])
    xticks(r*pi*[-1 -0.5 0 0.5 1]*f)
    xticklabels({'-$\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'})
    xlabel('$\varphi$','Interpreter','latex','FontSize', 12)
    ylabel(['$z$' unity],'Interpreter','latex','FontSize', 12)
else
    xlim(r*phi([1 end]))
    xlabel(['$z$' unity],'Interpreter','latex','FontSize', 12)
    ylabel(['$y$' unity],'Interpreter','latex','FontSize', 12)
end


% figure options
set(gca,'TickLabelInterpreter','latex')
if ct.equal == 1
    axis('equal')
end
if isempty(get(gcf,'Name'))
    set(gcf,'color','white');   
end

ylim(l_lb/2*[-1 1]*f)

box off
clear b
hold off