    w = ispc; % 1 --> windows
    p = h.flowopt.ax + 1; % pointer
    rx = 'xr'; zy = 'yz'; % coordinate: radial, axial
    mk = 'km';            % name of wavenumber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% geometry panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(h.gy.as.rc); cla(h.gy.as.rc);
    text(h.tp_l,h.tp_y,['$' rx(p) '_c=$'],h.latex_l{:})
    text(h.tp_r,h.tp_y,'mm',h.latex_r{:})
    % ------------------------------------------------------------------- %
    axes(h.gy.as.ri); cla(h.gy.as.ri);
    text(h.tp_l,h.tp_y,['$' rx(p) '_i=$'],h.latex_l{:})
    text(h.tp_r,h.tp_y,'mm',h.latex_r{:})
    % ------------------------------------------------------------------- %
    axes(h.gy.as.ro); cla(h.gy.as.ro);
    text(h.tp_l,h.tp_y,['$' rx(p) '_o=$'],h.latex_l{:})
    text(h.tp_r,h.tp_y,'mm',h.latex_r{:})
    % ------------------------------------------------------------------- %
    if h.flowopt.ax == 1
        set(h.gy.pl.radii,'title','Radii')
    else
        set(h.gy.pl.radii,'title','Widths')
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mesh panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if h.flowopt.ax == 1
        set(h.mh.pm.headline,'string',{'radial';'axial'})
    else
        set(h.mh.pm.headline,'string',{'horizontal';'vertical'})
    end
% --------------------------------------------------------------------- %
% ------------------------radial-Liquid-Bridge------------------------- %
% --------------------------------------------------------------------- %
    for i = 1:3
        axes(h.mh.as.radialLB.spacing(i)); cla(h.mh.as.radialLB.spacing(i))
        if i == 2
            text(h.tp_l,h.tp_y,'$\Delta_\mathrm{max}=$',h.latex_l{:})
        else
            text(h.tp_l,h.tp_y,'$\Delta_\mathrm{min}=$',h.latex_l{:})
        end
        text(h.tp_r,h.tp_y,'$\%\,$ of $\,d$',h.latex_r{:})
        %---------------------------------------------------------------%
        axes(h.mh.as.radialLB.stretchingC(i)); cla(h.mh.as.radialLB.stretchingC(i))
        text(h.tp_l,h.tp_y,'$f=$',h.latex_l{:})
        %---------------------------------------------------------------%
        cla(h.mh.as.radialLB.stretchingF(i))
    end
    axes(h.mh.as.radialLB.spacing(3));
    text(-5.2,h.tp_y,['At $\,' rx(p) '=' rx(p) '_i$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.radialLB.stretchingC(3));
    text(-5.2,h.tp_y,['At $\,' rx(p) ' \cong ' rx(p) '_i$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.radialLB.stretchingF(3));
    text(0.3,0.52,['At $\,' rx(p) ' \cong ' rx(p) '_i$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    if h.r_c == 0
        axes(h.mh.as.radialLB.spacing(1));
        text(-5.2,h.tp_y,['At $\,' rx(p) '=0$:'],h.latex_r{:})
        % --------------------------------------------------------------- %
        axes(h.mh.as.radialLB.spacing(2));
        text(-5.2,h.tp_y,['At $\,' rx(p) ' = ' rx(p) '_i/2$:'],h.latex_r{:})
        % --------------------------------------------------------------- %
        axes(h.mh.as.radialLB.stretchingC(1));
        text(-5.2,h.tp_y,['At $\,' rx(p) ' \cong 0$:'],h.latex_r{:})
        % --------------------------------------------------------------- %
        axes(h.mh.as.radialLB.stretchingC(2));
        text(-5.2,h.tp_y,['At $\,' rx(p) ' \cong ' rx(p) '_i/2$:'],h.latex_r{:})
        % --------------------------------------------------------------- %
        axes(h.mh.as.radialLB.stretchingF(1));
        text(0.3,0.52,['At $\,' rx(p) ' \cong 0$:'],h.latex_r{:})
        % --------------------------------------------------------------- %
        axes(h.mh.as.radialLB.stretchingF(2));
        text(0.3,0.52,['At $\,' rx(p) ' \cong ' rx(p) '_i/2$:'],h.latex_r{:})
    else
        axes(h.mh.as.radialLB.spacing(1));
        text(-5.2,h.tp_y,['At $\,' rx(p) '=' rx(p) '_c$:'],h.latex_r{:})
        % --------------------------------------------------------------- %
        axes(h.mh.as.radialLB.spacing(2));
        text(-5.2,h.tp_y,['At $\,' rx(p) ' = (' rx(p) '_c+' rx(p) '_i)/2$:'],h.latex_r{:})
        % --------------------------------------------------------------- %
        axes(h.mh.as.radialLB.stretchingC(1));
        text(-5.2,h.tp_y,['At $\,' rx(p) ' \cong ' rx(p) '_c$:'],h.latex_r{:})
        % --------------------------------------------------------------- %
        axes(h.mh.as.radialLB.stretchingC(2));
        text(-5.2,h.tp_y,['At $\,' rx(p) ' \cong (' rx(p) '_c+' rx(p) '_i)/2$:'],h.latex_r{:})
        % --------------------------------------------------------------- %
        axes(h.mh.as.radialLB.stretchingF(1));
        text(0.3,0.52,['At $\,' rx(p) ' \cong ' rx(p) '_c$:'],h.latex_r{:})
        % --------------------------------------------------------------- %
        axes(h.mh.as.radialLB.stretchingF(2));
        text(0.3,0.52,['At $\,' rx(p) ' \cong (' rx(p) '_c+' rx(p) '_i)/2$:'],h.latex_r{:})
    end
% --------------------------------------------------------------------- %
% -------------------------axial-Liquid-Bridge------------------------- %
% --------------------------------------------------------------------- %
    for i = 1:3
        axes(h.mh.as.axialLB.spacing(i)); cla(h.mh.as.axialLB.spacing(i))
        if i == 2
            text(h.tp_l,h.tp_y,'$\Delta_\mathrm{max}=$',h.latex_l{:})
        else
            text(h.tp_l,h.tp_y,'$\Delta_\mathrm{min}=$',h.latex_l{:})
        end
        text(h.tp_r,h.tp_y,'$\%\,$ of $\,d$',h.latex_r{:})
        %---------------------------------------------------------------%
        axes(h.mh.as.axialLB.stretchingC(i)); cla(h.mh.as.axialLB.stretchingC(i))
        text(h.tp_l,h.tp_y,'$f=$',h.latex_l{:})
        %---------------------------------------------------------------%
        cla(h.mh.as.axialLB.stretchingF(i))
    end
    axes(h.mh.as.axialLB.spacing(1));
    text(-5.2,h.tp_y,['At $\,' zy(p) '=d/2$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.axialLB.spacing(2));
    text(-5.2,h.tp_y,['At $\,' zy(p) ' = 0$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.axialLB.spacing(3));
    text(-5.2,h.tp_y,['At $\,' zy(p) ' = -d/2$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.axialLB.stretchingC(1));
    text(-5.2,h.tp_y,['At $\,' zy(p) ' \cong d/2$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.axialLB.stretchingC(2));
    text(-5.2,h.tp_y,['At $\,' zy(p) ' \cong 0$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.axialLB.stretchingC(3));
    text(-5.2,h.tp_y,['At $\,' zy(p) ' \cong -d/2$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.axialLB.stretchingF(1));
    text(0.3,0.52,['At $\,' zy(p) ' \cong d/2$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.axialLB.stretchingF(2));
    text(0.3,0.52,['At $\,' zy(p) ' \cong 0$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.axialLB.stretchingF(3));
    text(0.3,0.52,['At $\,' zy(p) ' \cong -d/2$:'],h.latex_r{:})
% --------------------------------------------------------------------- %
% -----------------------radial-Surrounding-Gas------------------------ %
% --------------------------------------------------------------------- %
    for i = 1:3
        axes(h.mh.as.radialSG.spacing(i)); cla(h.mh.as.radialSG.spacing(i))
        if i == 2
            text(h.tp_l,h.tp_y,'$\Delta_\mathrm{max}=$',h.latex_l{:})
        else
            text(h.tp_l,h.tp_y,'$\Delta_\mathrm{min}=$',h.latex_l{:})
        end
        %---------------------------------------------------------------%
        axes(h.mh.as.radialSG.stretchingC(i)); cla(h.mh.as.radialSG.stretchingC(i))
        text(h.tp_l,h.tp_y,'$f=$',h.latex_l{:})
        %---------------------------------------------------------------%
        cla(h.mh.as.radialSG.stretchingF(i))
    end
    axes(h.mh.as.radialSG.spacing(1));
    text(-5.2,h.tp_y,['At $\,' rx(p) '=' rx(p) '_i$:'],h.latex_r{:})
    text(h.tp_r,h.tp_y,'$\%\,$ of $\,d$',h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.radialSG.spacing(2));
    text(-5.2,h.tp_y,['At $\,' rx(p) ' = (' rx(p) '_i+' rx(p) '_o)/2$:'],h.latex_r{:})
    text(h.tp_r,h.tp_y,['$\%\,$ of $' rx(p) '_o-' rx(p) '_i$'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.radialSG.spacing(3));
    text(-5.2,h.tp_y,['At $\,' rx(p) '=' rx(p) '_o$:'],h.latex_r{:})
    text(h.tp_r,h.tp_y,['$\%\,$ of $' rx(p) '_o-' rx(p) '_i$'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.radialSG.stretchingC(1));
    text(-5.2,h.tp_y,['At $\,' rx(p) ' \cong ' rx(p) '_i$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.radialSG.stretchingC(2));
    text(-5.2,h.tp_y,['At $\,' rx(p) ' \cong (' rx(p) '_i+' rx(p) '_o)/2$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.radialSG.stretchingC(3));
    text(-5.2,h.tp_y,['At $\,' rx(p) ' \cong ' rx(p) '_o$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.radialSG.stretchingF(1));
    text(0.3,0.52,['At $\,' rx(p) ' \cong ' rx(p) '_i$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.radialSG.stretchingF(2));
    text(0.3,0.52,['At $\,' rx(p) ' \cong (' rx(p) '_i+' rx(p) '_o)/2$:'],h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.radialSG.stretchingF(3));
    text(0.3,0.52,['At $\,' rx(p) ' \cong ' rx(p) '_o$:'],h.latex_r{:})
% --------------------------------------------------------------------- %
% ------------------------axial-Surrounding-Gas------------------------ %
% --------------------------------------------------------------------- %
    for i = 1:7
        axes(h.mh.as.axialSG.spacing(i)); cla(h.mh.as.axialSG.spacing(i))
        if mod(i,2) == 0
            text(h.tp_l,h.tp_y,'$\Delta_\mathrm{max}=$',h.latex_l{:})
        else
            text(h.tp_l,h.tp_y,'$\Delta_\mathrm{min}=$',h.latex_l{:})
        end
        %---------------------------------------------------------------%
        axes(h.mh.as.axialSG.stretchingC(i)); cla(h.mh.as.axialSG.stretchingC(i))
        text(h.tp_l,h.tp_y,'$f=$',h.latex_l{:})
        %---------------------------------------------------------------%
        cla(h.mh.as.axialSG.stretchingF(i))
    end
    axes(h.mh.as.axialSG.spacing(1));
    text(-3.54,h.tp_y,['At $\,' zy(p) '=d/2+d_1$:'],h.latex_r{:})
    text(h.tp_r,h.tp_y,'$\%\,$ of $\,d_1$',h.latex_r{:})
    % --------------------------------------------------------------- %
    axes(h.mh.as.axialSG.spacing(2));
    text(-3.54,h.tp_y,['At $\,' zy(p) '= d/2+d_1/2$:'],h.latex_r{:})
    text(h.tp_r,h.tp_y,'$\%\,$ of $\,d_1$',h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.spacing(3));
    text(-3.54,h.tp_y,['At $\,' zy(p) '=d/2$:'],h.latex_r{:})
    text(h.tp_r,h.tp_y,'$\%\,$ of $\,d$',h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.spacing(4));
    text(-3.54,h.tp_y,['At $\,' zy(p) '= 0$:'],h.latex_r{:})
    text(h.tp_r,0.42,'$\%\,$ of $\,d$',h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.spacing(5));
    text(-3.46,h.tp_y,['At $\,' zy(p) '= -d/2$:'],h.latex_r{:})
    text(h.tp_r,h.tp_y,'$\%\,$ of $\,d$',h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.spacing(6));
    text(-3.46,h.tp_y,['At $\,' zy(p) '= -d/2-d_2/2$:'],h.latex_r{:})
    text(h.tp_r,h.tp_y,'$\%\,$ of $\,d_2$',h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.spacing(7));
    text(-3.46,h.tp_y,['At $\,' zy(p) '= -d/2-d_2$:'],h.latex_r{:})
    text(h.tp_r,h.tp_y,'$\%\,$ of $\,d_2$',h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingC(1));
    text(-3.54,h.tp_y,['At $\,' zy(p) ' \cong d/2+d_1$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingC(2));
    text(-3.54,h.tp_y,['At $\,' zy(p) ' \cong d/2+d_1/2$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingC(3));
    text(-3.54,h.tp_y,['At $\,' zy(p) ' \cong d/2$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingC(4));
    text(-3.54,h.tp_y,['At $\,' zy(p) ' \cong 0$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingC(5));
    text(-3.46,h.tp_y,['At $\,' zy(p) ' \cong -d/2$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingC(6));
    text(-3.46,h.tp_y,['At $\,' zy(p) ' \cong -d/2-d_2/2$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingC(7));
    text(-3.46,h.tp_y,['At $\,' zy(p) ' \cong -d/2-d_2$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingF(1));
    text(-0.05,0.52,['At $\,' zy(p) ' \cong d/2+d_1$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingF(2));
    text(-0.05,0.52,['At $\,' zy(p) ' \cong d/2+d_1/2$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingF(3));
    text(-0.05,0.52,['At $\,' zy(p) ' \cong d/2$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingF(4));
    text(-0.05,0.52,['At $\,' zy(p) ' \cong 0$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingF(5));
    text(-0.05,0.52,['At $\,' zy(p) ' \cong -d/2$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingF(6));
    text(-0.05,0.52,['At $\,' zy(p) ' \cong -d/2-d_2/2$:'],h.latex_r{:})
    % ------------------------------------------------------------------%
    axes(h.mh.as.axialSG.stretchingF(7));
    text(-0.05,0.52,['At $\,' zy(p) ' \cong -d/2-d_2$:'],h.latex_r{:})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% equations panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(h.es.as.continuity); cla(h.es.as.continuity);
    if h.NS==1 || get(h.es.cb.energy,'Value')==0
        text(0,0,'$\nabla \cdot \vec u=0$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
    elseif h.NS == 2
        text(0,0,'$-\rho_0\alpha_{\rho}\partial_t T+\nabla \cdot (\rho \vec u)=0$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
    else
        text(0,0,'$\partial_t \rho+\nabla \cdot (\rho \vec u)=0$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
    end
    % ------------------------------------------------------------------%
    axes(h.es.as.momentum); cla(h.es.as.momentum);
    if get(h.es.cb.energy,'Value')==0
        if get(h.es.cb.stokes,'Value')
            text(0,0,'$\partial_t \vec u = -\frac{1}{\rho_0} \nabla  p+\frac{\mu_0}{\rho_0} \nabla^2 \vec u$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
        else
            text(0,0,'$\partial_t \vec u + \vec u\cdot\nabla  \vec u = -\frac{1}{\rho_0} \nabla  p+\frac{\mu_0}{\rho_0} \nabla^2 \vec u$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
        end
    else
        if h.NS == 1
            if get(h.es.cb.stokes,'Value')
                text(0,0,'$\partial_t \vec u = -\frac{1}{\rho_0} \nabla  p+\frac{\mu_0}{\rho_0} \nabla^2 \vec u -\vec{g}\alpha_\rho (T-T_0)$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
            else
                text(0,0,'$\partial_t \vec u + \vec u\cdot\nabla\vec u = -\frac{1}{\rho_0} \nabla  p+\frac{\mu_0}{\rho_0} \nabla^2 \vec u -\vec{g}\alpha_\rho (T-T_0)$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
            end
        elseif h.NS == 2
            if get(h.es.cb.stokes,'Value')
                text(0,0,'$\rho \partial_t \vec u  = -\nabla  p+ \nabla \cdot {\mu \boldmath \tau}-\vec{g}\alpha_\rho (T-T_0)$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
            else
                text(0,0,'$\rho(\partial_t \vec u + \vec u\cdot\nabla\vec u) = -\nabla  p+ \nabla \cdot (\mu{ \boldmath \tau})-\vec{g}\alpha_\rho (T-T_0)$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
            end
        else
            if get(h.es.cb.stokes,'Value')
                text(0,0,'$\partial_t(\rho \vec u)  = -\nabla  p + \nabla \cdot (\mu{ \boldmath \tau})+\rho\vec{g}$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
            else
                text(0,0,'$\partial_t(\rho \vec u) + \nabla\cdot(\rho \vec u \vec u) = -\nabla  p + \nabla \cdot (\mu{ \boldmath \tau})+\rho\vec{g}$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
            end
        end
    end
    % ------------------------------------------------------------------%
    axes(h.es.as.energy); cla(h.es.as.energy);
    if h.NS==1 || get(h.es.cb.energy,'Value')==0
        tt_energy = text(0,0,'$\partial_t T+\vec u\cdot\nabla T =\frac{\lambda_0}{\rho_0 c_{p0}}\nabla^2 T$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center');
        if get(h.es.cb.energy,'value') == 0
            tt_energy.Color = 0.7*[1 1 1];
        end
    elseif h.NS==2
        text(0,0,'$\rho c_p(\partial_t T+\vec u\cdot\nabla T)=\lambda\nabla^2 T+\lambda_0\alpha_\lambda(\nabla T)^2$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center');
    else
        text(0,0,'$\rho c_p(\partial_t T+\vec u\cdot\nabla T)=\nabla \cdot (\lambda \nabla T)$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center');
    end
    %-------------------------------------------------------------------%
    if ispc
        strt = 4.2;
    else
        strt = 7;
    end
    dlt = 1;
    axes(h.es.as.properties); cla(h.es.as.properties);
    if h.NS == 1
        text(0,strt+4*dlt,'$\rho_0$',h.latex_r{:}); text(1-0.28*w,strt+4*dlt,'$\dots\;$ const.\ reference density',h.latex_r{:})
        text(0,strt+3*dlt,'$\mu_0$',h.latex_r{:});  text(1-0.28*w,strt+3*dlt,'$\dots\;$ const.\ dynamic viscosity',h.latex_r{:})
        if get(h.es.cb.energy,'Value')==1
            text(0,strt+2*dlt,'$c_{p0}$',h.latex_r{:});    text(1-0.28*w,strt+2*dlt,'$\dots\;$ const.\ specific heat capacity',h.latex_r{:})
            text(0,strt+1*dlt,'$\lambda_0$',h.latex_r{:}); text(1-0.28*w,strt+1*dlt,'$\dots\;$ const.\ thermal conductivity',h.latex_r{:})
            text(0,strt,'$\alpha_\rho$',h.latex_r{:});     text(1-0.28*w,strt,'$\dots\;$ const.\ thermal expansion coefficient',h.latex_r{:})
        end
    elseif h.NS == 2
        text(0,strt,'${\boldmath \tau} = \nabla \vec u + (\nabla \vec u)^\mathrm{T} -\frac{2}{3} (\nabla \cdot \vec u ){\boldmath I}$',h.latex_r{:})
        text(0,strt+1*dlt,'$\lambda=\lambda_0[1-\alpha_\lambda(T-T_0)]$',h.latex_r{:})
        text(0,strt+2*dlt,'$c_p=c_{p0}[1-\alpha_{cp}(T-T_0)]$',h.latex_r{:})
        text(0,strt+3*dlt,'$\mu=\mu_0[1-\alpha_\mu(T-T_0)]$',h.latex_r{:})
        text(0,strt+4*dlt,'$\rho=\rho_0[1-\alpha_\rho(T-T_0)]$',h.latex_r{:})
        text(5.85-0.28*w,strt,'$\dots\;$ stress tensor',h.latex_r{:})
        text(5.85-0.28*w,strt+1*dlt,'$\dots\;$ thermal conductivity',h.latex_r{:})
        text(5.85-0.28*w,strt+2*dlt,'$\dots\;$ specific heat capacity',h.latex_r{:})
        text(5.85-0.28*w,strt+3*dlt,'$\dots\;$ dynamic viscosity',h.latex_r{:})
        text(5.85-0.28*w,strt+4*dlt,'$\dots\;$ density',h.latex_r{:})
    else
        text(0,strt,'${\boldmath \tau} = \nabla \vec u + (\nabla \vec u)^\mathrm{T} -\frac{2}{3} (\nabla \cdot \vec u ){\boldmath I}$',h.latex_r{:})
        text(0,strt+1*dlt,'$\lambda=\lambda(T)$',h.latex_r{:})
        text(0,strt+2*dlt,'$c_p=c_p(T)$',h.latex_r{:})
        text(0,strt+3*dlt,'$\mu=\mu(T)$',h.latex_r{:})
        text(0,strt+4*dlt,'$\rho=\rho(T)$',h.latex_r{:})
        text(5.85-0.28*w,strt,'$\dots\;$ stress tensor',h.latex_r{:})
        text(5.85-0.28*w,strt+1*dlt,'$\dots\;$ thermal conductivity',h.latex_r{:})
        text(5.85-0.28*w,strt+2*dlt,'$\dots\;$ specific heat capacity',h.latex_r{:})
        text(5.85-0.28*w,strt+3*dlt,'$\dots\;$ dynamic viscosity',h.latex_r{:})
        text(5.85-0.28*w,strt+4*dlt,'$\dots\;$ density',h.latex_r{:})
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% boundary conditions panel %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if h.flowopt.ax == 1
        set(h.bc.pm.headline,'string',{'inner/outer','top/bottom'})
    else
        set(h.bc.pm.headline,'string',{'left/right','top/bottom'})
    end
% --------------------------------------------------------------------- %
% -------------------------radial-Liquid-Bridge------------------------ %
% --------------------------------------------------------------------- %
    axes(h.bc.as.radialLB.rc(1)); cla(h.bc.as.radialLB.rc(1));
    if h.r_c == 0
        text(0.4,0.52,['At $\,' rx(p) ' = 0: \, \, \, \, -d/2 \leq ' zy(p) ' \leq d/2$'],h.latex_r{:})
    else
        text(0.4,0.52,['At $\,' rx(p) ' = ' rx(p) '_c: \, \, \, \, -d/2 \leq ' zy(p) ' \leq d/2$'],h.latex_r{:})
    end
        if get(h.es.cb.energy,'value') && strcmp(h.b1.bc.z(1,3),'c')
            axes(h.bc.as.radialLB.rc(2)); cla(h.bc.as.radialLB.rc(2));
            text(0.94*h.tp_r,h.tp_y,'$^\circ$C',h.latex_r{:})
            if get(h.bc.pm.radialLB.rc.temperature,'value') == 1
                text(0.4*h.tp_l,h.tp_y,['$T_{' rx(p) '_c} =$'],h.latex_l{:})
            else
                text(0.4*h.tp_l,h.tp_y,['$T_{' rx(p) '_c}(' zy(p) ') =$'],h.latex_l{:})
            end
            %-------------------------------------------------------------------%
            displayTemperature = get(h.bc.et.radialLB.rc.temperature,'String');
            idx = strfind(displayTemperature,zy(3-p));
            if ~isempty(idx)
                displayTemperature(idx) = zy(p);
                set(h.bc.et.radialLB.rc.temperature,'String',displayTemperature);
            end
        end
        %-------------------------------------------------------------------%
    axes(h.bc.as.radialLB.ri); cla(h.bc.as.radialLB.ri);
    text(0.4,0.52,['At $\,' rx(p) ' = h(' zy(p) '): \, \, \, \, -d/2 \leq ' zy(p) ' \leq d/2$'],h.latex_r{:})
% --------------------------------------------------------------------- %
% --------------------------axial-Liquid-Bridge------------------------ %
% --------------------------------------------------------------------- %
    axes(h.bc.as.axialLB.d1(1)); cla(h.bc.as.axialLB.d1(1));
    if h.r_c == 0
        text(0.4,0.52,['At $\,' zy(p) ' = d/2: \, \, \, \, 0 \leq ' rx(p) ' \leq ' rx(p) '_i$'],h.latex_r{:})
    else
        text(0.4,0.52,['At $\,' zy(p) ' = d/2: \, \, \, \, ' rx(p) '_c \leq ' rx(p) ' \leq ' rx(p) '_i$'],h.latex_r{:})
    end
    %-------------------------------------------------------------------%
    if get(h.es.cb.energy,'value') && strcmp(h.b1.bc.r(1,3),'c')
        axes(h.bc.as.axialLB.d1(2)); cla(h.bc.as.axialLB.d1(2));
        text(0.94*h.tp_r,h.tp_y,'$^\circ$C',h.latex_r{:})
        if get(h.bc.pm.axialLB.d1.temperature,'value') == 2
            text(0.4*h.tp_l,h.tp_y,['$T_{d1}(' rx(p) ') =$'],h.latex_l{:})
        else
            text(0.4*h.tp_l,h.tp_y,'$T_{d1} =$',h.latex_l{:})
        end
        %-------------------------------------------------------------------%
        displayTemperature = get(h.bc.et.axialLB.d1.temperature,'String');
        idx = strfind(displayTemperature,rx(3-p));
        if ~isempty(idx)
            displayTemperature(idx) = rx(p);
            set(h.bc.et.axialLB.d1.temperature,'String',displayTemperature);
        end
    end
    %-------------------------------------------------------------------%
    axes(h.bc.as.axialLB.d2(1)); cla(h.bc.as.axialLB.d2(1));
    if h.r_c == 0
        text(0.4,0.52,['At $\,' zy(p) ' = -d/2: \, \, \, \, 0 \leq ' rx(p) ' \leq ' rx(p) '_i$'],h.latex_r{:})
    else
        text(0.4,0.52,['At $\,' zy(p) ' = -d/2: \, \, \, \, ' rx(p) '_c \leq ' rx(p) ' \leq ' rx(p) '_i$'],h.latex_r{:})
    end
    %-------------------------------------------------------------------%
    if get(h.es.cb.energy,'value') && strcmp(h.b1.bc.r(2,3),'c')
        axes(h.bc.as.axialLB.d2(2)); cla(h.bc.as.axialLB.d2(2));
        text(0.94*h.tp_r,h.tp_y,'$^\circ$C',h.latex_r{:})
        if get(h.bc.pm.axialLB.d2.temperature,'value') == 2
            text(0.4*h.tp_l,h.tp_y,['$T_{d2}(' rx(p) ') =$'],h.latex_l{:})
        else
            text(0.4*h.tp_l,h.tp_y,'$T_{d2} =$',h.latex_l{:})
        end
        %-------------------------------------------------------------------%
        displayTemperature = get(h.bc.et.axialLB.d2.temperature,'String');
        idx = strfind(displayTemperature,rx(3-p));
        if ~isempty(idx)
            displayTemperature(idx) = rx(p);
            set(h.bc.et.axialLB.d2.temperature,'String',displayTemperature);
        end
    end
% ---------------------------------------------------------------------- %
% -------------------------radial-Surrounding-Gas----------------------- %
% ---------------------------------------------------------------------- %
    axes(h.bc.as.radialSG.d1(1)); cla(h.bc.as.radialSG.d1(1));
    text(0.4,0.52,['At $\,' rx(p) ' = ' rx(p) '_i: \, \, \, \, d/2 \leq ' zy(p) ' \leq d/2+d_1$'],h.latex_r{:})
    %-------------------------------------------------------------------%
    if get(h.es.cb.energy,'value') && strcmp(h.b2.bc.z(1,3),'c')
        axes(h.bc.as.radialSG.d1(2)); cla(h.bc.as.radialSG.d1(2));
        text(0.94*h.tp_r,h.tp_y,'$^\circ$C',h.latex_r{:})
        if get(h.bc.pm.radialSG.d1.temperature,'value') == 2
            text(0.4*h.tp_l,h.tp_y,['$T_{' rx(p) '_i}(' zy(p) ') =$'],h.latex_l{:})
        else
            text(0.4*h.tp_l,h.tp_y,['$T_{' rx(p) '_i} =$'],h.latex_l{:})
        end
        %-------------------------------------------------------------------%
        displayTemperature = get(h.bc.et.radialSG.d1.temperature,'String');
        idx = strfind(displayTemperature,zy(3-p));
        if ~isempty(idx)
            displayTemperature(idx) = zy(p);
            set(h.bc.et.radialSG.d1.temperature,'String',displayTemperature);
        end
    end
    %-------------------------------------------------------------------%
    axes(h.bc.as.radialSG.d2(1)); cla(h.bc.as.radialSG.d2(1));
    text(0.4,0.52,['At $\,' rx(p) ' = ' rx(p) '_i: \, \, \, \, -d/2-d_2 \leq ' zy(p) ' \leq -d/2$'],h.latex_r{:})
    %-------------------------------------------------------------------%
    if get(h.es.cb.energy,'value') && strcmp(h.b2.bc.z(3,3),'c')
        axes(h.bc.as.radialSG.d2(2)); cla(h.bc.as.radialSG.d2(2));
        text(0.94*h.tp_r,h.tp_y,'$^\circ$C',h.latex_r{:})
        if get(h.bc.pm.radialSG.d2.temperature,'value') == 2
            text(0.4*h.tp_l,h.tp_y,['$T_{' rx(p) '_i}(' zy(p) ') =$'],h.latex_l{:})
        else
            text(0.4*h.tp_l,h.tp_y,['$T_{' rx(p) '_i} =$'],h.latex_l{:})
        end
        %-------------------------------------------------------------------%
        displayTemperature = get(h.bc.et.radialSG.d2.temperature,'String');
        idx = strfind(displayTemperature,zy(3-p));
        if ~isempty(idx)
            displayTemperature(idx) = zy(p);
            set(h.bc.et.radialSG.d2.temperature,'String',displayTemperature);
        end
    end
    %-------------------------------------------------------------------%
    axes(h.bc.as.radialSG.ro(1)); cla(h.bc.as.radialSG.ro(1));
    text(0.4,0.52,['At $\,' rx(p) ' = ' rx(p) '_o: \, \, \, \, -d/2-d_2 \leq ' zy(p) ' \leq d/2+d_1$'],h.latex_r{:})
    %-------------------------------------------------------------------%
    if get(h.es.cb.energy,'value') && strcmp(h.b2.bc.z(4,3),'c')
        axes(h.bc.as.radialSG.ro(2)); cla(h.bc.as.radialSG.ro(2));
        text(0.94*h.tp_r,h.tp_y,'$^\circ$C',h.latex_r{:})
        if get(h.bc.pm.radialSG.ro.temperature,'value') == 2
            text(0.4*h.tp_l,h.tp_y,['$T_{' rx(p) '_o}(' zy(p) ') =$'],h.latex_l{:})
        else
            text(0.4*h.tp_l,h.tp_y,['$T_{' rx(p) '_o} =$'],h.latex_l{:})
        end
        %-------------------------------------------------------------------%
        displayTemperature = get(h.bc.et.radialSG.ro.temperature,'String');
        idx = strfind(displayTemperature,zy(3-p));
        if ~isempty(idx)
            displayTemperature(idx) = zy(p);
            set(h.bc.et.radialSG.ro.temperature,'String',displayTemperature);
        end
    end
% ---------------------------------------------------------------------- %
% ---------------------------axial-Surrounding-Gas---------------------- %
% ---------------------------------------------------------------------- %
    axes(h.bc.as.axialSG.d1(1)); cla(h.bc.as.axialSG.d1(1));
    text(0.4,0.52,['At $\,' zy(p) ' = d/2+d_1: \, \, \, \, ' rx(p) '_i \leq ' rx(p) ' \leq ' rx(p) '_o$'],h.latex_r{:})
    %-------------------------------------------------------------------%
    if get(h.es.cb.energy,'value') && strcmp(h.b2.bc.r(1,3),'c')
        axes(h.bc.as.axialSG.d1(2)); cla(h.bc.as.axialSG.d1(2));
        text(0.94*h.tp_r,h.tp_y,'$^\circ$C',h.latex_r{:})
        if get(h.bc.pm.axialSG.d1.temperature,'value') == 2
            text(0.4*h.tp_l,h.tp_y,['$T_{1}(' rx(p) ') =$'],h.latex_l{:})
        else
            text(0.4*h.tp_l,h.tp_y,'$T_{1} =$',h.latex_l{:})
        end
        %-------------------------------------------------------------------%
        displayTemperature = get(h.bc.et.axialSG.d1.temperature,'String');
        idx = strfind(displayTemperature,rx(3-p));
        if ~isempty(idx)
            displayTemperature(idx) = rx(p);
            set(h.bc.et.axialSG.d1.temperature,'String',displayTemperature);
        end
    end
    %-------------------------------------------------------------------%
    if strcmp(h.b2.bc.r(1,1),'i')
        axes(h.bc.as.axialSG.d1(3)); cla(h.bc.as.axialSG.d1(3));
        text(0.94*h.tp_r,h.tp_y,'m/s',h.latex_r{:})
        if get(h.bc.pm.axialSG.d1.velocity,'value') == 2
            text(0.67*h.tp_l,h.tp_y,['$w_\mathrm{g}(' rx(p) ') =$'],h.latex_l{:})
        else
            text(0.67*h.tp_l,h.tp_y,'$\bar{w}_\mathrm{g} =$',h.latex_l{:})
        end
        %-------------------------------------------------------------------%
        displayVelocity = get(h.bc.et.axialSG.d1.velocity,'String');
        idx = strfind(displayVelocity,rx(3-p));
        if ~isempty(idx)
            displayVelocity(idx) = rx(p);
            set(h.bc.et.axialSG.d1.velocity,'String',displayVelocity);
        end
    end
    %-------------------------------------------------------------------%
    axes(h.bc.as.axialSG.d2(1)); cla(h.bc.as.axialSG.d2(1));
    text(0.4,0.52,['At $\,' zy(p) ' = -d/2-d_2: \, \, \, \, ' rx(p) '_i \leq ' rx(p) ' \leq ' rx(p) '_o$'],h.latex_r{:})
    %-------------------------------------------------------------------%
    if get(h.es.cb.energy,'value') && strcmp(h.b2.bc.r(2,3),'c')
        axes(h.bc.as.axialSG.d2(2)); cla(h.bc.as.axialSG.d2(2));
        text(0.94*h.tp_r,h.tp_y,'$^\circ$C',h.latex_r{:})
        if get(h.bc.pm.axialSG.d2.temperature,'value') == 2
            text(0.4*h.tp_l,h.tp_y,['$T_{2}(' rx(p) ') =$'],h.latex_l{:})
        else
            text(0.4*h.tp_l,h.tp_y,'$T_{2} =$',h.latex_l{:})
        end
        %-------------------------------------------------------------------%
        displayTemperature = get(h.bc.et.axialSG.d2.temperature,'String');
        idx = strfind(displayTemperature,rx(3-p));
        if ~isempty(idx)
            displayTemperature(idx) = rx(p);
            set(h.bc.et.axialSG.d2.temperature,'String',displayTemperature);
        end
    end
    %-------------------------------------------------------------------%
    if strcmp(h.b2.bc.r(2,1),'i')
        axes(h.bc.as.axialSG.d2(3)); cla(h.bc.as.axialSG.d2(3));
        text(0.94*h.tp_r,h.tp_y,'m/s',h.latex_r{:})
        if get(h.bc.pm.axialSG.d2.velocity,'value') == 2
            text(0.67*h.tp_l,h.tp_y,['$w_\mathrm{g}(' rx(p) ') =$'],h.latex_l{:})
        else
            text(0.67*h.tp_l,h.tp_y,'$\bar{w}_\mathrm{g} =$',h.latex_l{:})
        end
        %-------------------------------------------------------------------%
        displayVelocity = get(h.bc.et.axialSG.d1.velocity,'String');
        idx = strfind(displayVelocity,rx(3-p));
        if ~isempty(idx)
            displayVelocity(idx) = rx(p);
            set(h.bc.et.axialSG.d1.velocity,'String',displayVelocity);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% stability panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(h.sy.as.m_start); cla(h.sy.as.m_start);
    text(h.tp_l,h.tp_y,['$' mk(p) '_\mathrm{start} = $'],h.latex_l{:})
    % ----------------------------------------------------------------- %
    axes(h.sy.as.m_delta); cla(h.sy.as.m_delta);
    text(h.tp_l,h.tp_y,['$\Delta ' mk(p) ' = $'],h.latex_l{:})
    % ----------------------------------------------------------------- %
    axes(h.sy.as.m_end); cla(h.sy.as.m_end);
    text(h.tp_l,h.tp_y,['$' mk(p) '_\mathrm{end} = $'],h.latex_l{:})