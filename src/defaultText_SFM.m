%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% general panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(h.gl.as.g);
    text(h.tp_l,h.tp_y,'$g=$',h.latex_l{:})
    text(h.tp_r,h.tp_y,'$\cdot \, 9.81 \, \mathrm{m/s}^2$',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.gl.as.Vr);
    text(h.tp_l,h.tp_y,'$\mathcal{V}=V_\mathrm{l} / V_\mathrm{c} =$',h.latex_l{:})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% geometry panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(h.gy.as.llb);
    text(h.tp_l,h.tp_y,'$d=$',h.latex_l{:})
    text(h.tp_r,h.tp_y,'mm',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.gy.as.rc);
    text(h.tp_l,h.tp_y,'$r_c=$',h.latex_l{:})
    text(h.tp_r,h.tp_y,'mm',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.gy.as.ri);
    text(h.tp_l,h.tp_y,'$r_i=$',h.latex_l{:})
    text(h.tp_r,h.tp_y,'mm',h.latex_r{:})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% mesh panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------- %
% ------------------------radial-Liquid-Bridge------------------------- %
% --------------------------------------------------------------------- %
    for i = 1:3
        axes(h.mh.as.radialLB.spacing(i))
        if i == 2
            text(h.tp_l,h.tp_y,'$\Delta_\mathrm{max}=$',h.latex_l{:})
        else
            text(h.tp_l,h.tp_y,'$\Delta_\mathrm{min}=$',h.latex_l{:})
        end
        text(h.tp_r,h.tp_y,'$\%\,$ of $\,d$',h.latex_r{:})
        %---------------------------------------------------------------%
        axes(h.mh.as.radialLB.stretchingC(i));
        text(h.tp_l,h.tp_y,'$f=$',h.latex_l{:})
    end
    axes(h.mh.as.radialLB.spacing(1))
    text(-5.2,h.tp_y,'At $\,r=0$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.radialLB.spacing(2))
    text(-5.2,h.tp_y,'At $\,r=r_i/2$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.radialLB.spacing(3))
    text(-5.2,h.tp_y,'At $\,r=r_i$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.radialLB.stretchingC(1));
    text(-5.2,h.tp_y,'At $\,r \cong 0$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.radialLB.stretchingC(2));
    text(-5.2,h.tp_y,'At $\,r \cong r_i/2$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.radialLB.stretchingC(3));
    text(-5.2,h.tp_y,'At $\,r \cong r_i$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.radialLB.stretchingF(1));
    text(0.3,0.52,'At $\,r \cong 0$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.radialLB.stretchingF(2));
    text(0.3,0.52,'At $\,r \cong r_i/2$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.radialLB.stretchingF(3));
    text(0.3,0.52,'At $\,r \cong r_i$:',h.latex_r{:})
% --------------------------------------------------------------------- %
% -------------------------axial-Liquid-Bridge------------------------- %
% --------------------------------------------------------------------- %
    for i = 1:3
        axes(h.mh.as.axialLB.spacing(i))
        if i == 2
            text(h.tp_l,h.tp_y,'$\Delta_\mathrm{max}=$',h.latex_l{:})
        else
            text(h.tp_l,h.tp_y,'$\Delta_\mathrm{min}=$',h.latex_l{:})
        end
        text(h.tp_r,h.tp_y,'$\%\,$ of $\,d$',h.latex_r{:})
        %---------------------------------------------------------------%
        axes(h.mh.as.axialLB.stretchingC(i));
        text(h.tp_l,h.tp_y,'$f=$',h.latex_l{:})
    end
    axes(h.mh.as.axialLB.spacing(1))
    text(-5.2,h.tp_y,'At $\,z=d/2$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.axialLB.spacing(2))
    text(-5.2,h.tp_y,'At $\,z=0$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.axialLB.spacing(3))
    text(-5.2,h.tp_y,'At $\,z = -d/2$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.axialLB.stretchingC(1));
    text(-5.2,h.tp_y,'At $\,z \cong d/2$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.axialLB.stretchingC(2));
    text(-5.2,h.tp_y,'At $\,z \cong 0$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.axialLB.stretchingC(3));
    text(-5.2,h.tp_y,'At $\,z \cong -d/2$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.axialLB.stretchingF(1));
    text(0.3,0.52,'At $\,z \cong d/2$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.axialLB.stretchingF(2));
    text(0.3,0.52,'At $\,z \cong 0$:',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.mh.as.axialLB.stretchingF(3));
    text(0.3,0.52,'At $\,z \cong -d/2$:',h.latex_r{:})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% equations panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(h.es.as.continuity);
    text(0,0,'$\partial_t \rho+\nabla \cdot (\rho \vec u)=0$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
    %-------------------------------------------------------------------%
    axes(h.es.as.momentum);
    text(0,0,'$\partial_t(\rho \vec u) + \nabla\cdot(\rho \vec u \vec u) = -\nabla  p + \nabla \cdot (\mu{ \boldmath \tau})+\rho\vec{g}$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center')
    %-------------------------------------------------------------------%
    axes(h.es.as.energy);
    text(0,0,'$\partial_t(\rho c_p T)+\nabla \cdot (\rho c_p T \vec u)=\nabla \cdot (\lambda \nabla T)$','interpreter','Latex','FontSize',16,'HorizontalAlignment','Center');
    %-------------------------------------------------------------------%
    if ispc
        strt = 4.2;
    else
        strt = 7;
    end
    dlt = 1;
    axes(h.es.as.properties);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% boundary conditions panel %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ----------------------------------------------------------------- %
    % ------------------------radial-Liquid-Bridge--------------------- %
    % ----------------------------------------------------------------- %
    axes(h.bc.as.radialLB.rc(1));
    text(0.4,0.52,'At $\,r = 0: \, \, \, \, -d/2 \leq z \leq d/2$',h.latex_r{:})
    % ----------------------------------------------------------------- %
    axes(h.bc.as.radialLB.ri);
    text(0.4,0.52,'At $\,r = h(z): \, \, \, \, -d/2 \leq z \leq d/2$',h.latex_r{:})
    % --------------------------------------------------------------------- %
    % ------------------------axial-Liquid-Bridge-------------------------- %
    % --------------------------------------------------------------------- %
    axes(h.bc.as.axialLB.d1(1));
    text(0.4,0.52,'At $\,z = d/2: \, \, \, \, 0 \leq r \leq r_i$',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.bc.as.axialLB.d1(2));
    text(0.4*h.tp_l,h.tp_y,'$T_{d1}=$',h.latex_l{:})
    text(0.94*h.tp_r,h.tp_y,'$^\circ$C',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.bc.as.axialLB.d2(1));
    text(0.4,0.52,'At $\,z = -d/2: \, \, \, \, 0 \leq r \leq r_i$',h.latex_r{:})
    %-------------------------------------------------------------------%
    axes(h.bc.as.axialLB.d2(2));
    text(0.4*h.tp_l,h.tp_y,'$T_{d2}=$',h.latex_l{:})
    text(0.94*h.tp_r,h.tp_y,'$^\circ$C',h.latex_r{:})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% simulation panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(h.sn.as.res);
    text(h.tp_l,h.tp_y,'$f\left(\mbox{\boldmath{q}}_0^{(k)}\right)= $',h.latex_l{:})
    %-------------------------------------------------------------------%
    axes(h.sn.as.newton);
    text(h.tp_l,h.tp_y,'$\delta \mbox{\boldmath{q}}= $',h.latex_l{:})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% optical ray panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(h.or.as.ray(1));
    text(h.tp_l,h.tp_y,'$r_0=$',h.latex_l{:})
    text(h.tp_r,h.tp_y,'$\cdot\, r_i$',h.latex_r{:})
    % ----------------------------------------------------------------- %
    axes(h.or.as.ray(2));
    text(-0.51,h.tp_y,'$z_0=d/2$',h.latex_r{:})
    % ----------------------------------------------------------------- %
    axes(h.or.as.particle(1));
    text(-0.53,h.tp_y,'$r_\mathrm{p}=$',h.latex_r{:})
    % ----------------------------------------------------------------- %
    axes(h.or.as.particle(2));
    text(h.tp_l,h.tp_y,'$z_\mathrm{p}=$',h.latex_l{:})
    text(h.tp_r,h.tp_y,'$\cdot\, d$',h.latex_r{:})
    % ----------------------------------------------------------------- %
    axes(h.or.as.N(1));
    text(h.tp_l,h.tp_y,'$\mathcal{N}(T)=$',h.latex_l{:})
    text(h.tp_r-0.05,h.tp_y,'$-$ ',h.latex_r{:})
    % ----------------------------------------------------------------- %
    axes(h.or.as.N(2)); cla(h.or.as.N(2));
    text(h.tp_r-0.05,h.tp_y,'$\cdot \Big(T-$ ',h.latex_r{:})
    % ----------------------------------------------------------------- %
    axes(h.or.as.N(3)); cla(h.or.as.N(3));
    text(h.tp_r-0.05,h.tp_y,'$^\circ \mathrm{C} \Big)$',h.latex_r{:})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% stability panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(h.sy.as.n);
    text(h.tp_l,h.tp_y,'$n = $',h.latex_l{:})
    % ----------------------------------------------------------------- %
    axes(h.sy.as.n_cay);
    text(h.tp_l,h.tp_y,'$n_\mathrm{cayley} = $',h.latex_l{:})
    % ----------------------------------------------------------------- %
    axes(h.sy.as.m_start);
    text(h.tp_l,h.tp_y,'$m_\mathrm{start} = $',h.latex_l{:})
    % ----------------------------------------------------------------- %
    axes(h.sy.as.m_delta);
    text(h.tp_l,h.tp_y,'$\Delta m = $',h.latex_l{:})
    % ----------------------------------------------------------------- %
    axes(h.sy.as.m_end);
    text(h.tp_l,h.tp_y,'$m_\mathrm{end} = $',h.latex_l{:})
    % ----------------------------------------------------------------- %
    axes(h.sy.as.delta_T);
    text(0,0,'$\Delta T$',h.latex_r{:})
    % ----------------------------------------------------------------- %
    axes(h.sy.as.T_d1);
    text(0,0,'$T_{d1}$',h.latex_r{:})
    % ----------------------------------------------------------------- %
    axes(h.sy.as.T_d2);
    text(0,0,'$T_{d2}$',h.latex_r{:})