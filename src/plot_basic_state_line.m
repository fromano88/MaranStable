if isempty(get(gcf,'Name'))
    close gcf
    set(gcf,'color','white');
    % scaling
        unity = 'm'; dim = 1;
    % line plot: (1) leftPlot, (2) rightPlot
        lp.n = 1; lp.axis = 'z'; lp.zoom = [0 1]; lp.coord = 'curved';
        lp.domain   = [1 1]; % [1 0]-->LB, [0 1]-->SG, [1 1]-->all
        lp.quantity = ["w", "p"];
    % plotting position
        % z = e_z(1)*l_d1 + e_z(2)*l_lb + e_z(3)*l_d2
            lp.e_z = [1 0.5 0];
        % r = r_c + e_r(1)*r_i + e_r(2)*r_o
            lp.e_r = [1 0.01];
end

clear tl; hold on

% needed to non-dimensionialize
    mu      = str2func(b1.mu);
    dsigma  = str2func(b1.dsigma);
    lambda  = str2func(b1.lambda);
    T_0     = (T_d1l+T_d2l)/2;
    delta_T = abs(T_d1l-T_d2l);

% coordinate shift
    z_shift = b1.geom.z(1)+0.5*diff(b1.geom.z);

% colors
    lp.colors = [0      0.4470 0.7410;
                 0.8500 0.3250 0.0980];
    
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

% adjust domain for SFM
    if length(blocks) == 1
        lp.domain = [1 0];
        [r_o, l_d1, l_d2] = deal(0);
        [lp.e_r(2), lp.e_z(1), lp.e_z(3)] = deal(0,1,0);
    end

% adjust coordinates if h, Bi or n_nabla_T is plotted
if ismember(lp.quantity(1), {'h', 'n_nabla_T', 'Bi'}) || ismember(lp.quantity(lp.n), {'h', 'n_nabla_T', 'Bi'})
    lp.coord = 'curved';
end

% adjust domain according to epsilon and compute position
    if strcmp(lp.axis,'z')
        lp.pos = r_c + lp.e_r(1)*(r_i-r_c) + lp.e_r(2)*(r_o-r_i);
        z_fs = []; % coordinates of intersection between h_fs and straight line

        if strcmp(lp.coord,'curved')
            if lp.e_r(2) > 0
                lp.domain = [0 1];
            else
                lp.domain = [1 0];
            end
        else
            h_fs = b1.R.v(end,:);
            if prod(lp.pos<h_fs) == 1
                lp.domain = [1 0];
            elseif prod(lp.pos>h_fs) == 1
                lp.domain = [0 1];
            end

            % find index shift between LB and SG
            if length(blocks) == 2
                i_shift(1) = find(b2.Z.v(1,:)==b1.geom.z(1),1)-1;
                i_shift(2) = size(b2.Z.v,2)-find(b2.Z.v(1,:)==b1.geom.z(end),1,"last");
            end

            % find intersection of h_fs and lp.pos
            if max(abs(h_fs-r_i)/r_i) < 1e-12
                h_fs = r_i*ones(size(h_fs));
            end
            idx = lp.pos<=h_fs;
            if sum(idx)==2 && idx(1)*idx(end)==1 % only the attached points belong to the free surface (V_r<1, pos=r_i)
                j_fs = [1 length(idx)];
                z_fs = b1.Z.v(1,[1 end]);
            else
                h_fs_interp = @(z) interp1(b1.Z.v(end,:),h_fs,z,'spline');
                j_fs = find(abs(diff(idx))==1);
                for j = j_fs % find every intersection point
                    z_fs = [z_fs fsolve(@(x)h_fs_interp(x)-lp.pos,b1.Z.v(1,j),optimset('Display','off'))]; %#ok<AGROW> 
                end
            end
            if ~isempty(j_fs) && length(blocks) == 2
                j_fs = [j_fs; j_fs+i_shift(1)];
            end
        end

        tl = tiledlayout(1,1);
        if lp.n == 1
            ax = axes(tl);
        else
            ax = [axes(tl) axes(tl)];
        end
        
    else    
        if sum(lp.e_z) < 1 || sum(lp.e_z) > 2
            lp.domain = [0 1];
        end
    end

% transform lp.domain into 1, 2 or [1 2] if not already transformed by GUI
    if length(lp.domain)==2 && sum(lp.domain) < 3
        lp.domain = nonzeros(lp.domain.*[1 2])';
    end
    
    
% ----------- plotting loop -------------
lp.merged_value{1} = []; lp.merged_value{2} = []; % needed only for straight plot over z over 2 domains
p = flowopt.ax+1; % pointer
for i = lp.domain
    b = eval(blocks{i});

    for n = 1:lp.n
        switch lp.quantity(n)
            case 'T'
                lp.value = b.T.v;
                if dim == 1
                    lp.label = '$T\; [^{\circ}$C]';
                else
                    lp.value = (lp.value-T_0)/delta_T;
                    lp.label = '$\vartheta$';
                end
    
            case 'p'
                if isfield(b1,'dp_lb') && i==1
                    lp.value = b.p.v+b1.dp_lb;
                else
                    lp.value = b.p.v;
                end
    
                if dim == 1
                    lp.label = '$p$ [Pa]';
                else
                    lp.value = lp.value*l_lb/(-dsigma(T_0)*delta_T);
                    lp.label = '$p$';
                end
                
            case 'w'
                lp.value = -1./b.ETA_z.v.*b.w.v;
                if dim == 1
                    lp.value = lp.value*f;
                    lp.label = strcat('$w$ [',unity,'/s]');
                else
                    lp.value = lp.value*mu(T_0)/(-dsigma(T_0)*delta_T);
                    lp.label = '$w$';
                end
                
            case 'u'
                lp.value = 1./b.XI_r.v.*b.u.v-b.JA.v.*b.XI_z.v.*b.w.v;
                if dim == 1
                    lp.value = lp.value*f;
                    lp.label = strcat('$u$ [',unity,'/s]');
                else
                    lp.value = lp.value*mu(T_0)/(-dsigma(T_0)*delta_T);
                    lp.label = '$u$';
                end
                
            case 'norm_u'
                lp.value = sqrt((1./b.ETA_z.v.*b.w.v).^2+(1./b.XI_r.v.*b.u.v-b.JA.v.*b.XI_z.v.*b.w.v).^2);
                if dim == 1
                    lp.value = lp.value*f;
                    lp.label = strcat('$| \vec{u} |$ [',unity,'/s]');
                else
                    lp.value = lp.value*mu(T_0)/(-dsigma(T_0)*delta_T);
                    lp.label = '$| \vec{u} |$';
                end
                
            case 'h'
                lp.value = (b1.R.v(end,:)-r_i)*f;
                if dim == 1
                    lp.label = ['h [' unity ']'];
                else
                    lp.label = '$h$';
                end
                
            case {'n_nabla_T','Bi'}
                N = sqrt(1+b1.R_z.u(end,:).^2);
                n_nabla_T = N(2:end-1).^2.*b1.XI_r.T(end,2:end-1).*(b1.T.T(end,2:end-1)-b1.T.T(end-1,2:end-1))./b1.DXI.c(end,:).*2-b1.R_z.T(end,2:end-1).*b1.ETA_z.T(end,2:end-1).*(b1.T.w(end,2:end)-b1.T.w(end,1:end-1))./b1.DETA.c(end,:);
                if strcmp(lp.quantity(n),'n_nabla_T')
                    if dim == 1
                        lp.value = n_nabla_T*lambda(T_0);
                        lp.label = '$\lambda\vec{n}\cdot \nabla T\; [\mathrm{W/m}^2$]';
                    else
                        lp.value = n_nabla_T*l_lb*f/delta_T;
                        lp.label = '$\vec{n}\cdot \nabla \vartheta$';
                    end
                else
                    lp.value = -n_nabla_T.*l_lb./(b1.T.T(end,2:end-1)-T_d2l);
                    lp.label = '$\mathrm{Bi}$';
                end
                lp.value = interp1(b1.ETA.u(end,2:end-1),lp.value,b1.ETA.v(end,2:end-1),'spline');
                lp.value = [NaN lp.value NaN];
        end

        if strcmp(lp.axis,'z')
            lp.xplot = (z_shift-b.Z.v(end,:))*f;

            if strcmp(lp.coord,'straight') % straight plot
                lp.value_interp = scatteredInterpolant(b.R.v(:),b.Z.v(:),lp.value(:),'linear','none');
                lp.value = lp.value_interp(lp.pos*ones(size(b.R.v(end,:))),b.Z.v(end,:));

                if i == 1 && prod(lp.pos<=h_fs)==0
                    h_fs = b1.R.v(end,:);
                    idx = lp.pos<=h_fs;
                elseif i == 2 && prod(lp.pos>h_fs)==0
                    h_fs = b2.R.v(1,:);
                    idx = lp.pos>h_fs;
                else
                    idx = true(size(b.R.v(1,:)));
                end
                lp.true_value = lp.value.*idx;
                lp.true_value(isnan(lp.true_value)) = 0;

                if length(lp.domain) == 1 % only 1 domain shall be plotted
                    if prod(idx) ~= 1
                        lp.value(~lp.true_value) = NaN;
                    else
                        lp.value = lp.true_value;
                    end
                    if ~isempty(j_fs)     % add another value is straight line intersects with free surface
                        lp.c.xplot  = lp.xplot;                                      % copy of x-vector
                        lp.c.value  = lp.value;                                      % copy of value to plot
                        lp.value_fs = lp.value_interp(lp.pos*ones(size(z_fs)),z_fs); % values on the free surface
                        j_fs_end = [j_fs(i,:) length(lp.xplot)];
                        
                        lp.xplot = lp.c.xplot(1:j_fs(1));
                        lp.value = lp.c.value(1:j_fs(1));
                        for j = 1:length(j_fs_end)-1
                            lp.value = [lp.value lp.value_fs(j) lp.c.value(j_fs_end(j)+1:j_fs_end(j+1))];
                            lp.xplot = [lp.xplot (z_shift-z_fs(j))*f lp.c.xplot(j_fs_end(j)+1:j_fs_end(j+1))];
                        end
                    end
                else
                    if lp.e_r(2) > 0 % SG is the main domain
                        lp.xplot = (z_shift-b2.Z.v(end,:))*f;
                        lp.value_fs(n,:) = lp.value_interp(lp.pos*ones(size(z_fs)),z_fs); % values on the free surface

                        if isempty(lp.merged_value{n})
                            lp.merged_value{n} = [zeros(1,i_shift(1)) lp.true_value zeros(1,i_shift(2))];
                        else
                            lp.merged_value{n} = lp.merged_value{n}+lp.true_value;
                        end
                    else % LB is the main domain
                        lp.xplot = (z_shift-b1.Z.v(end,:))*f;
                        if sum(idx)==2 && idx(1)*idx(end)==1 % only the attached points belong to the free surface (V_r<1, pos=r_i)
                            lp.value_fs(n,:) = lp.true_value(idx);
                        else
                            lp.value_fs(n,:) = lp.value_interp(lp.pos*ones(size(z_fs)),z_fs); % values on the free surface
                        end

                        if isempty(lp.merged_value{n})
                            lp.merged_value{n} = lp.true_value;
                        else
                            lp.merged_value{n} = lp.merged_value{n}+lp.true_value(i_shift(1)+(1:size(b1.R.v,2)));
                        end
                    end
                    % add another value if straight line intersects with free surface
                    if i==2 && ~isempty(j_fs) && (j_fs(1,1)~=1 || j_fs(1,2)~=size(b1.R.v,2))
                        lp.value_fs(n,:) = lp.value_interp(lp.pos*ones(size(z_fs)),z_fs);
                        lp.c.value = lp.merged_value{n}; % copy of value to plot
                        lp.c.xplot = lp.xplot;           % copy of x-vector
                        j_fs_end = [j_fs(ceil(sum(lp.e_r)),:) length(lp.xplot)];
                        
                        lp.xplot = lp.c.xplot(1:j_fs_end(1));
                        lp.merged_value{n} = lp.c.value(1:j_fs_end(1));
                        for j = 1:length(j_fs_end)-1
                            lp.merged_value{n} = [lp.merged_value{n} lp.value_fs(n,j) lp.c.value(j_fs_end(j)+1:j_fs_end(j+1))];
                            lp.xplot = [lp.xplot (z_shift-z_fs(j))*f lp.c.xplot(j_fs_end(j)+1:j_fs_end(j+1))];
                        end
                    end
                end

            else % curved plot
                if size(lp.value,1) ~= 1
                    [~, idx] = min(abs(b.R.v(:,1)-lp.pos));
                    lp.value = lp.value(idx,:);
                end
            end

            % plotting
            hold on
            if isempty(lp.merged_value{n})
                plot(ax(n),lp.value,lp.xplot,'Color',lp.colors(n,:));
            elseif i == 2 % plot merged_value only at end of domain loop because it needs to be merged, first
                plot(ax(n),lp.merged_value{n},lp.xplot,'Color',lp.colors(n,:));
            end
            ax(n).XColor = lp.colors(n,:);
            ax(n).XLabel.String = lp.label;
            ax(n).XLabel.Interpreter = "latex";
            if n == 1
                zy = 'yz';
                if dim == 1
                    ax(1).YLabel.String = ['$' zy(p) '$ [' unity ']'];
                else
                    ax(1).YLabel.String = ['$' zy(p) '$'];
                end
                ax(1).YLabel.Interpreter = "latex";
                % horizontal lines
                if i == lp.domain(end) && ~isempty(z_fs)
                    yline((z_shift-z_fs)*f,'k--','LineWidth',0.1)
                end
                ax(1).Box = 'on';
            else
                ax(2).XAxisLocation = 'top';
                ax(2).YAxisLocation = 'right';
                ax(2).YTickLabel = {};
                ax(2).Color = 'none';
                ax(1).Box   = 'off';
                ax(2).Box   = 'off';
            end
            % vertical lines
            if i == lp.domain(end)
                x_lim = get(ax(n),'xlim');
                if x_lim(1)*x_lim(2) < 0
                    xline(0,'Color',lp.colors(n,:),'LineStyle',':','LineWidth',0.25,'Parent',ax(n))
                end
            end
            % limits
            if lp.e_r(2) > 0
                set(ax(n),'ylim',(z_shift-b2.geom.z(end))*f+lp.zoom*(b2.geom.z(end)-b2.geom.z(1))*f)
            else
                set(ax(n),'ylim',(-l_lb*f/2+lp.zoom*l_lb*f))
            end

        else % plot over r
            box on
            if i == 1
                lp.pos = lp.e_z(2);
            else
                lp.pos = (lp.e_z(1)*l_d1 + lp.e_z(2)*l_lb + lp.e_z(3)*l_d2)/(l_d1+l_lb+l_d2);
            end
            lp.pos = lp.e_z(1)*l_d1 + lp.e_z(2)*l_lb + lp.e_z(3)*l_d2;
            [~, idx] = min(abs(b.Z.v(1,:)-lp.pos));

            lp.value = lp.value(:,idx);
            lp.xplot = b.R.v(:,idx)*f;

            % plotting
            if lp.n == 2
                if n == 1
                    yyaxis left; hold on
                else
                    yyaxis right; hold on
                end
            end
            plot(lp.xplot,lp.value,'Color',lp.colors(n,:),'LineStyle','-')
            ylabel(lp.label,'Interpreter','latex','FontSize', 12)
            rx = 'xr';
            if dim == 1
                xlabel(['$' rx(p) '$ [' unity ']'],'Interpreter','latex','FontSize', 12)
            else
                xlabel(['$' rx(p) '$'],'Interpreter','latex','FontSize', 12)
            end

            % creating horizontal and vertical lines
            if i == lp.domain(end)
                y_lim = ylim;
                [~, idx] = min(abs(b1.Z.v(1,:)-lp.pos));
                h_fs     = b1.R.v(end,idx);
                x_lim = [h_fs*f b.geom.r(end)*f];
                if i == 2
                    xline(h_fs*f,'k--','LineWidth',0.1)
                end
                if y_lim(1)*y_lim(2) < 0
                    yline(0,'Color',lp.colors(n,:),'LineStyle',':','LineWidth',0.25)
                end
                xlim(r_c*f+lp.zoom*(x_lim(i)-r_c*f))
            end

        end

    end

end

clear tl

hold off