% uploading sketch
    axes(h.as.sketch); cla(h.as.sketch);
    drawSketch_TFM
    
% labeling sketch
    text(3.6,7.7,'$d_{1}$',h.Latex{:})
    text(3.6,2.7,'$d_{2}$',h.Latex{:})
    if strcmp(h.b2.bc.r(1,1),'i')
        text(5.5,9.8,'$\bar{w}_\mathrm{g}$',h.Latex{:})
    elseif strcmp(h.b2.bc.r(1,1),'o')
        text(5.5,-0.3,'$\bar{w}_\mathrm{g}$',h.Latex{:})
    end
    
    if h.flowopt.g > 0
        text(6.2,6.2,'$g$',h.Latex{:})
    elseif h.flowopt.g < 0
        text(6.2,6.2,'$g$',h.Latex{:})
    end
    
    if h.flowopt.ax == 1
        text(-0.4,5.8,'$z$',h.LATEX{:})
        text(1,4.9,'$r$',h.LATEX{:})
        text(1.2,2.15,'$r_i$',h.Latex{:})
        text(1.2,1.2,'$r_o$',h.Latex{:})
    else
        text(-0.4,5.8,'$y$',h.LATEX{:})
        text(1,4.9,'$x$',h.LATEX{:})
        text(1.2,2.15,'$x_i$',h.Latex{:})
        text(1.2,1.2,'$x_o$',h.Latex{:})
    end

    if h.r_c == 0
        text(0.1,3.7,'$d$',h.Latex{:})
    else
        text(0.7,3.7,'$d$',h.Latex{:})
        if h.flowopt.ax == 1
            text(0.1,5.3,'$r_c$',h.Latex{:})
        else
            text(0.1,5.3,'$x_c$',h.Latex{:})
        end
    end
    

color = [255 215 0]./255;
t_pause = 0.5;

if exist('updated','var') == 1
    switch updated
        case 'g'
            if h.flowopt.g ~= 0
                if g ~= 0
                    plt = plot(r_i+0.8*delta_r*[1 1],l_d2+l_lb+0.2*l_d1-[0 l.g],'lineWidth',1.2*w,'color',color);
                    if g > 0
                        plt2 = fill(r_i+0.8*delta_r+d*[-1 0 1],l_d2+l_lb+0.2*l_d1-l.g-[0 3*d 0],color,'lineWidth',1.2*w,'EdgeColor',color);
                    else
                        plt2 = fill(r_i+0.8*delta_r+d*[-1 0 1],l_d2+l_lb+0.2*l_d1+[0 3*d 0],color,'lineWidth',1.2*w,'EdgeColor',color);
                    end
                end
            end

        case 'l_d1'
            plt = plot(r_i*[1 1],l_ges-[l_d1 0],'color',color,'lineWidth',2*w);

        case 'l_lb'
            if h.r_c == 0
                plt = plot([0 0],l_d2+[0 l_lb],'color',color,'lineWidth',2*w);
            else
                plt = plot(r_C+[0 0],l_d2+[0 l_lb],'color',color,'lineWidth',2*w);
            end

        case 'l_d2'
            plt = plot(r_i*[1 1],[0 l_d2],'color',color,'lineWidth',2*w);

        case 'r_c'
            plt = plot([0 r_C-3*d],l_d2+vPos.rc*[1 1],'color',color,'lineWidth',w);
            plt2 = fill(r_C-3*d*[1 0 1],l_d2+vPos.rc+d*[-1 0 1],color,'LineStyle','none');

        case 'r_i'
            plt = plot([0 r_i-3*d],vPos.ri*[1 1],'lineWidth',w,'color',color);
            plt2 = fill(r_i-3*d*[1 0 1],vPos.ri+d*[-1 0 1],color,'LineStyle','none');

        case 'r_o'
            plt = plot([0 r_i+delta_r-3*d],vPos.ro*[1 1],'lineWidth',w,'color',color);
            plt2 = fill(r_i+delta_r-3*d*[1 0 1],vPos.ro+d*[-1 0 1],color,'LineStyle','none');
            
        case 'w_in'            
            if w_in > 0
                plt = plot(r_i+delta_r/2*[1 1],l_ges+l.in/2*[-1 1],'lineWidth',1.5*w,'color',color);
                plt2 = fill(r_i+0.5*delta_r+d*[-1 0 1],l_ges-l.in/2-[0 3*d 0],color,'lineWidth',1.5*w,'EdgeColor',color);
            else
                plt = plot(r_i+delta_r/2*[1 1],l.in/2*[-1 1],'lineWidth',1.5*w,'color',color);
                plt2 = fill(r_i+0.5*delta_r+d*[-1 0 1],l.in/2+[0 3*d 0],color,'lineWidth',1.5*w,'EdgeColor',color);
            end
            
        case 'jump'
            for i = 1:numberOfJumps
                pt(i) = rectangle('Position',[eval(pos_jump(i,1))-d eval(pos_jump(i,2))-d 2*d 2*d],'Curvature',[1 1],'FaceColor',[1 1 0],'lineWidth',w);
            end
    end
    pause(t_pause)
    if exist('plt','var') == 1
        set(plt,'visible','off')
    end
    if exist('plt2','var') == 1
        set(plt2,'visible','off')
    end
end