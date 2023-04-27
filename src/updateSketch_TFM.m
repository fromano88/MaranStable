% uploading sketch
    axes(h.as.sketch); cla(h.as.sketch);
    drawSketch_TFM
    
% labeling sketch
    text(1.76,4,'$d$',h.Latex{:},'FontSize',20)
    text(4.1,7.7,'$d_{1}$',h.Latex{:},'FontSize',20)
    text(4.1,1.7,'$d_{2}$',h.Latex{:},'FontSize',20)
    if strcmp(h.b2.bc.r(1,1),'i')
        text(5.5,9.8,'$\bar{w}_\mathrm{g}$',h.Latex{:},'FontSize',20)
    elseif strcmp(h.b2.bc.r(1,1),'o')
        text(5.5,-0.3,'$\bar{w}_\mathrm{g}$',h.Latex{:},'FontSize',20)
    end
    
    if h.flowopt.g > 0
        text(6.2,6.2,'$g$',h.Latex{:},'FontSize',20)
    elseif h.flowopt.g < 0
        text(6.2,6.2,'$g$',h.Latex{:},'FontSize',20)
    end
    
    if h.flowopt.ax == 1
        text(0.12,5.8,'$z$',h.Latex{:})
        text(1,4.9,'$r$',h.Latex{:})
        text(1.5,2.12,'$r_i$',h.Latex{:},'FontSize',20)
        text(1.5,1.17,'$r_o$',h.Latex{:},'FontSize',20)
    else
        text(0.12,5.8,'$y$',h.Latex{:})
        text(1,4.9,'$x$',h.Latex{:})
        text(1.5,2.12,'$x_i$',h.Latex{:},'FontSize',20)
        text(1.5,1.17,'$x_o$',h.Latex{:},'FontSize',20)
    end

    if h.r_c > 0
        if h.flowopt.ax == 1
            text(0.17,3.85,'$r_c$',h.Latex{:},'FontSize',20)
        else
            text(0.17,3.85,'$x_c$',h.Latex{:},'FontSize',20)
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
            plt  = plot(hPos.l_d1*[1 1],l_d2+l_lb+[0 l_d1],'lineWidth',w*0.7,'color',color);
            plt2 = fill(hPos.l_d1+0.75*d*[-1 0 1 0],l_ges-2.5*d*[1 0 1 0.6],color,'LineStyle','none');
            plt3 = fill(hPos.l_d1+0.75*d*[-1 0 1 0],l_d2+l_lb+2.5*d*[1 0 1 0.6],color,'LineStyle','none');

        case 'l_lb'
            plt = plot(hPos.l_lb*[1 1],l_d2+[0 l_lb],'lineWidth',w*0.7,'color',color);
            plt2 = fill(hPos.l_lb+0.75*d*[-1 0 1 0],l_d2+l_lb-2.5*d*[1 0 1 0.6],color,'LineStyle','none');
            plt3 = fill(hPos.l_lb+0.75*d*[-1 0 1 0],l_d2+2.5*d*[1 0 1 0.6],color,'LineStyle','none');

        case 'l_d2'
            plt  = plot(hPos.l_d1*[1 1],[0 l_d1],'lineWidth',w*0.7,'color',color);
            plt2 = fill(hPos.l_d1+0.75*d*[-1 0 1 0],l_d2-2.5*d*[1 0 1 0.6],color,'LineStyle','none');
            plt3 = fill(hPos.l_d1+0.75*d*[-1 0 1 0],2.5*d*[1 0 1 0.6],color,'LineStyle','none');

        case 'r_c'
            plt = plot([0 r_C],l_d2+vPos.rc*[1 1],'color',color,'lineWidth',w);
            plt2 = fill(r_C-2.0*d*[1 0 1 0.6],l_d2+vPos.rc+0.75*d*[-1 0 1 0],color,'LineStyle','none');
            if ax == 0 && strcmp(bc_r0,'w')
                plt3 = fill(2.0*d*[1 0 1 0.6],l_d2+vPos.rc+0.75*d*[-1 0 1 0],color,'LineStyle','none');
            end

        case 'r_i'
            plt = plot([0 r_i],vPos.ri*[1 1],'lineWidth',0.7*w,'color',color);
            plt2 = fill(r_i-2.5*d*[1 0 1 0.6],vPos.ri+0.75*d*[-1 0 1 0],color,'LineStyle','none');
            if ax == 0 && strcmp(bc_r0,'w')
                plt3 = fill(2.5*d*[1 0 1 0.6],vPos.ri-0.75*d*[-1 0 1 0],color,'LineStyle','none');
            end

        case 'r_o'
            plt = plot([0 r_i+delta_r],vPos.ro*[1 1],'lineWidth',0.7*w,'color',color);
            plt2 = fill(r_i+delta_r-2.5*d*[1 0 1 0.6],vPos.ro+0.75*d*[-1 0 1 0],color,'LineStyle','none');
            if ax == 0 && strcmp(bc_r0,'w')
                plt3 = fill(2.5*d*[1 0 1 0.6],vPos.ro-0.75*d*[-1 0 1 0],color,'LineStyle','none');
            end
            
        case 'w_in'            
            if w_in > 0
                plt = plot(r_i+delta_r/2*[1 1],l_ges+l.in/2*[-1 1],'lineWidth',1.5*w,'color',color);
                plt2 = fill(r_i+0.5*delta_r+d*[-1 0 1],l_ges-l.in/2-[0 3*d 0],color,'lineWidth',1.5*w,'EdgeColor',color);
            else
                plt = plot(r_i+delta_r/2*[1 1],l.in/2*[-1 1],'lineWidth',1.5*w,'color',color);
                plt2 = fill(r_i+0.5*delta_r+d*[-1 0 1],l.in/2+[0 3*d 0],color,'lineWidth',1.5*w,'EdgeColor',color);
            end
            
    end
    pause(t_pause)
    if exist('plt','var') == 1
        set(plt,'visible','off')
    end
    if exist('plt2','var') == 1
        set(plt2,'visible','off')
    end
    if exist('plt3','var') == 1
        set(plt3,'visible','off')
    end
end