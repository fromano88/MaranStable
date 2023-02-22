res.u_vec = [res.u_vec res.u];
res.w_vec = [res.w_vec res.w];
res.conti_vec = [res.conti_vec res.p];
res.newton_vec = [res.newton_vec norm(dx)];
if flowopt.energy
    res.T_vec = [res.T_vec res.T];
end

plot(flowopt.tolerance.residuals*ones(length(res.u_vec)+1),'--','Color','k'); hold on;

if flowopt.tolerance.residuals ~= flowopt.tolerance.newton
    plot(flowopt.tolerance.newton*ones(length(res.u_vec)+1),'-.','Color','k'); hold on;
end

if length(res.u_vec) == 1
    h(1) = plot(res.u_vec,'Color',mycolor.green,'Marker','.','MarkerSize',8); hold on;
    h(2) = plot(res.w_vec,'Color',mycolor.ocher,'Marker','.','MarkerSize',8); hold on;
    h(3) = plot(res.conti_vec,'Color',mycolor.violet,'Marker','.','MarkerSize',8); hold on;
    h(4) = plot(res.newton_vec,'Color',mycolor.blue,'Marker','.','MarkerSize',8); hold on;
    if flowopt.energy
        h(5) = plot(res.T_vec,'Color',mycolor.red,'Marker','.','MarkerSize',8); hold off;
        legend(h(1:5),{'u-velocity','w-velocity','Continuity','$\delta \mbox{\boldmath{q}}$','Temperature'},'Interpreter','latex');
    else
        legend(h(1:4),{'u-velocity','w-velocity','Continuity','$\delta \mbox{\boldmath{q}}$'},'Interpreter','latex');
    end
else
    h(1) = plot(res.u_vec,'Color',mycolor.green,'Marker','.','MarkerSize',8,'LineWidth',1); hold on;
    h(2) = plot(res.w_vec,'Color',mycolor.ocher,'Marker','.','MarkerSize',8,'LineWidth',1); hold on;
    h(3) = plot(res.conti_vec,'Color',mycolor.violet,'Marker','.','MarkerSize',8,'LineWidth',1); hold on;
    if flowopt.energy
        h(4) = plot(res.newton_vec,'Color',mycolor.blue,'Marker','.','MarkerSize',8,'LineWidth',1); hold on;
        h(5) = plot(res.T_vec,'Color',mycolor.red,'Marker','.','MarkerSize',8,'LineWidth',1);
        if length(res.u_vec) < 4
            legend(h(1:5),{'u-velocity','w-velocity','Continuity','$\delta \mbox{\boldmath{q}}$','Temperature'},'Interpreter','latex');
        else
            legend(h(1:5),{'u-velocity','w-velocity','Continuity','$\delta \mbox{\boldmath{q}}$','Temperature'},'Location','southwest','Interpreter','latex');
        end
    else
        h(4) = plot(res.newton_vec,'Color',mycolor.blue,'Marker','.','MarkerSize',8,'LineWidth',1);
        if length(res.u_vec) < 4
            legend(h(1:4),{'u-velocity','w-velocity','Continuity','$\delta \mbox{\boldmath{q}}$'},'Interpreter','latex');
        else
            legend(h(1:4),{'u-velocity','w-velocity','Continuity','$\delta \mbox{\boldmath{q}}$'},'Location','southwest','Interpreter','latex');
        end
    end
        
end
set(gca, 'YScale', 'log', 'xtick', 1:length(res.u_vec)+1)
if length(res.u_vec) > 10
    xlim([length(res.u_vec)-10 length(res.u_vec)+1]);
else
    xlim([1 length(res.u_vec)+1]);
end
xlabel('Iterations','Interpreter','latex','FontSize',12); ylabel('Residuals','Interpreter','latex','FontSize',12);
drawnow