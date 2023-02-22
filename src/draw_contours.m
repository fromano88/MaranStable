hold on

for i = 1:length(blocks)
    b = eval(blocks{i});
    plot(b.R.v(:,[1 end]) ,z_shift-b.Z.v(:,[1 end]),'k');
    plot(b.R.v([1 end],:)',z_shift-b.Z.v([1 end],:)','k')
end

if strcmp(lp.axis,'r')
    if length(blocks) == 2
        if sum(lp.e_z)<1 || sum(lp.e_z)>2
            plot([r_i r_o],z_shift-lp.pos*[1 1],'Color',[0.9290 0.6940 0.1250])
        else
            plot([r_c r_o],z_shift-lp.pos*[1 1],'Color',[0.9290 0.6940 0.1250])
        end
    else
        [~, idx] = min(abs(b1.Z.v(1,:)-lp.pos));
        plot(b1.R.v(:,idx),z_shift-b1.Z.v(:,idx),'Color',[0.9290 0.6940 0.1250])
    end
else
    if lp.e_r(2) > 0
        b = b2;
    else
        b = b1;
    end
    if strcmp(lp.coord,'straight')
        plot(lp.pos*[1 1],z_shift-b.geom.z([1 end]),'Color',[0.9290 0.6940 0.1250])
    else
        [~, idx] = min(abs(b.R.v(:,1)-lp.pos));
        plot(b.R.v(idx,:),z_shift-b.Z.v(idx,:),'Color',[0.9290 0.6940 0.1250])
    end
end

% coordinate system
w = 1.5; d = min([l_lb r_i])/30;
plot([0 r_i/4],[0 0],'lineWidth',w,'color','k')
plot([0 0],[0 r_i/4],'lineWidth',w,'color','k')
rectangle('Position',[-d/2 -d/2 d d],'Curvature',[1 1],'FaceColor',[0 0 0],'lineWidth',w)
fill(r_i/4+[0 3*d 0],[-d 0 d],'k','lineWidth',w)
fill([-d 0 d],r_i/4+[0 3*d 0],'k','lineWidth',w)

axx = gca;
axx.Toolbar.Visible = 'off';
disableDefaultInteractivity(axx)

clear axx

axis equal
axis off