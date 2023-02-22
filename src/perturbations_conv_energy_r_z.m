if ~exist('pc','var')
    % perturbation contour plot
        pc.draw = 1; pc.domain = 'all'; pc.equal = 1;

    % isolines
        pc.streamlines = 1; pc.isolines = 0;
        pc.numStreamlines = 8; pc.numIsolines = 20;

    % colormaps: jet(500), parula(500), bone(500), myColormap('Cold2Warm'), myColormap('BlueWhiteRed'), myColormap('thermal2'), myColormap('thermal2')
    pc.colormap = myColormap('Cold2Warm');
end

phi=phi_T_max+pi/8;%30;
time=time_conv_energy_max;%zeile/250*2*pi/imag(gamma(1,1));

if isempty(get(gcf,'Name')), clf, end; hold on

if pc.equal == 1
    axis('equal')
end

for n=1:length(blocks)
    if strcmp(pc.domain,'SG')
        n = 2;
    end
    if pc.draw == 1
        [conv_energy,linie] = eval(['contourf(' blocks{n} '.Z.p-b1.geom.z(1)-0.5*(b1.geom.z(2)-b1.geom.z(1)),' blocks{n} '.R.p,' blocks{n} '.conv_energy_hat.*exp(-gamma(1,1)*time+1i*m*phi)+conj(' blocks{n} '.conv_energy_hat.*exp(-gamma(1,1)*time+1i*m*phi)),500)']);
        set(linie,'Edgecolor','none');
    end
    eval(['plot(' blocks{n} '.Z.v(:,[1 end])-b1.geom.z(1)-0.5*(b1.geom.z(2)-b1.geom.z(1)),' blocks{n} '.R.v(:,[1 end]),''k'')'])
    eval(['plot((' blocks{n} '.Z.v([1 end],:)-b1.geom.z(1))''-0.5*(b1.geom.z(2)-b1.geom.z(1)),' blocks{n} '.R.v([1 end],:)'',''k'')'])
    if pc.isolines == 1
        eval(['contour(' blocks{n} '.Z.p-b1.geom.z(1)-0.5*(b1.geom.z(2)-b1.geom.z(1)),' blocks{n} '.R.p,(' blocks{n} '.T.p-(T_d2l+T_d1l)/2)/abs(T_d2l-T_d1l),' num2str(pc.numIsolines) ',''k'');'])
    end
    sf = eval(['full(-cumsum((' blocks{n} '.R.w(2:end-1,1)).*' blocks{n} '.JA.w(2:end-1,1).*rho(' blocks{n} '.T.w(2:end-1,1)).*' blocks{n} '.w.w(2:end-1,1).*' blocks{n} '.DXI.c(:,1))*ones(1,' blocks{n} '.J) + cumsum(' blocks{n} '.R.p.*' blocks{n} '.JA.p.*rho(' blocks{n} '.T.p).*' blocks{n} '.u.p.*' blocks{n} '.DETA.c,2))']);
    if pc.streamlines == 1
        eval(['contour(' blocks{n} '.Z.p-b1.geom.z(1)-0.5*(b1.geom.z(2)-b1.geom.z(1)),' blocks{n} '.R.p,sf,' num2str(pc.numStreamlines) ', ''black'');'])
    end
    
    if strcmp(pc.domain,'LB')
        break
    end
end
caxis(0.01*[-conv_energy_max,conv_energy_max])
caxis manual
colormap(pc.colormap)

p = flowopt.ax + 1;
rx = 'xr'; zy = 'yz';
xlabel(['$' zy(p) '$ in [m]'],'Interpreter','latex','FontSize', 14)
ylabel(['$' rx(p) '$ in [m]'],'Interpreter','latex','FontSize', 14)
pc.cb = colorbar('YTick', 0.8*0.01*[-conv_energy_max,conv_energy_max],'YTickLabel',{'low','high'});
set(get(pc.cb,'ylabel'),'String','perturbation convective energy','Interpreter','latex','FontSize', 14)

if exist('independent_2','var')
    ind2=[independent_2 '=' num2str(eval(independent_2)) ', '];
    
else
    ind2=[];
    
end

if exist('movie','var') && movie == 0 % document
    title([file_name(1:end-4) ': ' independent '=' num2str(eval(independent)) ', ' ind2 dependent '=' num2str(eval(dependent))],'interpreter','none');

end

set(gca,'TickLabelInterpreter','latex')
if isempty(get(gcf,'Name'))
    set(gcf, 'color', 'white');
    
end

hold off