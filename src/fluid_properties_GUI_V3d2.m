function[b]=fluid_properties_GUI_V3d2(substance, fluid, T_fluid, p_fluid, b)
% -  this file originates from fluid_properties.m and should only be
%       used for creating a standalone application of MaranStable
% -  reason: standalone applications don't support sympbolic programming (such as syms)
% -  moreover, for standalone applications only 1 directory should be used (to avoid errors)
% -  fortunately, fluid_properties.m is only file which uses symbolic programming

if ~isempty(T_fluid)
    theta = T_fluid; %#ok<*NASGU> 
end

R = 8.3143;

if isfile([tempdir 'my_fluids.mat']) && ismember(substance,who('-file', [tempdir 'my_fluids.mat'])) % user defined fluids
    load([tempdir 'my_fluids.mat'])
    mu      = eval([substance '.mu']);
    rho     = eval([substance '.rho']);
    lambda  = eval([substance '.lambda']);
    cp      = eval([substance '.cp']);
    dmu     = eval([substance '.dmu']);
    drho    = eval([substance '.drho']);
    dlambda = eval([substance '.dlambda']);
    dcp     = eval([substance '.dcp']);
    if strcmp(fluid,'liquid')
        sigma  = eval([substance '.sigma']);
        dsigma = eval([substance '.dsigma']);
    end

else % already implemented fluids
    if strcmp(fluid, 'liquid')==1
        % liquids
            if strcmp(substance, 'KF-96L-5cs')==1
                %Shin Etsu silicon fluid KF-96-5cs 25°C
                nu_0 = 5e-6; nu_exp = 5.892; sigma_1 = -6.23e-5; % MEIS-2
                Coef_rho    = polyfit([-20 0 25 40 80],[961 939 915 900 862],2);
                Coef_lambda = polyfit([-60 -40 -20 0 25 40 80 120 160 200 240],[0.1560    0.1500    0.1450    0.1390    0.1330    0.1280    0.1160    0.1050    0.0930    0.0820    0.0710],2);
                Coef_cp     = polyfit([-60 -40 -20 0 25 40 80 120 160 200 240],[1505	1535	1565	1595	1630	1650	1713	1773	1831	1890	1950],2);
                
                nu     = [num2str(nu_0) '*exp(' num2str(nu_exp) '*(25-theta)./(273.15+theta))'];
                rho    = [num2str(Coef_rho(1),16) '*theta.^2+' num2str(Coef_rho(2),16) '*theta+' num2str(Coef_rho(3),16)];
                mu     = [nu '.*(' rho ')'];
                lambda = [num2str(Coef_lambda(1),16) '*theta.^2+' num2str(Coef_lambda(2),16) '*theta+' num2str(Coef_lambda(3),16) '-0.133+0.12'];
                cp     = [num2str(Coef_cp(1),16) '*theta.^2+' num2str(Coef_cp(2),16) '*theta+' num2str(Coef_cp(3),16) '-1630+1800'];
                sigma  = ['19.7e-3+' num2str(sigma_1) '*(theta-25)'];
                
                dnu     = [nu '.*-' num2str(nu_exp) '*(273.15+25)./(273.15+theta).^2'];
                drho    = [num2str(2*Coef_rho(1),16) '*theta+' num2str(Coef_rho(2),16)];
                dmu     = [dnu '.*(' rho ') +' nu '.*(' drho ')'];
                dlambda = [num2str(2*Coef_lambda(1),16) '*theta+' num2str(Coef_lambda(2),16)];
                dcp     = [num2str(2*Coef_cp(1),16) '*theta+' num2str(Coef_cp(2),16)];
                ddcp    = num2str(2*Coef_cp(1),16);
                dsigma  = num2str(sigma_1);

            elseif strcmp(substance, 'KF-96L-2cs')==1
                %Shin Etsu silicon fluid KF-96-2cs 25°C
                nu_0 = 2e-6; nu_exp = 5.892; sigma_1 = -7.0e-5; % Romano 2017 Physics of fluids
                Coef_rho    = polyfit([-20 0 25 40 80],[917 898 873 859 823],2);
                Coef_lambda = polyfit([-60 -40 -20 0 25 40 80 120 160 200 240],[0.1560    0.1500    0.1450    0.1390    0.1330    0.1280    0.1160    0.1050    0.0930    0.0820    0.0710],2);
                Coef_cp     = polyfit([-60 -40 -20 0 25 40 80 120 160 200 240],[1505	1535	1565	1595	1630	1650	1713	1773	1831	1890	1950],2);
                
                nu     = [num2str(nu_0) '*exp(' num2str(nu_exp) '*(25-theta)./(273.15+theta))'];
                rho    = [num2str(Coef_rho(1),16) '*theta.^2+' num2str(Coef_rho(2),16) '*theta+' num2str(Coef_rho(3),16)];
                mu     = [nu '.*(' rho ')'];
                lambda = [num2str(Coef_lambda(1),16) '*theta.^2+' num2str(Coef_lambda(2),16) '*theta+' num2str(Coef_lambda(3),16) '-0.133+0.11'];
                cp     = [num2str(Coef_cp(1),16) '*theta.^2+' num2str(Coef_cp(2),16) '*theta+' num2str(Coef_cp(3),16) '-1630+1800'];
                sigma  = ['18.3e-3+' num2str(sigma_1) '*(theta-25)'];
                
                dnu     = [nu '.*-' num2str(nu_exp) '*(273.15+25)./(273.15+theta).^2'];
                drho    = [num2str(2*Coef_rho(1),16) '*theta+' num2str(Coef_rho(2),16)];
                dmu     = [dnu '.*(' rho ') +' nu '.*(' drho ')'];
                dlambda = [num2str(2*Coef_lambda(1),16) '*theta+' num2str(Coef_lambda(2),16)];
                dcp     = [num2str(2*Coef_cp(1),16) '*theta+' num2str(Coef_cp(2),16)];
                ddcp    = num2str(2*Coef_cp(1),16);
                dsigma  = num2str(sigma_1);

            elseif strcmp(substance, 'KF-96L-20cs')==1
                %Shin Etsu silicon fluid KF-96-20cs 25°C
                nu_0 = 20e-6; nu_exp = 5.892; sigma_1 = -6.24e-5; % Yano 2011 J. Phys.: Conf. Ser. 327 012029
                Coef_rho    = polyfit([-40 0 25 50 100],[1008 974 950 930 885],2);
                Coef_lambda = polyfit([-60 -40 -20 0 25 40 80 120 160 200 240],[0.174    0.168    0.163    0.157    0.151    0.146    0.134    0.123    0.112    0.1    0.089],2);
                Coef_cp     = polyfit([-60 -40 -20 0 25 40 80 120 160 200 240],[1435	1465	1493	1525	1560	1593	1642	1701	1760	1820	1880],2);
                
                nu     = [num2str(nu_0) '*exp(' num2str(nu_exp) '*(25-theta)./(273.15+theta))'];
                rho    = [num2str(Coef_rho(1),16) '*theta.^2+' num2str(Coef_rho(2),16) '*theta+' num2str(Coef_rho(3),16)];
                mu     = [nu '.*(' rho ')'];
                lambda = [num2str(Coef_lambda(1),16) '*theta.^2+' num2str(Coef_lambda(2),16) '*theta+' num2str(Coef_lambda(3),16) '-0.151+0.15'];
                cp     = [num2str(Coef_cp(1),16) '*theta.^2+' num2str(Coef_cp(2),16) '*theta+' num2str(Coef_cp(3),16) '-1560+1600'];
                sigma  = ['20.6e-3+' num2str(sigma_1) '*(theta-25)'];
                
                dnu     = [nu '.*-' num2str(nu_exp) '*(273.15+25)./(273.15+theta).^2'];
                drho    = [num2str(2*Coef_rho(1),16) '*theta+' num2str(Coef_rho(2),16)];
                dmu     = [dnu '.*(' rho ') +' nu '.*(' drho ')'];
                dlambda = [num2str(2*Coef_lambda(1),16) '*theta+' num2str(Coef_lambda(2),16)];
                dcp     = [num2str(2*Coef_cp(1),16) '*theta+' num2str(Coef_cp(2),16)];
                ddcp    = num2str(2*Coef_cp(1),16);
                dsigma  = num2str(sigma_1);

            else
                % VDI Heat Atlas 2010 2nd Edition
                load('Coef_liquid.mat')
                load('caloric_critical_data.mat', 'rho_crit', 'T_crit', 'M')
                rho_crit_l = rho_crit.(genvarname(substance));
                T_crit_l   = T_crit.(genvarname(substance));
                M_l        = M.(genvarname(substance))*10^-3;
                Coef_rho     = Coef_rho_liquid.(genvarname(substance)).';
                Coef_mu      = Coef_mu_liquid.(genvarname(substance)).';
                Coef_lambda  = Coef_lambda_liquid.(genvarname(substance)).'.*[1 10^-2 10^-4 10^-7  10^-10];
                Coef_cp      = Coef_cp_liquid.(genvarname(substance)).';
                Coef_sigma_l = Coef_sigma.(genvarname(substance)).';

                rho    = [num2str(rho_crit_l,10) '+' num2str(Coef_rho(1),10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^0.35 +' num2str(Coef_rho(2),10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^(2/3) +' num2str(Coef_rho(3),10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ') +' num2str(Coef_rho(4),10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^(4/3)'];
                mu     = [num2str(Coef_mu(5),10) '*exp(' num2str(Coef_mu(1),10) '.*((' num2str(Coef_mu(3),10) '-273.15-theta)./(273.15+theta-' num2str(Coef_mu(4),10) ')).^(1/3)+' num2str(Coef_mu(2),10) '.*((' num2str(Coef_mu(3),10) '-273.15-theta)./(273.15+theta-' num2str(Coef_mu(4),10) ')).^(4/3))'];
                lambda = [num2str(Coef_lambda(1),10) '+' num2str(Coef_lambda(2),10) '*(273.15+theta) +' num2str(Coef_lambda(3),10) '*(273.15+theta).^2 +' num2str(Coef_lambda(4),10) '*(273.15+theta).^3 +' num2str(Coef_lambda(5),10) '*(273.15+theta).^4'];
                cp     = [num2str(R) '/' num2str(M_l,10) '*(' num2str(Coef_cp(1),10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^(-1) +' num2str(Coef_cp(2),10) '+' num2str(Coef_cp(3),10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ') +' num2str(Coef_cp(4),10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^2 +' num2str(Coef_cp(5),10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^3 +' num2str(Coef_cp(6),10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^4)'];
                sigma  = [num2str(Coef_sigma_l(1),10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^(' num2str(Coef_sigma_l(2),10) '+' num2str(Coef_sigma_l(3),10) '*((273.15+theta)/' num2str(T_crit_l,10) ') +' num2str(Coef_sigma_l(4),10) '*((273.15+theta)/' num2str(T_crit_l,10) ').^2 +' num2str(Coef_sigma_l(5),10) '*((273.15+theta)/' num2str(T_crit_l,10) ').^3)'];
                
                drho    = ['0.35*' num2str(Coef_rho(1),10) '/' num2str(-T_crit_l) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^(-0.65) + 2/3*' num2str(Coef_rho(2),10) '/' num2str(-T_crit_l) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^(-1/3) +' num2str(Coef_rho(3),10) '/' num2str(-T_crit_l) '+ 4/3*' num2str(Coef_rho(4),10) '/' num2str(-T_crit_l) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^(1/3)'];
                dmu     = [mu '.*(1/3*' num2str(Coef_mu(1),10) '.*((' num2str(Coef_mu(3),10) '-273.15-theta)./(273.15+theta-' num2str(Coef_mu(4),10) ')).^(-2/3).*'  num2str(Coef_mu(4)-Coef_mu(3),10) './(273.15+theta-' num2str(Coef_mu(4),10) ').^2 + 4/3*' num2str(Coef_mu(2),10) '.*((' num2str(Coef_mu(3),10) '-273.15-theta)./(273.15+theta-' num2str(Coef_mu(4),10) ')).^(1/3).*' num2str(Coef_mu(4)-Coef_mu(3),10) './(273.15+theta-' num2str(Coef_mu(4),10) ').^2)'];
                dlambda = [num2str(Coef_lambda(2),10) '+ 2*' num2str(Coef_lambda(3),10) '*(273.15+theta) + 3*' num2str(Coef_lambda(4),10) '*(273.15+theta).^2 + 4*' num2str(Coef_lambda(5),10) '*(273.15+theta).^3'];
                dcp     = [num2str(R) '/' num2str(M_l,10) '*(' num2str(Coef_cp(1),10) '/' num2str(T_crit_l,10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^(-2) -' num2str(Coef_cp(3),10) '/' num2str(T_crit_l,10) '-2*' num2str(Coef_cp(4),10) '/' num2str(T_crit_l,10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ') -3*' num2str(Coef_cp(5),10) '/' num2str(T_crit_l,10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^2 -4*' num2str(Coef_cp(6),10) '/' num2str(T_crit_l,10) '*(1-(273.15+theta)/' num2str(T_crit_l,10) ').^3)'];
                dsigma  = [num2str(Coef_sigma_l(1),10) '.*(log(1 - (theta + 273.15)./' num2str(T_crit_l,10) ').*(1 - (theta + 273.15)./' num2str(T_crit_l,10) ').^(' num2str(Coef_sigma_l(2),10) '+ (' num2str(Coef_sigma_l(3),10) '.*(theta + 273.15))./' num2str(T_crit_l,10) '+ (' num2str(Coef_sigma_l(4),10) '.*(20.*theta + 5463).^2)./(400.*' num2str(T_crit_l^2) ') + (' num2str(Coef_sigma_l(5),10) '.*(20.*theta + 5463).^3)./(8000.*' num2str(T_crit_l^3) ')).*(' num2str(Coef_sigma_l(3),10) './' num2str(T_crit_l,10) '+ (' num2str(Coef_sigma_l(4),10) '.*(800.*theta + 218520))./(400.*' num2str(T_crit_l^2) ') + (3.*' num2str(Coef_sigma_l(5),10) '.*(20.*theta + 5463).^2)./(400.*' num2str(T_crit_l^3) ')) - ((1 - (theta + 273.15)./' num2str(T_crit_l,10) ').^(' num2str(Coef_sigma_l(2),10) '+ (' num2str(Coef_sigma_l(3),10) '.*(theta + 273.15))./' num2str(T_crit_l,10) '+ (' num2str(Coef_sigma_l(4),10) '.*(20.*theta + 5463).^2)./(400.*' num2str(T_crit_l^2) ') + (' num2str(Coef_sigma_l(5),10) '.*(20.*theta + 5463).^3)./(8000.*' num2str(T_crit_l^3) ') - 1).*(' num2str(Coef_sigma_l(2),10) '+ (' num2str(Coef_sigma_l(3),10) '.*(theta + 273.15))./' num2str(T_crit_l,10) '+ (' num2str(Coef_sigma_l(4),10) '.*(20.*theta + 5463).^2)./(400.*' num2str(T_crit_l^2) ') + (' num2str(Coef_sigma_l(5),10) '.*(20.*theta + 5463).^3)./(8000.*' num2str(T_crit_l^3) ')))./' num2str(T_crit_l,10) ')'];

                ddcp    = [num2str(R) '/' num2str(M_l,10) '.*(1.0./' num2str(T_crit_l,10) '.^2.*' num2str(Coef_cp(4),10) '.*2.0-1.0./' num2str(T_crit_l,10) '.^2.*' num2str(Coef_cp(1),10) '.*1.0./((theta+2.7315e+2)./' num2str(T_crit_l,10) '-1.0).^3.*2.0+1.0./' num2str(T_crit_l,10) '.^2.*' num2str(Coef_cp(6),10) '.*((theta+2.7315e+2)./' num2str(T_crit_l,10) '-1.0).^2.*1.2e+1-1.0./' num2str(T_crit_l,10) '.^2.*' num2str(Coef_cp(5),10) '.*((theta+2.7315e+2)./' num2str(T_crit_l,10) '-1.0).*6.0)'];
                
            end
    
    else % gases
            load('Coef_gas.mat')
            load('caloric_critical_data.mat', 'M')
            Coef_mu     = Coef_mu_gas.(genvarname(substance)).'.*[10^-5 10^-7 10^-10 10^-12 10^-15];
            Coef_lambda = Coef_lambda_gas.(genvarname(substance)).'.*[10^-3 10^-3 10^-6  10^-9 10^-12];
            Coef_cp     = Coef_cp_gas.(genvarname(substance)).';
            M_g = M.(genvarname(substance))*10^-3;

            rho    = [num2str(p_fluid*M_g,10) '/' num2str(R) './(273.15+theta)'];
            mu     = [num2str(Coef_mu(1),10)     '+' num2str(Coef_mu(2),10)     '*(273.15+theta) +' num2str(Coef_mu(3),10)     '*(273.15+theta).^2 +' num2str(Coef_mu(4),10)     '*(273.15+theta).^3 +' num2str(Coef_mu(5),10)     '*(273.15+theta).^4'];
            lambda = [num2str(Coef_lambda(1),10) '+' num2str(Coef_lambda(2),10) '*(273.15+theta) +' num2str(Coef_lambda(3),10) '*(273.15+theta).^2 +' num2str(Coef_lambda(4),10) '*(273.15+theta).^3 +' num2str(Coef_lambda(5),10) '*(273.15+theta).^4'];
            cp     = [num2str(R) '/' num2str(M_g,10) '*(' num2str(Coef_cp(2),10) '+' num2str((Coef_cp(3)-Coef_cp(2)),10) '*((273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta)).^2.*(1-' num2str(Coef_cp(1),10) './(' num2str(Coef_cp(1),10) '+273.15+theta).*('  num2str(Coef_cp(4),10) '+' num2str(Coef_cp(5),10) '.*(273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta)+' num2str(Coef_cp(6),10) '.*((273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta)).^2+' num2str(Coef_cp(7),10) '.*((273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta)).^3)))'];
            
            drho    = [num2str(-p_fluid*M_g,10) '/' num2str(R) './(273.15+theta).^2'];
            dmu     = [num2str(Coef_mu(2),10)     '+ 2*' num2str(Coef_mu(3),10)     '*(273.15+theta) + 3*' num2str(Coef_mu(4),10)     '*(273.15+theta).^2 + 4*' num2str(Coef_mu(5),10)     '*(273.15+theta).^3'];
            dlambda = [num2str(Coef_lambda(2),10) '+ 2*' num2str(Coef_lambda(3),10) '*(273.15+theta) + 3*' num2str(Coef_lambda(4),10) '*(273.15+theta).^2 + 4*' num2str(Coef_lambda(5),10) '*(273.15+theta).^3'];
            dcp     = [num2str(R) '/' num2str(M_g,10) '*(' num2str(2*Coef_cp(1)*(Coef_cp(3)-Coef_cp(2)),10) '.*(273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta).^3.*(1-' num2str(Coef_cp(1),10) './(' num2str(Coef_cp(1),10) '+273.15+theta).*('  num2str(Coef_cp(4),10) '+' num2str(Coef_cp(5),10) '.*(273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta)+' num2str(Coef_cp(6),10) '.*((273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta)).^2+' num2str(Coef_cp(7),10) '.*((273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta)).^3))     +     ' num2str(Coef_cp(1).*(Coef_cp(3)-Coef_cp(2)),10) '.*(273.15+theta).^2./(' num2str(Coef_cp(1),10) '+273.15+theta).^4.*( ('  num2str(Coef_cp(4),10) '+' num2str(Coef_cp(5),10) '.*(273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta)+' num2str(Coef_cp(6),10) '.*((273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta)).^2+' num2str(Coef_cp(7),10) '.*((273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta)).^3)-' num2str(Coef_cp(1),10) './(' num2str(Coef_cp(1),10) '+273.15+theta).*(' num2str(Coef_cp(5),10) '+' num2str(2*Coef_cp(6),10) '.*((273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta))+' num2str(3*Coef_cp(7),10) '.*((273.15+theta)./(' num2str(Coef_cp(1),10) '+273.15+theta)).^2) ) )'];

            ddcp    = [num2str(R) '/' num2str(M_g,10) '.*(' num2str(Coef_cp(2)-Coef_cp(3),10) '.*((' num2str(Coef_cp(1),10) '.*(' num2str(Coef_cp(4),10) '+(' num2str(Coef_cp(5),10) '.*(theta+2.7315e+2))./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)+' num2str(Coef_cp(6),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2+' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^3.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3))./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)-1.0).*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2.*2.0+(' num2str(Coef_cp(2),10) '-' num2str(Coef_cp(3),10) ').*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2.*((' num2str(Coef_cp(1),10) '.*(' num2str(Coef_cp(5),10) '.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2.*-2.0+' num2str(Coef_cp(6),10) '.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2.*2.0+' num2str(Coef_cp(5),10) '.*(theta+2.7315e+2).*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3.*2.0-' num2str(Coef_cp(6),10) '.*(theta.*2.0+5.463e+2).*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3.*4.0+' num2str(Coef_cp(6),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^4.*6.0+' num2str(Coef_cp(7),10) '.*(theta.*2.0+5.463e+2).*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3.*3.0-' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^4.*1.8e+1+' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^3.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^5.*1.2e+1))./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)+' num2str(Coef_cp(1),10) '.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3.*(' num2str(Coef_cp(4),10) '+(' num2str(Coef_cp(5),10) '.*(theta+2.7315e+2))./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)+' num2str(Coef_cp(6),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2+' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^3.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3).*2.0-' num2str(Coef_cp(1),10) '.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2.*(' num2str(Coef_cp(5),10) './(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)-' num2str(Coef_cp(5),10) '.*(theta+2.7315e+2).*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2+' num2str(Coef_cp(6),10) '.*(theta.*2.0+5.463e+2).*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2-' num2str(Coef_cp(6),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3.*2.0+' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3.*3.0-' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^3.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^4.*3.0).*2.0)-(theta.*2.0+5.463e+2).*(' num2str(Coef_cp(2),10) '-' num2str(Coef_cp(3),10) ').*((' num2str(Coef_cp(1),10) '.*(' num2str(Coef_cp(4),10) '+(' num2str(Coef_cp(5),10) '.*(theta+2.7315e+2))./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)+' num2str(Coef_cp(6),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2+' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^3.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3))./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)-1.0).*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3.*4.0+(' num2str(Coef_cp(2),10) '-' num2str(Coef_cp(3),10) ').*((' num2str(Coef_cp(1),10) '.*(' num2str(Coef_cp(4),10) '+(' num2str(Coef_cp(5),10) '.*(theta+2.7315e+2))./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)+' num2str(Coef_cp(6),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2+' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^3.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3))./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)-1.0).*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^4.*6.0-(' num2str(Coef_cp(1),10) '.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2.*(' num2str(Coef_cp(4),10) '+(' num2str(Coef_cp(5),10) '.*(theta+2.7315e+2))./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)+' num2str(Coef_cp(6),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2+' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^3.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3)-(' num2str(Coef_cp(1),10) '.*(' num2str(Coef_cp(5),10) './(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)-' num2str(Coef_cp(5),10) '.*(theta+2.7315e+2).*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2+' num2str(Coef_cp(6),10) '.*(theta.*2.0+5.463e+2).*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2-' num2str(Coef_cp(6),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3.*2.0+' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3.*3.0-' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^3.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^4.*3.0))./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)).*(theta.*2.0+5.463e+2).*(' num2str(Coef_cp(2),10) '-' num2str(Coef_cp(3),10) ').*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2.*2.0+(' num2str(Coef_cp(1),10) '.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2.*(' num2str(Coef_cp(4),10) '+(' num2str(Coef_cp(5),10) '.*(theta+2.7315e+2))./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)+' num2str(Coef_cp(6),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2+' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^3.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3)-(' num2str(Coef_cp(1),10) '.*(' num2str(Coef_cp(5),10) './(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)-' num2str(Coef_cp(5),10) '.*(theta+2.7315e+2).*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2+' num2str(Coef_cp(6),10) '.*(theta.*2.0+5.463e+2).*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^2-' num2str(Coef_cp(6),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3.*2.0+' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3.*3.0-' num2str(Coef_cp(7),10) '.*(theta+2.7315e+2).^3.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^4.*3.0))./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2)).*(' num2str(Coef_cp(2),10) '-' num2str(Coef_cp(3),10) ').*(theta+2.7315e+2).^2.*1.0./(' num2str(Coef_cp(1),10) '+theta+2.7315e+2).^3.*4.0)'];
        
    end
end
    
if ~isempty(T_fluid)
    rho    = num2str(eval(rho),16);
    mu     = num2str(eval(mu),16);
    lambda = num2str(eval(lambda),16);
    cp     = num2str(eval(cp),16);
    dmu = '0.0'; drho = '0.0'; dlambda = '0.0'; dcp = '0.0'; ddcp = '0.0';
end

if strcmp(fluid,'liquid')
    if ~isempty(T_fluid)
        sigma  = num2str(eval(sigma),16);
        dsigma = num2str(eval(dsigma),16);
        sigma  = [sigma '+' dsigma '*(theta-' num2str(T_fluid) ')'];
    end
    
    b.sigma  = ['@(theta)' sigma];
    b.dsigma = ['@(theta)' dsigma];
end

b.mu      = ['@(theta)' mu];
b.rho     = ['@(theta)' rho];
b.lambda  = ['@(theta)' lambda];
b.cp      = ['@(theta)' cp];
b.dmu     = ['@(theta)' dmu];
b.drho    = ['@(theta)' drho];
b.dlambda = ['@(theta)' dlambda];
b.dcp     = ['@(theta)' dcp];
b.ddcp    = ['@(theta)' dcp];