function [lz, lz_inter, eta_z, eta_inter_z] = segmentation(geomz, z)
% if     z.sf='tanh' -> tanh streching function
% elseif z.sf='gp' -> geometric progression streching function
    l=1;
    z.delta_start=z.delta_start./(geomz(2:end)-geomz(1:end-1));
    z.delta_fit=z.delta_fit./(geomz(2:end)-geomz(1:end-1));
    z.delta_end=z.delta_end./(geomz(2:end)-geomz(1:end-1));
            
    for n=1:length(geomz)-1
        %l=geomz(n+1)-geomz(n);
        
        if l<=0 || z.delta_start(n)<=0 || z.delta_fit(n)<=0 || z.delta_end(n)<=0 || z.f_start(n)<=0 || z.f_end(n)<=0

                    error(['########## Grid parameters have to be larger than 0 ' num2str(n) '!!! ########## -> l(' num2str(n) ')=' num2str(l) ', delta_start(' num2str(n) ')=' num2str(z.delta_start(n)) ', delta_fit(' num2str(n) ')=' num2str(z.delta_fit(n)) ', delta_end(' num2str(n) ')=' num2str(z.delta_end(n)) ', f_start(' num2str(n) ')=' num2str(z.f_start(n)) ', f_end(' num2str(n) ')=' num2str(z.f_end(n))])
        end
        
        if z.delta_start(n)==z.delta_end(n) && z.delta_start(n)==z.delta_fit(n) && z.f_start(n)==1 &&  z.f_end(n)==1
            z.n_start(n)=0;
            z.n_fit(n)=round(l/z.delta_fit(n));
            z.n_end(n)=0;
            z.delta_start(n)=l/z.n_fit(n);
            z.delta_fit(n)=l/z.n_fit(n);
            z.delta_end(n)=l/z.n_fit(n);
            l_start(n)=0;
            l_fit(n)=l;
            l_end(n)=0;
            slope_start(n)=z.delta_start(n);
            slope_end(n)=slope_start(n);
                
        elseif strcmp(z.sf(n),'tanh')==1 && z.f_start(n)~=1 && z.f_end(n)~=1
            alpha_start(n)=acosh(sqrt(z.delta_fit(n)/z.delta_start(n)));
            x = max(-(2*alpha_start(n))/log((2*z.f_start(n) + exp(2*alpha_start(n)) + exp(2*alpha_start(n))*((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + ((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + 1)/(exp(2*alpha_start(n)) - exp(2*alpha_start(n))*((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + 2*z.f_start(n)*exp(2*alpha_start(n)) - ((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + 1)), -(2*alpha_start(n))/log((2*z.f_start(n) + exp(2*alpha_start(n)) - exp(2*alpha_start(n))*((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) - ((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + 1)/(exp(2*alpha_start(n)) + exp(2*alpha_start(n))*((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + 2*z.f_start(n)*exp(2*alpha_start(n)) + ((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + 1)));
            z.n_start(n)=round(x(1));
            alpha_end(n)=acosh(sqrt(z.delta_fit(n)/z.delta_end(n)));
            x = max(-(2*alpha_end(n))/log((2*z.f_end(n) + exp(2*alpha_end(n)) + exp(2*alpha_end(n))*((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + ((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + 1)/(exp(2*alpha_end(n)) - exp(2*alpha_end(n))*((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + 2*z.f_end(n)*exp(2*alpha_end(n)) - ((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + 1)), -(2*alpha_end(n))/log((2*z.f_end(n) + exp(2*alpha_end(n)) - exp(2*alpha_end(n))*((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) - ((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + 1)/(exp(2*alpha_end(n)) + exp(2*alpha_end(n))*((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + 2*z.f_end(n)*exp(2*alpha_end(n)) + ((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + 1)));
            z.n_end(n)=round(x(1));
            l_start(n)=z.delta_start(n)/alpha_start(n)*z.n_start(n)*tanh(alpha_start(n));
            l_end(n)=z.delta_end(n)/alpha_end(n)*z.n_end(n)*tanh(alpha_end(n));
            l_fit(n)=l-l_start(n)-l_end(n);
            z.n_fit(n)=l_fit(n)/z.delta_fit(n);
            
            x = fsolve(@(x)[l-x(1)-x(2)-x(3);
                        z.delta_start(n)-x(1)*(1-tanh(x(4)*(1-1/x(6)))/tanh(x(4)));
                        z.delta_end(n)-x(3)+x(3)*tanh(x(5)*((x(8)-1)/x(8)))/tanh(x(5));
                        z.delta_fit(n)-x(1)*x(4)/x(6)/tanh(x(4));
                        z.delta_fit(n)-x(3)*x(5)/x(8)/tanh(x(5));
                        z.delta_fit(n)-x(2)/x(7);
                        z.f_start(n)-(tanh(x(4)*(1-1/x(6)))-tanh(x(4)*(1-2/x(6))))/(tanh(x(4))-tanh(x(4)*(1-1/x(6))));
                        z.f_end(n)-(tanh(x(5)*((x(8)-1)/x(8)))-tanh(x(5)*((x(8)-2)/x(8))))/(tanh(x(5))-tanh(x(5)*((x(8)-1)/x(8))))],[l_start(n);l_fit(n);l_end(n);alpha_start(n);alpha_end(n);z.n_start(n);z.n_fit(n);z.n_end(n)],optimset('Display','off'));

            l_start(n)=x(1);
            l_fit(n)=x(2);
            l_end(n)=x(3);
            alpha_start(n)=x(4);
            alpha_end(n)=x(5);
            z.n_start(n)=round(x(6));
            z.n_fit(n)=round(x(7));
            z.n_end(n)=round(x(8));
            slope_start(n)=l_start(n)*2/sinh(2*alpha_start(n))*alpha_start(n)/z.n_start(n);
            slope_end(n)=l_end(n)*2/sinh(2*alpha_end(n))*alpha_end(n)/z.n_end(n);

            if l_start(n)<=0 || l_fit(n)<=0 || l_end(n)<=0 || z.n_start(n)<1 || z.n_fit(n)<1 || z.n_end(n)<1 || alpha_start(n)>10 || alpha_end(n)>10 || max(abs(imag(x)))>0

                error(['########## Change grid parameters of line ' num2str(n) '!!! ########## -> n_fit(' num2str(n) ')=' num2str(z.n_fit(n)) ', n_start(' num2str(n) ')=' num2str(z.n_start(n)) ', n_end(' num2str(n) ')=' num2str(z.n_end(n)) ', delta_start(' num2str(n) ')=' num2str(z.delta_start(n)) ', delta_fit(' num2str(n) ')=' num2str(z.delta_fit(n)) ', delta_end(' num2str(n) ')=' num2str(z.delta_end(n)) ', l_start(' num2str(n) ')=' num2str(l_start(n)) ', l_fit(' num2str(n) ')=' num2str(l_fit(n)) ', l_end(' num2str(n) ')=' num2str(l_end(n)) ', alpha_start(' num2str(n) ')=' num2str(alpha_start(n)) ', alpha_end(' num2str(n) ')=' num2str(alpha_end(n))])

            end
            
            
        elseif strcmp(z.sf(n),'tanh')==1 && z.f_start(n)==1 && z.f_end(n)~=1 && z.delta_start(n)==z.delta_fit(n)
            alpha_end(n)=acosh(sqrt(z.delta_fit(n)/z.delta_end(n)));
            x = max(-(2*alpha_end(n))/log((2*z.f_end(n) + exp(2*alpha_end(n)) + exp(2*alpha_end(n))*((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + ((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + 1)/(exp(2*alpha_end(n)) - exp(2*alpha_end(n))*((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + 2*z.f_end(n)*exp(2*alpha_end(n)) - ((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + 1)), -(2*alpha_end(n))/log((2*z.f_end(n) + exp(2*alpha_end(n)) - exp(2*alpha_end(n))*((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) - ((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + 1)/(exp(2*alpha_end(n)) + exp(2*alpha_end(n))*((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + 2*z.f_end(n)*exp(2*alpha_end(n)) + ((- 4*exp(2*alpha_end(n))*z.f_end(n)^2 + 2*exp(2*alpha_end(n)) + exp(4*alpha_end(n)) + 1)/(exp(2*alpha_end(n)) + 1)^2)^(1/2) + 1)));
            z.n_end(n)=round(x(1));
            l_end(n)=z.delta_end(n)/alpha_end(n)*z.n_end(n)*tanh(alpha_end(n));
            l_start(n)=l-l_end(n);
            z.n_start(n)=l_start(n)/z.delta_start(n);
            
            x = fsolve(@(x)[l-x(1)-x(2);
                        z.delta_end(n)-x(2)+x(2)*tanh(x(3)*((x(5)-1)/x(5)))/tanh(x(3));
                        z.delta_start(n)-x(2)*x(3)/x(5)/tanh(x(3));
                        z.delta_start(n)-x(1)/x(4);
                        z.f_end(n)-(tanh(x(3)*((x(5)-1)/x(5)))-tanh(x(3)*((x(5)-2)/x(5))))/(tanh(x(3))-tanh(x(3)*((x(5)-1)/x(5))))],[l_start(n);l_end(n);alpha_end(n);z.n_start(n);z.n_end(n)],optimset('Display','off'));

            l_start(n)=x(1);
            l_end(n)=x(2);
            alpha_end(n)=x(3);
            z.n_start(n)=round(x(4));
            z.n_end(n)=round(x(5));
            slope_start(n)=z.delta_start(n);
            slope_end(n)=l_end(n)*2/sinh(2*alpha_end(n))*alpha_end(n)/z.n_end(n);
            z.n_fit(n)=0;
            l_fit(n)=0;
            
            if l_start(n)<=0 || l_end(n)<=0 || z.n_start(n)<1 || z.n_end(n)<1 || alpha_end(n)>10 || max(abs(imag(x)))>0

                error(['########## Change grid parameters of line ' num2str(n) '!!! ########## -> n_start(' num2str(n) ')=' num2str(z.n_start(n)) ', n_end(' num2str(n) ')=' num2str(z.n_end(n)) ', delta_start(' num2str(n) ')=' num2str(z.delta_start(n)) ', delta_end(' num2str(n) ')=' num2str(z.delta_end(n)) ', l_start(' num2str(n) ')=' num2str(l_start(n)) ', l_end(' num2str(n) ')=' num2str(l_end(n)) ', alpha_end(' num2str(n) ')=' num2str(alpha_end(n))])

            end
                        
        elseif strcmp(z.sf(n),'tanh')==1 && z.f_start(n)~=1 && z.f_end(n)==1 && z.delta_end(n)==z.delta_fit(n)
            alpha_start(n)=acosh(sqrt(z.delta_fit(n)/z.delta_start(n)));
            x = max(-(2*alpha_start(n))/log((2*z.f_start(n) + exp(2*alpha_start(n)) + exp(2*alpha_start(n))*((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + ((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + 1)/(exp(2*alpha_start(n)) - exp(2*alpha_start(n))*((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + 2*z.f_start(n)*exp(2*alpha_start(n)) - ((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + 1)), -(2*alpha_start(n))/log((2*z.f_start(n) + exp(2*alpha_start(n)) - exp(2*alpha_start(n))*((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) - ((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + 1)/(exp(2*alpha_start(n)) + exp(2*alpha_start(n))*((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + 2*z.f_start(n)*exp(2*alpha_start(n)) + ((- 4*exp(2*alpha_start(n))*z.f_start(n)^2 + 2*exp(2*alpha_start(n)) + exp(4*alpha_start(n)) + 1)/(exp(2*alpha_start(n)) + 1)^2)^(1/2) + 1)));
            z.n_start(n)=round(x(1));
            l_start(n)=z.delta_start(n)/alpha_start(n)*z.n_start(n)*tanh(alpha_start(n));
            l_end(n)=l-l_start(n);
            z.n_end(n)=l_end(n)/z.delta_end(n);
            
            x = fsolve(@(x)[l-x(1)-x(2);
                        z.delta_start(n)-x(1)*(1-tanh(x(3)*(1-1/x(4)))/tanh(x(3)));
                        z.delta_end(n)-x(1)*x(3)/x(4)/tanh(x(3));
                        z.delta_end(n)-x(2)/x(5);
                        z.f_start(n)-(tanh(x(3)*(1-1/x(4)))-tanh(x(3)*(1-2/x(4))))/(tanh(x(3))-tanh(x(3)*(1-1/x(4))))],[l_start(n);l_end(n);alpha_start(n);z.n_start(n);z.n_end(n)],optimset('Display','off'));

            l_start(n)=x(1);
            l_end(n)=x(2);
            alpha_start(n)=x(3);
            z.n_start(n)=round(x(4));
            z.n_end(n)=round(x(5));
            slope_start(n)=l_start(n)*2/sinh(2*alpha_start(n))*alpha_start(n)/z.n_start(n);
            slope_end(n)=z.delta_end(n);
            z.n_fit(n)=0;
            l_fit(n)=0;
            
            if l_start(n)<=0 || l_end(n)<=0 || z.n_start(n)<1 || z.n_end(n)<1 || alpha_start(n)>10 || max(abs(imag(x)))>0

                error(['########## Change grid parameters of line ' num2str(n) '!!! ########## -> n_start(' num2str(n) ')=' num2str(z.n_start(n)) ', n_end(' num2str(n) ')=' num2str(z.n_end(n)) ', delta_start(' num2str(n) ')=' num2str(z.delta_start(n)) ', delta_end(' num2str(n) ')=' num2str(z.delta_end(n)) ', l_start(' num2str(n) ')=' num2str(l_start(n)) ', l_end(' num2str(n) ')=' num2str(l_end(n)) ', alpha_start(' num2str(n) ')=' num2str(alpha_start(n))])

            end
                        
        elseif strcmp(z.sf(n),'gp')==1 && z.f_start(n)~=1 && z.f_end(n)~=1
            z.n_start(n)=round(log(z.delta_fit(n)/z.delta_start(n))/log(z.f_start(n)));
            z.n_end(n)=round(log(z.delta_fit(n)/z.delta_end(n))/log(z.f_end(n)));
            l_start(n)=z.delta_start(n)*(z.f_start(n)^z.n_start(n)-1)/(z.f_start(n)-1);
            l_end(n)=z.delta_end(n)*(z.f_end(n)^z.n_end(n)-1)/(z.f_end(n)-1);
            l_fit(n)=l-l_start(n)-l_end(n);
            z.n_fit(n)=round(l_fit(n)/z.delta_fit(n));

            if l_start(n)<=0 || l_fit(n)<=0 || l_end(n)<=0 || z.n_start(n)<1 || z.n_fit(n)<1 || z.n_end(n)<1

                error(['########## Change grid parameters of line ' num2str(n) '!!! ########## -> n_fit(' num2str(n) ')=' num2str(z.n_fit(n)) ', n_start(' num2str(n) ')=' num2str(z.n_start(n)) ', n_end(' num2str(n) ')=' num2str(z.n_end(n)) ', l_start(' num2str(n) ')=' num2str(l_start(n)) ', l_fit(' num2str(n) ')=' num2str(l_fit(n)) ', l_end(' num2str(n) ')=' num2str(l_end(n))])

            end
            x = fsolve(@(x)[l-x(1)-x(2)-x(3);
                            z.delta_start(n)-x(1)*(z.f_start(n)-1)/(z.f_start(n)^x(4)-1);
                            z.delta_end(n)-x(3)*(z.f_end(n)-1)/(z.f_end(n)^x(6)-1);
                            z.delta_fit(n)-x(1)*z.f_start(n)^x(4)*log(z.f_start(n))/(z.f_start(n)^x(4)-1);
                            z.delta_fit(n)-x(3)*z.f_end(n)^x(6)*log(z.f_end(n))/(z.f_end(n)^x(6)-1);
                            z.delta_fit(n)-x(2)/x(5)],[l_start(n);l_fit(n);l_end(n);z.n_start(n);z.n_fit(n);z.n_end(n)],optimset('Display','off','TolFun',eps,'TolX',eps));

            l_start(n)=x(1);
            l_fit(n)=x(2);
            l_end(n)=x(3);
            z.n_start(n)=round(x(4));
            z.n_fit(n)=round(x(5));
            z.n_end(n)=round(x(6));
            slope_start(n)=l_start(n)*log(z.f_start(n))/(z.f_start(n)^z.n_start(n)-1);
            slope_end(n)=l_end(n)*log(z.f_end(n))/(z.f_end(n)^z.n_end(n)-1);
            
        elseif strcmp(z.sf(n),'gp')==1 && z.f_start(n)==1 && z.f_end(n)~=1 && z.delta_start(n)==z.delta_fit(n)
            z.n_end(n)=round(log(z.delta_fit(n)/z.delta_end(n))/log(z.f_end(n)));
            l_end(n)=z.delta_end(n)*(z.f_end(n)^z.n_end(n)-1)/(z.f_end(n)-1);
            l_fit(n)=0;
            l_start(n)=l-l_end(n);
            z.n_start(n)=round(l_start(n)/z.delta_start(n));
            z.n_fit(n)=0;

            if l_start(n)<=0 || l_end(n)<=0 || z.n_start(n)<1 || z.n_end(n)<1

                error(['########## Change grid parameters of line ' num2str(n) '!!! ########## -> n_start(' num2str(n) ')=' num2str(z.n_start(n)) ', n_end(' num2str(n) ')=' num2str(z.n_end(n)) ', l_start(' num2str(n) ')=' num2str(l_start(n)) ', l_end(' num2str(n) ')=' num2str(l_end(n))])

            end
            x = fsolve(@(x)[l-x(1)-x(2);
                            z.delta_end(n)-x(2)*(z.f_end(n)-1)/(z.f_end(n)^x(4)-1);
                            z.delta_start(n)-x(2)*z.f_end(n)^x(4)*log(z.f_end(n))/(z.f_end(n)^x(4)-1);
                            z.delta_start(n)-x(1)/x(3)],[l_start(n);l_end(n);z.n_start(n);z.n_end(n)],optimset('Display','off','TolFun',eps,'TolX',eps));

            l_start(n)=x(1);
            l_end(n)=x(2);
            z.n_start(n)=round(x(3));
            z.n_end(n)=round(x(4));
            slope_start(n)=z.delta_start(n);
            slope_end(n)=l_end(n)*log(z.f_end(n))/(z.f_end(n)^z.n_end(n)-1);
            
        elseif strcmp(z.sf(n),'gp')==1 && z.f_start(n)~=1 && z.f_end(n)==1 && z.delta_end(n)==z.delta_fit(n)
            z.n_start(n)=round(log(z.delta_fit(n)/z.delta_start(n))/log(z.f_start(n)));
            l_start(n)=z.delta_start(n)*(z.f_start(n)^z.n_start(n)-1)/(z.f_start(n)-1);
            l_fit(n)=0;
            l_end(n)=l-l_start(n);
            z.n_end(n)=round(l_end(n)/z.delta_end(n));
            z.n_fit(n)=0;

            if l_start(n)<=0 || l_end(n)<=0 || z.n_start(n)<1 || z.n_end(n)<1

                error(['########## Change grid parameters of line ' num2str(n) '!!! ########## -> n_start(' num2str(n) ')=' num2str(z.n_start(n)) ', n_end(' num2str(n) ')=' num2str(z.n_end(n)) ', l_start(' num2str(n) ')=' num2str(l_start(n)) ', l_end(' num2str(n) ')=' num2str(l_end(n))])

            end
            x = fsolve(@(x)[l-x(1)-x(2);
                            z.delta_start(n)-x(1)*(z.f_start(n)-1)/(z.f_start(n)^x(3)-1);
                            z.delta_end(n)-x(1)*z.f_start(n)^x(3)*log(z.f_start(n))/(z.f_start(n)^x(3)-1);
                            z.delta_end(n)-x(2)/x(4)],[l_start(n);l_end(n);z.n_start(n);z.n_end(n)],optimset('Display','off','TolFun',eps,'TolX',eps));

            l_start(n)=x(1);
            l_end(n)=x(2);
            z.n_start(n)=round(x(3));
            z.n_end(n)=round(x(4));
            slope_start(n)=l_start(n)*log(z.f_start(n))/(z.f_start(n)^z.n_start(n)-1);
            slope_end(n)=z.delta_end(n);
            
        else
            error(['########## Change grid parameters of line ' num2str(n) '!!! ##########'])
        end
    end
    
    slope=[slope_start(1), (slope_start(2:end)+slope_end(1:end-1))/2, slope_end(end)];
    for n=1:length(geomz)-1
        if z.f_start(n)==1 && z.f_end(n)==1
            slope(n)=slope_start(n);
            slope(n+1)=slope_end(n);
        end
    end
     
    for n=1:length(geomz)-1
        %l=geomz(n+1)-geomz(n);
        
        if strcmp(z.sf(n),'tanh')==1 && z.f_start(n)==1 && z.f_end(n)~=1
             x = fsolve(@(x)[l-x(1)-x(2);
                        x(5)-x(2)+x(2)*tanh(x(3)*((z.n_end(n)-1)/z.n_end(n)))/tanh(x(3));
                        x(4)-x(2)*x(3)/z.n_end(n)/tanh(x(3));
                        x(4)-x(1)/z.n_start(n);
                        x(4)-slope(n)],[l_start(n);l_end(n);alpha_end(n);z.delta_start(n);z.delta_end(n)],optimset('Display','off','TolFun',eps,'TolX',eps));
                        
            l_start(n)=x(1);
            l_end(n)=x(2);
            alpha_end(n)=x(3);
            z.delta_start(n)=x(4);
            z.delta_end(n)=x(5);
            slope(n+1)=l_end(n)*2/sinh(2*alpha_end(n))*alpha_end(n)/z.n_end(n);
            z.delta_fit(n)=z.delta_start(n);
            z.n_fit(n)=0;
            l_fit(n)=0;
            z.f_end(n)=(tanh(alpha_end(n)*((z.n_end(n)-1)/z.n_end(n)))-tanh(alpha_end(n)*((z.n_end(n)-2)/z.n_end(n))))/(tanh(alpha_end(n))-tanh(alpha_end(n)*((z.n_end(n)-1)/z.n_end(n))));
            
        elseif strcmp(z.sf(n),'tanh')==1 && z.f_start(n)~=1 && z.f_end(n)==1
            x = fsolve(@(x)[l-x(1)-x(2);
                        x(4)-x(1)*(1-tanh(x(3)*(1-1/z.n_start(n)))/tanh(x(3)));
                        x(5)-x(1)*x(3)/z.n_start(n)/tanh(x(3));
                        x(5)-x(2)/z.n_end(n);
                        x(5)-slope(n+1)],[l_start(n);l_end(n);alpha_start(n);z.delta_start(n);z.delta_end(n)],optimset('Display','off','TolFun',eps,'TolX',eps));

            l_start(n)=x(1);
            l_end(n)=x(2);
            alpha_start(n)=x(3);
            z.delta_start(n)=x(4);
            z.delta_end(n)=x(5);
            slope(n)=l_start(n)*2/sinh(2*alpha_start(n))*alpha_start(n)/z.n_start(n);
            z.delta_fit(n)=z.delta_end(n);
            z.n_fit(n)=0;
            l_fit(n)=0;
            z.f_start(n)=(tanh(alpha_start(n)*(1-1/z.n_start(n)))-tanh(alpha_start(n)*(1-2/z.n_start(n))))/(tanh(alpha_start(n))-tanh(alpha_start(n)*(1-1/z.n_start(n))));
            
        elseif strcmp(z.sf(n),'gp')==1 && z.f_start(n)==1 && z.f_end(n)~=1
            x = fsolve(@(x)[l-x(1)-x(2);
                            slope(n)-x(2)*x(5)^z.n_end(n)*log(x(5))/(x(5)^z.n_end(n)-1);%x(4)-slope(n)/(x(5)^(z.n_end(n)-1));
                            x(4)-x(2)*(x(5)-1)/(x(5)^z.n_end(n)-1);
                            slope(n)-x(1)/z.n_start(n);
                            x(3)-slope(n)],[l_start(n);l_end(n);z.delta_start(n);z.delta_end(n);z.f_end(n)],optimset('Display','off','TolFun',eps,'TolX',eps));
            
            l_start(n)=x(1);
            l_end(n)=x(2);
            z.delta_start(n)=x(3);
            z.delta_fit(n)=z.delta_start(n);
            z.delta_end(n)=x(4);
            z.f_end(n)=x(5);
            slope(n+1)=l_end(n)*log(z.f_end(n))/(z.f_end(n)^z.n_end(n)-1);
            
        elseif strcmp(z.sf(n),'gp')==1 && z.f_start(n)~=1 && z.f_end(n)==1
            x = fsolve(@(x)[l-x(1)-x(2);
                            slope(n+1)-x(1)*x(5)^z.n_start(n)*log(x(5))/(x(5)^z.n_start(n)-1);%x(3)-slope(n+1)/(x(5)^(z.n_start(n)-1));
                            x(3)-x(1)*(x(5)-1)/(x(5)^z.n_start(n)-1);
                            slope(n+1)-x(2)/z.n_end(n);
                            x(4)-slope(n+1)],[l_start(n);l_end(n);z.delta_start(n);z.delta_end(n);z.f_start(n)],optimset('Display','off','TolFun',eps,'TolX',eps));
            
            l_start(n)=x(1);
            l_end(n)=x(2);
            z.delta_start(n)=x(3);
            z.delta_end(n)=x(4);
            z.delta_fit(n)=z.delta_end(n);
            z.f_start(n)=x(5);
            slope(n)=l_start(n)*log(z.f_start(n))/(z.f_start(n)^z.n_start(n)-1);
            
        end
    end
    
    lz=geomz(1);
    lz_inter=geomz(1);
    eta_z=0;
    eta_inter_z=[];
            
    for n=1:length(geomz)-1
        %l=geomz(n+1)-geomz(n);
        eta_z(end)=[];
        
        if z.delta_start(n)==z.delta_end(n) && z.delta_start(n)==z.delta_fit(n) && z.f_start(n)==1 &&  z.f_end(n)==1
            
            lz=[lz geomz(n)+l_fit(n)*(1:z.n_fit(n)-1)/z.n_fit(n)*(geomz(n+1)-geomz(n)) geomz(n+1)];
            lz_inter=[lz_inter geomz(n)+l_fit(n)*(0.5:z.n_fit(n)-0.5)/z.n_fit(n)*(geomz(n+1)-geomz(n))];
            eta_z=[eta_z z.n_fit(n)/l_fit(n)*ones(1,z.n_fit(n)+1)/(geomz(n+1)-geomz(n))];
            eta_inter_z=[eta_inter_z z.n_fit(n)/l_fit(n)*ones(1,z.n_fit(n))/(geomz(n+1)-geomz(n))];
                        
        elseif strcmp(z.sf(n),'tanh')==1 && z.f_start(n)~=1 && z.f_end(n)~=1    
            x = fsolve(@(x)[l-x(1)-x(2)-x(3);
                        x(6)-x(1)*(1-tanh(x(4)*(1-1/z.n_start(n)))/tanh(x(4)));
                        x(8)-x(3)+x(3)*tanh(x(5)*((z.n_end(n)-1)/z.n_end(n)))/tanh(x(5));
                        x(7)-x(1)*x(4)/z.n_start(n)/tanh(x(4));
                        x(7)-x(3)*x(5)/z.n_end(n)/tanh(x(5));
                        x(7)-x(2)/z.n_fit(n);
                        slope(n)-x(1)*2/sinh(2*x(4))*x(4)/z.n_start(n);
                        slope(n+1)-x(3)*2/sinh(2*x(5))*x(5)/z.n_end(n)],[l_start(n);l_fit(n);l_end(n);alpha_start(n);alpha_end(n);z.delta_start(n);z.delta_fit(n);z.delta_end(n)],optimset('Display','off','TolFun',eps,'TolX',eps));

            l_start(n)=x(1);
            l_fit(n)=x(2);
            l_end(n)=x(3);
            alpha_start(n)=x(4);
            alpha_end(n)=x(5);
            z.delta_start(n)=x(6);
            z.delta_fit(n)=x(7);
            z.delta_end(n)=x(8);
            z.f_start(n)=(tanh(alpha_start(n)*(1-1/z.n_start(n)))-tanh(alpha_start(n)*(1-2/z.n_start(n))))/(tanh(alpha_start(n))-tanh(alpha_start(n)*(1-1/z.n_start(n))));
            z.f_end(n)=(tanh(alpha_end(n)*((z.n_end(n)-1)/z.n_end(n)))-tanh(alpha_end(n)*((z.n_end(n)-2)/z.n_end(n))))/(tanh(alpha_end(n))-tanh(alpha_end(n)*((z.n_end(n)-1)/z.n_end(n))));
            
            lz_start=l_start(n)*(1-tanh(alpha_start(n)*(1-(0:z.n_start(n))/z.n_start(n)))/tanh(alpha_start(n)));
            lz_fit=l_fit(n)*(0:z.n_fit(n))/z.n_fit(n);
            lz_end=l_end(n)*(tanh(alpha_end(n)*((0:z.n_end(n))/z.n_end(n)))/tanh(alpha_end(n)));
            %lz=[lz geomz(n)+[l_start(n)*(1-tanh(alpha_start(n)*(1-(1:z.n_start(n))/z.n_start(n)))/tanh(alpha_start(n))), l_start(n)+l_fit(n)*(1:z.n_fit(n))/z.n_fit(n), l_start(n)+l_fit(n)+l_end(n)*(tanh(alpha_end(n)*((1:z.n_end(n))/z.n_end(n)))/tanh(alpha_end(n)))]];
            lz=[lz geomz(n)+[lz_start(2:end), l_start(n)+lz_fit(2:end), l_start(n)+l_fit(n)+lz_end(2:end-1)]*(geomz(n+1)-geomz(n)) geomz(n+1)];
            
            eta_z_start=z.n_start(n)/alpha_start(n)*tanh(alpha_start(n))/l_start(n)./(1-((1-lz_start/l_start(n))*tanh(alpha_start(n))).^2);
            eta_z_fit=z.n_fit(n)/l_fit(n)*ones(1,length(lz_fit));
            eta_z_end=z.n_end(n)/alpha_end(n)*tanh(alpha_end(n))/l_end(n)./(1-((lz_end/l_end(n))*tanh(alpha_end(n))).^2);
            eta_z=[eta_z [eta_z_start eta_z_fit(2:end) eta_z_end(2:end)]/(geomz(n+1)-geomz(n))];
            
            lz_inter_start=l_start(n)*(1-tanh(alpha_start(n)*(1-(0.5:z.n_start(n)-0.5)/z.n_start(n)))/tanh(alpha_start(n)));
            lz_inter_fit=l_fit(n)*(0.5:z.n_fit(n)-0.5)/z.n_fit(n);
            lz_inter_end=l_end(n)*(tanh(alpha_end(n)*((0.5:z.n_end(n)-0.5)/z.n_end(n)))/tanh(alpha_end(n)));
            %lz_inter=[lz_inter geomz(n)+[l_start(n)*(1-tanh(alpha_start(n)*(1-(0.5:z.n_start(n)-0.5)/z.n_start(n)))/tanh(alpha_start(n))), l_start(n)+l_fit(n)*(0.5:z.n_fit(n)-0.5)/z.n_fit(n), l_start(n)+l_fit(n)+l_end(n)*(tanh(alpha_end(n)*((0.5:z.n_end(n)-0.5)/z.n_end(n)))/tanh(alpha_end(n)))]];
            lz_inter=[lz_inter geomz(n)+[lz_inter_start, l_start(n)+lz_inter_fit, l_start(n)+l_fit(n)+lz_inter_end]*(geomz(n+1)-geomz(n))];
            
            eta_inter_z_start=z.n_start(n)/alpha_start(n)*tanh(alpha_start(n))/l_start(n)./(1-((1-lz_inter_start/l_start(n))*tanh(alpha_start(n))).^2);
            eta_inter_z_fit=z.n_fit(n)/l_fit(n)*ones(1,length(lz_inter_fit));
            eta_inter_z_end=z.n_end(n)/alpha_end(n)*tanh(alpha_end(n))/l_end(n)./(1-((lz_inter_end/l_end(n))*tanh(alpha_end(n))).^2);
            eta_inter_z=[eta_inter_z [eta_inter_z_start eta_inter_z_fit eta_inter_z_end]/(geomz(n+1)-geomz(n))];
            
        elseif strcmp(z.sf(n),'tanh')==1 && z.f_start(n)==1 && z.f_end(n)~=1 && z.delta_start(n)==z.delta_fit(n)
            
            lz_start=l_start(n)*(0:z.n_start(n))/z.n_start(n);
            lz_end=l_end(n)*(tanh(alpha_end(n)*((0:z.n_end(n))/z.n_end(n)))/tanh(alpha_end(n)));
            %lz=[lz geomz(n)+[l_start(n)*(1:z.n_start(n))/z.n_start(n), l_start(n)+l_end(n)*(tanh(alpha_end(n)*((1:z.n_end(n))/z.n_end(n)))/tanh(alpha_end(n)))]];
            lz=[lz geomz(n)+[lz_start(2:end), l_start(n)+lz_end(2:end-1)]*(geomz(n+1)-geomz(n)) geomz(n+1)];
            
            eta_z_start=z.n_start(n)/l_start(n)*ones(1,length(lz_start));
            eta_z_end=z.n_end(n)/alpha_end(n)*tanh(alpha_end(n))/l_end(n)./(1-((lz_end/l_end(n))*tanh(alpha_end(n))).^2);
            eta_z=[eta_z [eta_z_start eta_z_end(2:end)]/(geomz(n+1)-geomz(n))];
            
            lz_inter_start=l_start(n)*(0.5:z.n_start(n)-0.5)/z.n_start(n);
            lz_inter_end=l_end(n)*(tanh(alpha_end(n)*((0.5:z.n_end(n)-0.5)/z.n_end(n)))/tanh(alpha_end(n)));
            %lz_inter=[lz geomz(n)+[l_start(n)*(0.5:z.n_start(n)-0.5)/z.n_start(n), l_start(n)+l_end(n)*(tanh(alpha_end(n)*((0.5:z.n_end(n)-0.5)/z.n_end(n)))/tanh(alpha_end(n)))]];
            lz_inter=[lz_inter geomz(n)+[lz_inter_start, l_start(n)+lz_inter_end]*(geomz(n+1)-geomz(n))];
            
            eta_inter_z_start=z.n_start(n)/l_start(n)*ones(1,length(lz_inter_start));
            eta_inter_z_end=z.n_end(n)/alpha_end(n)*tanh(alpha_end(n))/l_end(n)./(1-((lz_inter_end/l_end(n))*tanh(alpha_end(n))).^2);
            eta_inter_z=[eta_inter_z [eta_inter_z_start eta_inter_z_end]/(geomz(n+1)-geomz(n))];
            
        elseif strcmp(z.sf(n),'tanh')==1 && z.f_start(n)~=1 && z.f_end(n)==1 && z.delta_end(n)==z.delta_fit(n)
            
            lz_start=l_start(n)*(1-tanh(alpha_start(n)*(1-(0:z.n_start(n))/z.n_start(n)))/tanh(alpha_start(n)));
            lz_end=l_end(n)*(0:z.n_end(n))/z.n_end(n);
            %lz=[lz geomz(n)+[l_start(n)*(1-tanh(alpha_start(n)*(1-(1:z.n_start(n))/z.n_start(n)))/tanh(alpha_start(n))), l_start(n)+l_end(n)*(1:z.n_end(n))/z.n_end(n)]];
            lz=[lz geomz(n)+[lz_start(2:end), l_start(n)+lz_end(2:end-1)]*(geomz(n+1)-geomz(n)) geomz(n+1)];
            
            eta_z_start=z.n_start(n)/alpha_start(n)*tanh(alpha_start(n))/l_start(n)./(1-((1-lz_start/l_start(n))*tanh(alpha_start(n))).^2);
            eta_z_end=z.n_end(n)/l_end(n)*ones(1,length(lz_end));
            eta_z=[eta_z [eta_z_start eta_z_end(2:end)]/(geomz(n+1)-geomz(n))];
            
            lz_inter_start=l_start(n)*(1-tanh(alpha_start(n)*(1-(0.5:z.n_start(n)-0.5)/z.n_start(n)))/tanh(alpha_start(n)));
            lz_inter_end=l_end(n)*(0.5:z.n_end(n)-0.5)/z.n_end(n);
            %lz_inter=[lz geomz(n)+[l_start(n)*(1-tanh(alpha_start(n)*(1-(0.5:z.n_start(n)-0.5)/z.n_start(n)))/tanh(alpha_start(n))), l_start(n)+l_end(n)*(0.5:z.n_end(n)-0.5)/z.n_end(n)]];
            lz_inter=[lz_inter geomz(n)+[lz_inter_start, l_start(n)+lz_inter_end]*(geomz(n+1)-geomz(n))];
            
            eta_inter_z_start=z.n_start(n)/alpha_start(n)*tanh(alpha_start(n))/l_start(n)./(1-((1-lz_inter_start/l_start(n))*tanh(alpha_start(n))).^2);
            eta_inter_z_end=z.n_end(n)/l_end(n)*ones(1,length(lz_inter_end));
            eta_inter_z=[eta_inter_z [eta_inter_z_start eta_inter_z_end]/(geomz(n+1)-geomz(n))];
            
        elseif strcmp(z.sf(n),'gp')==1 && z.f_start(n)~=1 && z.f_end(n)~=1    
            x = fsolve(@(x)[l-x(1)-x(2)-x(3);
                            x(4)-x(1)*(x(7)-1)/(x(7)^z.n_start(n)-1);
                            x(6)-x(3)*(x(8)-1)/(x(8)^z.n_end(n)-1);
                            x(5)-x(1)*x(7)^z.n_start(n)*log(x(7))/(x(7)^z.n_start(n)-1);
                            x(5)-x(3)*x(8)^z.n_end(n)*log(x(8))/(x(8)^z.n_end(n)-1);
                            x(5)-x(2)/z.n_fit(n);
                            slope(n)-x(1)*log(x(7))/(x(7)^z.n_start(n)-1);
                            slope(n+1)-x(3)*log(x(8))/(x(8)^z.n_end(n)-1)],[l_start(n);l_fit(n);l_end(n);z.delta_start(n);z.delta_fit(n);z.delta_end(n);z.f_start(n);z.f_end(n)],optimset('Display','off','TolFun',eps,'TolX',eps));

            l_start(n)=x(1);
            l_fit(n)=x(2);
            l_end(n)=x(3);
            z.delta_start(n)=x(4);
            z.delta_fit(n)=x(5);
            z.delta_end(n)=x(6);
            z.f_start(n)=x(7);
            z.f_end(n)=x(8);
            
            lz_start=l_start(n)*(z.f_start(n).^(0:z.n_start(n))-1)/(z.f_start(n)^z.n_start(n)-1);
            lz_fit=l_fit(n)*(0:z.n_fit(n))/z.n_fit(n);
            lz_end=l_end(n)*(1-(z.f_end(n).^(z.n_end(n)-(0:z.n_end(n)))-1)/(z.f_end(n)^z.n_end(n)-1));
            %lz=[lz geomz(n)+[l_start(n)*(z.f_start(n).^(1:z.n_start(n))-1)/(z.f_start(n)^z.n_start(n)-1), l_start(n)+l_fit(n)*(1:z.n_fit(n))/z.n_fit(n), l_start(n)+l_fit(n)+l_end(n)*(1-(z.f_end(n).^(z.n_end(n)-(1:z.n_end(n)))-1)/(z.f_end(n)^z.n_end(n)-1))]];
            lz=[lz geomz(n)+[lz_start(2:end), l_start(n)+lz_fit(2:end), l_start(n)+l_fit(n)+lz_end(2:end-1)]*(geomz(n+1)-geomz(n)) geomz(n+1)];
            
            eta_z_start=(z.f_start(n)^z.n_start(n)-1)/l_start(n)/log(z.f_start(n))./(lz_start/l_start(n)*(z.f_start(n)^z.n_start(n)-1)+1);
            eta_z_fit=z.n_fit(n)/l_fit(n)*ones(1,length(lz_fit));
            eta_z_end=(z.f_end(n)^z.n_end(n)-1)/l_end(n)/log(z.f_end(n))./((1-lz_end/l_end(n))*(z.f_end(n)^z.n_end(n)-1)+1);
            eta_z=[eta_z [eta_z_start eta_z_fit(2:end) eta_z_end(2:end)]/(geomz(n+1)-geomz(n))];
            
            lz_inter_start=l_start(n)*(z.f_start(n).^(0.5:z.n_start(n)-0.5)-1)/(z.f_start(n)^z.n_start(n)-1);
            lz_inter_fit=l_fit(n)*(0.5:z.n_fit(n)-0.5)/z.n_fit(n);
            lz_inter_end=l_end(n)*(1-(z.f_end(n).^(z.n_end(n)-(0.5:z.n_end(n)-0.5))-1)/(z.f_end(n)^z.n_end(n)-1));
            %lz_inter=[lz_inter geomz(n)+[l_start(n)*(z.f_start(n).^(0.5:z.n_start(n)-0.5)-1)/(z.f_start(n)^z.n_start(n)-1), l_start(n)+l_fit(n)*(0.5:z.n_fit(n)-0.5)/z.n_fit(n), l_start(n)+l_fit(n)+l_end(n)*(1-(z.f_end(n).^(z.n_end(n)-(0.5:z.n_end(n)-0.5))-1)/(z.f_end(n)^z.n_end(n)-1))]];
            lz_inter=[lz_inter geomz(n)+[lz_inter_start, l_start(n)+lz_inter_fit, l_start(n)+l_fit(n)+lz_inter_end]*(geomz(n+1)-geomz(n))];
            
            eta_inter_z_start=(z.f_start(n)^z.n_start(n)-1)/l_start(n)/log(z.f_start(n))./(lz_inter_start/l_start(n)*(z.f_start(n)^z.n_start(n)-1)+1);
            eta_inter_z_fit=z.n_fit(n)/l_fit(n)*ones(1,length(lz_inter_fit));
            eta_inter_z_end=(z.f_end(n)^z.n_end(n)-1)/l_end(n)/log(z.f_end(n))./((1-lz_inter_end/l_end(n))*(z.f_end(n)^z.n_end(n)-1)+1);
            eta_inter_z=[eta_inter_z [eta_inter_z_start eta_inter_z_fit eta_inter_z_end]/(geomz(n+1)-geomz(n))];
            
        elseif strcmp(z.sf(n),'gp')==1 && z.f_start(n)==1 && z.f_end(n)~=1 && z.delta_start(n)==z.delta_fit(n)
            
            
            lz_start=l_start(n)*(0:z.n_start(n))/z.n_start(n);
            lz_end=l_end(n)*(1-(z.f_end(n).^(z.n_end(n)-(0:z.n_end(n)))-1)/(z.f_end(n)^z.n_end(n)-1));
            %lz=[lz geomz(n)+[l_start(n)*(1:z.n_start(n))/z.n_start(n), l_start(n)+l_end(n)*(1-(z.f_end(n).^(z.n_end(n)-(1:z.n_end(n)))-1)/(z.f_end(n)^z.n_end(n)-1))]];
            lz=[lz geomz(n)+[lz_start(2:end), l_start(n)+lz_end(2:end-1)]*(geomz(n+1)-geomz(n)) geomz(n+1)];
            
            eta_z_start=z.n_start(n)/l_start(n)*ones(1,length(lz_start));
            eta_z_end=(z.f_end(n)^z.n_end(n)-1)/l_end(n)/log(z.f_end(n))./((1-lz_end/l_end(n))*(z.f_end(n)^z.n_end(n)-1)+1);
            eta_z=[eta_z [eta_z_start eta_z_end(2:end)]/(geomz(n+1)-geomz(n))];
            
            lz_inter_start=l_start(n)*(0.5:z.n_start(n)-0.5)/z.n_start(n);
            lz_inter_end=l_end(n)*(1-(z.f_end(n).^(z.n_end(n)-(0.5:z.n_end(n)-0.5))-1)/(z.f_end(n)^z.n_end(n)-1));
            %lz_inter=[lz_inter geomz(n)+[l_start(n)*(0.5:z.n_start(n)-0.5)/z.n_start(n), l_start(n)+l_end(n)*(1-(z.f_end(n).^(z.n_end(n)-(0.5:z.n_end(n)-0.5))-1)/(z.f_end(n)^z.n_end(n)-1))]];
            lz_inter=[lz_inter geomz(n)+[lz_inter_start, l_start(n)+lz_inter_end]*(geomz(n+1)-geomz(n))];
            
            eta_inter_z_start=z.n_start(n)/l_start(n)*ones(1,length(lz_inter_start));
            eta_inter_z_end=(z.f_end(n)^z.n_end(n)-1)/l_end(n)/log(z.f_end(n))./((1-lz_inter_end/l_end(n))*(z.f_end(n)^z.n_end(n)-1)+1);
            eta_inter_z=[eta_inter_z [eta_inter_z_start eta_inter_z_end]/(geomz(n+1)-geomz(n))];
                        
        elseif strcmp(z.sf(n),'gp')==1 && z.f_start(n)~=1 && z.f_end(n)==1 && z.delta_end(n)==z.delta_fit(n)
            
            lz_start=l_start(n)*(z.f_start(n).^(0:z.n_start(n))-1)/(z.f_start(n)^z.n_start(n)-1);
            lz_end=l_end(n)*(0:z.n_end(n))/z.n_end(n);
            %lz=[lz geomz(n)+[l_start(n)*(z.f_start(n).^(1:z.n_start(n))-1)/(z.f_start(n)^z.n_start(n)-1), l_start(n)+l_end(n)*(1:z.n_end(n))/z.n_end(n)]];
            lz=[lz geomz(n)+[lz_start(2:end), l_start(n)+lz_end(2:end-1)]*(geomz(n+1)-geomz(n)) geomz(n+1)];
            
            eta_z_start=(z.f_start(n)^z.n_start(n)-1)/l_start(n)/log(z.f_start(n))./(lz_start/l_start(n)*(z.f_start(n)^z.n_start(n)-1)+1);
            eta_z_end=z.n_end(n)/l_end(n)*ones(1,length(lz_end));
            eta_z=[eta_z [eta_z_start eta_z_end(2:end)]/(geomz(n+1)-geomz(n))];
            
            lz_inter_start=l_start(n)*(z.f_start(n).^(0.5:z.n_start(n)-0.5)-1)/(z.f_start(n)^z.n_start(n)-1);
            lz_inter_end=l_end(n)*(0.5:z.n_end(n)-0.5)/z.n_end(n);
            %lz_inter=[lz_inter geomz(n)+[l_start(n)*(z.f_start(n).^(0.5:z.n_start(n)-0.5)-1)/(z.f_start(n)^z.n_start(n)-1), l_start(n)+l_end(n)*(0.5:z.n_end(n)-0.5)/z.n_end(n)]];
            lz_inter=[lz_inter geomz(n)+[lz_inter_start, l_start(n)+lz_inter_end]*(geomz(n+1)-geomz(n))];
            
            eta_inter_z_start=(z.f_start(n)^z.n_start(n)-1)/l_start(n)/log(z.f_start(n))./(lz_inter_start/l_start(n)*(z.f_start(n)^z.n_start(n)-1)+1);
            eta_inter_z_end=z.n_end(n)/l_end(n)*ones(1,length(lz_inter_end));
            eta_inter_z=[eta_inter_z [eta_inter_z_start eta_inter_z_end]/(geomz(n+1)-geomz(n))];
            
        end
        
    end
    
    lz_inter=[lz_inter geomz(end)];
    eta_inter_z=[eta_z(1) eta_inter_z eta_z(end)];