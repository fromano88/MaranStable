neutral_mode_time = tic;

disp([independent ' = ' num2str(eval(independent))])
disp([dependent ' = ' num2str(eval(dependent))])
basic_state

gamma_m = [];
for m = m_start:m_end
    disp(['m = ', num2str(m)]) 
    most_dangerous_mode
    gamma_m = [gamma_m; min(real(gamma))];
    
end

tol=flowopt.tolerance*10;
    
if min(real(gamma_m))<0
    if abs(old.delta_alpha)<floor(7*pi/8/tol)*tol
        if exist('new_delta_alpha','var')==0
            delta_alpha=floor(pi/8/tol)*tol;
            
        else
            delta_alpha=new_delta_alpha;
            clear new_delta_alpha
            
        end
        
    else
        delta_alpha=floor((pi-old.delta_alpha)/tol)*tol;
        new_delta_alpha=delta_alpha;
        
    end
    while min(real(gamma_m))<0
        gammaminus_m=gamma_m;
        if exist('k','var')
            alpha_minus=alpha;
            alpha=alpha-rotdir*delta_alpha;
            %delta_alpha=delta_alpha*floor(0.682327803828019/tol)*tol;
            if abs(alpha-old.alpha)>=pi
                eval([independent '=old.' independent ';'])
                eval([dependent '=old.' dependent ';'])
                alpha=old.alpha;
                rotdir=-1*rotdir;
                neutral_mode
                return
                
            end
            eval([independent '=old.' independent '+' independent '_resolution*cos(alpha);'])
            disp([independent ' = ' num2str(eval(independent))]) 
            eval([dependent '=old.' dependent '+' dependent '_resolution*sin(alpha);'])
            disp([dependent ' = ' num2str(eval(dependent))])
            
        else
            eval(['alpha_minus=' dependent ';'])
            eval([dependent '=' dependent '/1.05;'])
            eval(['alpha=' dependent ';'])
            disp([dependent ' = ' num2str(eval(dependent))])
                        
        end
        basic_state
        gamma_m=[];
        for m=m_start:m_end
            disp(['m = ', num2str(m)]) 
            most_dangerous_mode
            gamma_m=[gamma_m; min(real(gamma))];
        end
    end
    alpha_plus_m=alpha;
    gammaplus_m=gamma_m;
    [gamma_m,m_n]=sort(real(gammaminus_m));
    gamma_m=gamma_m(real(gamma_m)<0);
    gammaplus_m=gammaplus_m(m_n);
    gammaplus_m=gammaplus_m(real(gamma_m)<0);
        
else
    if abs(old.delta_alpha)<floor(7*pi/8/tol)*tol
        if exist('new_delta_alpha','var')==0
            delta_alpha=floor(pi/8/tol)*tol;
            
        else
            delta_alpha=new_delta_alpha;
            clear new_delta_alpha
            
        end
        
    else
        delta_alpha=floor((pi-old.delta_alpha)/tol)*tol;
        new_delta_alpha=delta_alpha;
        
    end
    while min(real(gamma_m))>0
        gammaplus_m=gamma_m;
        if exist('k','var')
            alpha_plus_m=alpha;
            alpha=alpha+rotdir*delta_alpha;
            %delta_alpha=delta_alpha*floor(0.682327803828019/tol)*tol;
            if abs(alpha-old.alpha)>=pi
                eval([independent '=old.' independent ';'])
                eval([dependent '=old.' dependent ';'])
                alpha=old.alpha;
                rotdir=-1*rotdir;
                neutral_mode
                return
                
            end 
            eval([independent '=old.' independent '+' independent '_resolution*cos(alpha);'])
            disp([independent ' = ' num2str(eval(independent))]) 
            eval([dependent '=old.' dependent '+' dependent '_resolution*sin(alpha);'])
            disp([dependent ' = ' num2str(eval(dependent))])
            
        else
            eval(['alpha_plus_m=' dependent ';'])
            eval([dependent '=' dependent '*1.05;'])
            eval(['alpha=' dependent ';'])
            disp([dependent ' = ' num2str(eval(dependent))])
                        
        end
        basic_state
        gamma_m=[];
        for m=m_start:m_end
            disp(['m = ', num2str(m)]) 
            most_dangerous_mode
            gamma_m=[gamma_m; min(real(gamma))];
        end
    end
    alpha_minus=alpha;
    [gamma_m,m_n]=sort(real(gamma_m));
    gamma_m=gamma_m(real(gamma_m)<0);
    gammaplus_m=gammaplus_m(m_n);
    gammaplus_m=gammaplus_m(real(gamma_m)<0);
    
end
eval([independent '_' dependent '_gamma_m=[];'])
while ~isempty(gamma_m)
    alpha_plus=alpha_plus_m;
    m=m_n(1)+m_start-1;
    disp(['m = ', num2str(m)])   
    gammaminus=gamma_m(1);
    gammaplus=gammaplus_m(1);
    gammaplus_m=gammaplus_m(2:end);
    alpha_n=alpha_minus;
    gamman=gammaminus;
    while (abs(cos(alpha_plus)-cos(alpha_minus))>=tol || abs(sin(alpha_plus)-sin(alpha_minus))>=tol) && min(abs(real(gamma)))>=tol
        %Gottlieb10
        D=(alpha_minus-alpha_plus)/2;
        alpha_n=(alpha_minus+alpha_plus)/2;
        alpha =alpha_n;
        if exist('k','var')
            disp(['alpha = ' num2str(alpha)])
            eval([independent '=old.' independent '+' independent '_resolution*cos(alpha);'])
            disp([independent ' = ' num2str(eval(independent))]) 
            eval([dependent '=old.' dependent '+' dependent '_resolution*sin(alpha);'])
            disp([dependent ' = ' num2str(eval(dependent))])
            
        else
            eval([dependent '=alpha;'])
            disp([dependent ' = ' num2str(eval(dependent))])
                        
        end
        basic_state
        most_dangerous_mode
        if (abs(cos(alpha_plus)-cos(alpha_minus))<=tol && abs(sin(alpha_plus)-sin(alpha_minus))<=tol) || min(abs(real(gamma)))<=tol, break, end
        gamman=min(real(gamma));
        a=(gammaminus+gammaplus-2*gamman)/(2*D^2);
        b=(gammaminus-gammaplus)/(2*D);
        alpha =alpha_n-2*gamman/(b*(1+sqrt(1-4*a*gamman/b^2)));
        if exist('k','var')
            disp(['alpha = ' num2str(alpha)])
            eval([independent '=old.' independent '+' independent '_resolution*cos(alpha);'])
            disp([independent ' = ' num2str(eval(independent))]) 
            eval([dependent '=old.' dependent '+' dependent '_resolution*sin(alpha);'])
            disp([dependent ' = ' num2str(eval(dependent))])
            
        else
            eval([dependent '=alpha;'])
            disp([dependent ' = ' num2str(eval(dependent))])
                        
        end
        basic_state
        most_dangerous_mode
        if min(real(gamma))<0
            alpha_minus=alpha;
            gammaminus=min(real(gamma));
            if gamman>0
                alpha_plus=alpha_n;
                gammaplus=gamman;
            end
        else
            alpha_plus=alpha;
            gammaplus=min(real(gamma));
            if gamman<0
                alpha_minus=alpha_n;
                gammaminus=gamman;
            end
        end
    end
    eval([independent '_' dependent '_gamma_m=[' independent '_' dependent '_gamma_m;' independent ' ' dependent ' gamma(1,1) m];'])
    eigenvector_m=[eigenvector; m];
    x_m=[x; m];
    if length(m_n)==1, break, end
    m_n=m_n(2:length(real(gamma_m)<0));
    gamma_m=[];
    for n=1:length(m_n)
        m=m_n(n)+m_start-1;
        disp(['m = ', num2str(m)]) 
        most_dangerous_mode
        gamma_m=[gamma_m; min(real(gamma))];
    end
    [gamma_m,index_gamma_m]=sort(real(gamma_m));
    gamma_m=gamma_m(real(gamma_m)<0);
    m_n=m_n(index_gamma_m);
    m_n=m_n(real(gamma_m)<0);
    gammaplus_m=gammaplus_m(index_gamma_m);
    gammaplus_m=gammaplus_m(real(gamma_m)<0);
end

neutral_mode_time = toc(neutral_mode_time);
print_time(neutral_mode_time,'neutral mode time')
