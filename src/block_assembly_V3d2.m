if flowopt.stability~=1
    flowopt.m=[];
    if exist('Rs','var')==0
        Rs=[];
    end
        
end
%clear x0 x_new

x0=[];
x_new=[];
vector_index=0;
    
for n=1:length(blocks)
    for o=1:length(eval([blocks{n} '.bc.z']))
        if ((eval([blocks{n} '.bc.z{' num2str(o) '}(1)'])=='s' && eval([blocks{n} '.bc.z{' num2str(o) '}(4)'])>'0') || eval([blocks{n} '.bc.z{' num2str(o) '}(1)'])=='c' || eval([blocks{n} '.bc.z{' num2str(o) '}(1)'])=='p') && exist(eval([blocks{n} '.bc.z{' num2str(o) '}']),'var')
            eval(['clear ' eval([blocks{n} '.bc.z{' num2str(o) '}'])])
            eval([blocks{n} '.' eval([blocks{n} '.bc.z{' num2str(o) '}']) '.o=0;'])

        elseif (eval([blocks{n} '.bc.z{' num2str(o) '}(1)'])=='s' || eval([blocks{n} '.bc.z{' num2str(o) '}(1)'])=='c' || eval([blocks{n} '.bc.z{' num2str(o) '}(1)'])=='p')
            eval([eval([blocks{n} '.bc.z{' num2str(o) '}']) '=1;'])
            eval([blocks{n} '.' eval([blocks{n} '.bc.z{' num2str(o) '}']) '.o=1;'])

        end

    end
    
    eval([blocks{n} '.vector_index=vector_index;'])
            
    % initial condition / residuals
        if exist('x','var')==1
            eval([blocks{n} '.Rs1=[];'])
            eval([blocks{n} '.Rsend=[];'])
            eval(['old' blocks{n} '.Rs1=[];'])
            eval(['old' blocks{n} '.Rsend=[];'])
            if eval(['max(strncmp(''ss'', {' blocks{n} '.bc.z{1:end}},2))==1'])
                if exist('old','var')==0 || isfield(old,'x')==0 || eval(['exist(''old' blocks{n} ''',''var'')'])==0
                    old.x=[];
                    eval(['old' blocks{n} '=[];'])
                    
                end
                eval(['[' blocks{n} ', old' blocks{n} ']=metric_update(' blocks{n} ', x, Rs, old, old' blocks{n} ');'])
                
            end
            
            if eval(['exist(''old' blocks{n} ''',''var'')==0 || isfield(old' blocks{n} ',''XI'')==0 || isequal(' blocks{n} '.Z,old' blocks{n} '.Z)'])
                eval([blocks{n} '.u.u = reshape(x(vector_index+1:vector_index+(' blocks{n} '.K+1)*(' blocks{n} '.J+2),1),[' blocks{n} '.J+2,' blocks{n} '.K+1]); ' blocks{n} '.u.u = ' blocks{n} '.u.u.'';'])
                eval([blocks{n} '.w.w = reshape(x(vector_index+(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+1:vector_index+(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1),1),[' blocks{n} '.J+1,' blocks{n} '.K+2]); ' blocks{n} '.w.w = ' blocks{n} '.w.w.'';'])
                eval([blocks{n} '.p.p = reshape(x(vector_index+(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+1:vector_index+(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+' blocks{n} '.K*' blocks{n} '.J,1),[' blocks{n} '.J,' blocks{n} '.K]); ' blocks{n} '.p.p = ' blocks{n} '.p.p.'';'])
                if flowopt.energy==1
                    eval([blocks{n} '.T.T = reshape(x(vector_index+(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+' blocks{n} '.K*' blocks{n} '.J+1:vector_index+(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+' blocks{n} '.K*' blocks{n} '.J+(' blocks{n} '.K+2)*(' blocks{n} '.J+2),1),[' blocks{n} '.J+2,' blocks{n} '.K+2]); ' blocks{n} '.T.T = ' blocks{n} '.T.T.'';'])
                    eval(['vector_index=vector_index+numel(' blocks{n} '.T.T);'])
                    
                end
                eval(['vector_index=vector_index+numel(' blocks{n} '.u.u)+numel(' blocks{n} '.w.w)+numel(' blocks{n} '.p.p)+numel(' blocks{n} '.Rs1)+numel(' blocks{n} '.Rsend);'])
                                
            else
                eval([blocks{n} '.u.u=interp2(old' blocks{n} '.ETA.u,old' blocks{n} '.XI.u,old' blocks{n} '.u.u,' blocks{n} '.ETA.u,' blocks{n} '.XI.u,''spline'');'])
                eval([blocks{n} '.w.w=interp2(old' blocks{n} '.ETA.w,old' blocks{n} '.XI.w,old' blocks{n} '.w.w,' blocks{n} '.ETA.w,' blocks{n} '.XI.w,''spline'');'])
                eval([blocks{n} '.p.p=interp2(old' blocks{n} '.ETA.p,old' blocks{n} '.XI.p,old' blocks{n} '.p.p,' blocks{n} '.ETA.p,' blocks{n} '.XI.p,''spline'');'])
                if flowopt.energy==1
                    eval([blocks{n} '.T.T=interp2(old' blocks{n} '.ETA.T,old' blocks{n} '.XI.T,old' blocks{n} '.T.T,' blocks{n} '.ETA.T,' blocks{n} '.XI.T,''spline'');'])
                    eval(['vector_index=vector_index+numel(old' blocks{n} '.T.T);'])
                    
                else
                    T=[];
                    
                end
                eval(['vector_index=vector_index+numel(old' blocks{n} '.u.u)+numel(old' blocks{n} '.w.w)+numel(old' blocks{n} '.p.p+numel(old' blocks{n} '.Rs1)+numel(old' blocks{n} '.Rsend));'])
                
            end
            
            eval(['u=' blocks{n} '.u.u.'';']); u=u(:);
            eval(['w=' blocks{n} '.w.w.'';']); w=w(:);
            eval(['p=' blocks{n} '.p.p.'';']); p=p(:);
            if flowopt.energy==1
                eval(['T=' blocks{n} '.T.T.'';']); T=T(:);

            else
                T=[];

            end

            eval(['x_new=[x_new;u;w;p;T;' blocks{n} '.Rs1;' blocks{n} '.Rsend];'])
                
            eval(['old' blocks{n} '.XI=' blocks{n} '.XI;'])
            eval(['old' blocks{n} '.ETA=' blocks{n} '.ETA;'])
            eval(['old' blocks{n} '.Z=' blocks{n} '.Z;'])
            eval(['old' blocks{n} '.geom=' blocks{n} '.geom;'])
        
            eval(['old' blocks{n} '.u.u = ' blocks{n} '.u.u;'])
            eval([blocks{n} '.u.w=interp2(' blocks{n} '.ETA.u,' blocks{n} '.XI.u,' blocks{n} '.u.u,' blocks{n} '.ETA.w,' blocks{n} '.XI.w);'])
            eval([blocks{n} '.u.p=interp2(' blocks{n} '.ETA.u,' blocks{n} '.XI.u,' blocks{n} '.u.u,' blocks{n} '.ETA.p,' blocks{n} '.XI.p);'])
            eval([blocks{n} '.u.v=interp2(' blocks{n} '.ETA.u,' blocks{n} '.XI.u,' blocks{n} '.u.u,' blocks{n} '.ETA.v,' blocks{n} '.XI.v);'])
            eval([blocks{n} '.u.T=interp2(' blocks{n} '.ETA.u,' blocks{n} '.XI.u,' blocks{n} '.u.u,' blocks{n} '.ETA.T,' blocks{n} '.XI.T);'])
            
            eval(['old' blocks{n} '.w.w = ' blocks{n} '.w.w;'])
            eval([blocks{n} '.w.u=interp2(' blocks{n} '.ETA.w,' blocks{n} '.XI.w,' blocks{n} '.w.w,' blocks{n} '.ETA.u,' blocks{n} '.XI.u);'])
            eval([blocks{n} '.w.p=interp2(' blocks{n} '.ETA.w,' blocks{n} '.XI.w,' blocks{n} '.w.w,' blocks{n} '.ETA.p,' blocks{n} '.XI.p);'])
            eval([blocks{n} '.w.v=interp2(' blocks{n} '.ETA.w,' blocks{n} '.XI.w,' blocks{n} '.w.w,' blocks{n} '.ETA.v,' blocks{n} '.XI.v);'])
            eval([blocks{n} '.w.T=interp2(' blocks{n} '.ETA.w,' blocks{n} '.XI.w,' blocks{n} '.w.w,' blocks{n} '.ETA.T,' blocks{n} '.XI.T);'])
            
            eval(['old' blocks{n} '.p.p = ' blocks{n} '.p.p;'])
            eval([blocks{n} '.p.v=interp1(' blocks{n} '.xi.w(2:end-1).'',' blocks{n} '.p.p,' blocks{n} '.xi.u.'',''linear'',''extrap'');'])
            eval([blocks{n} '.p.v=interp1(' blocks{n} '.eta.u(2:end-1).'',' blocks{n} '.p.v.'',' blocks{n} '.eta.w.'',''linear'',''extrap'').'';'])
                
            if flowopt.energy==1
                eval(['old' blocks{n} '.T.T = ' blocks{n} '.T.T;'])
                eval([blocks{n} '.T.u=interp2(' blocks{n} '.ETA.T,' blocks{n} '.XI.T,' blocks{n} '.T.T,' blocks{n} '.ETA.u,' blocks{n} '.XI.u,''linear'');'])
                eval([blocks{n} '.T.w=interp2(' blocks{n} '.ETA.T,' blocks{n} '.XI.T,' blocks{n} '.T.T,' blocks{n} '.ETA.w,' blocks{n} '.XI.w,''linear'');'])
                eval([blocks{n} '.T.p=interp2(' blocks{n} '.ETA.T,' blocks{n} '.XI.T,' blocks{n} '.T.T,' blocks{n} '.ETA.p,' blocks{n} '.XI.p,''linear'');'])
                eval([blocks{n} '.T.v=interp2(' blocks{n} '.ETA.T,' blocks{n} '.XI.T,' blocks{n} '.T.T,' blocks{n} '.ETA.v,' blocks{n} '.XI.v,''linear'');'])
                
            else
                T_init=str2func(eval([blocks{n} '.T_init.T']));
                eval([blocks{n} '.T.T=T_init(' blocks{n} '.Z.T,' blocks{n} '.R_cyl.T);']);
                eval([blocks{n} '.T.u=interp2(' blocks{n} '.ETA.T,' blocks{n} '.XI.T,' blocks{n} '.T.T,' blocks{n} '.ETA.u,' blocks{n} '.XI.u);']);
                eval([blocks{n} '.T.w=interp2(' blocks{n} '.ETA.T,' blocks{n} '.XI.T,' blocks{n} '.T.T,' blocks{n} '.ETA.w,' blocks{n} '.XI.w);']);
                eval([blocks{n} '.T.p=interp2(' blocks{n} '.ETA.T,' blocks{n} '.XI.T,' blocks{n} '.T.T,' blocks{n} '.ETA.p,' blocks{n} '.XI.p);']);
                eval([blocks{n} '.T.v=interp2(' blocks{n} '.ETA.T,' blocks{n} '.XI.T,' blocks{n} '.T.T,' blocks{n} '.ETA.v,' blocks{n} '.XI.v);']);    
                
            end
                                
        else
            eval(['[' blocks{n} ']=initial_condition(' blocks{n} ', flowopt);'])
            x0=eval(['[x0;' blocks{n} '.x]']);
            vector_index=numel(x0);
            
        end
    
    % boundary conditions
        eval(['[' blocks{n} ', Rs]=block_boundary_V3d2(' blocks{n} ', Rs, flowopt);'])
        
    % A
        eval(['[' blocks{n} ']=block_A(' blocks{n} ', flowopt);'])
        
    % J
        eval(['[' blocks{n} ']=block_J_V3d2(' blocks{n} ', flowopt);'])

    % coupling boundaries
        eval(['[' blocks{n} ']=block_coupling(' blocks{n} ', T_fluid, flowopt);'])
        
    if flowopt.stability==1
        % A_stability
            eval(['[' blocks{n} ']=block_A_stability_V3d2(' blocks{n} ', flowopt);'])

        % B_stability
            eval(['[' blocks{n} ']=block_B_stability(' blocks{n} ', flowopt);'])
    
    end
    
end

if isempty(x_new)~=1
    x=x_new;
    
end

% assemble matrices, vectors

clear A J A_stability B_stability b res.u res.w res.T res.p

if flowopt.stability~=1
    A=[];
    J=[];
    b=[];
    res.u=[];
    res.w=[];
    res.T=[];
    res.p=[];
    
else
    A_stability=[];
    B_stability=[];
    
end

for n=1:length(blocks)
    %clear Ar Jr Ar_stability Br_stability
    if flowopt.stability~=1
        Ar=[];
        Jr=[];
        
    else
        Ar_stability=[];
        Br_stability=[];
        
    end
    if flowopt.stability~=1
        for o=1:length(blocks)
            if strcmp(blocks{n},blocks{o})==0
                eval(['[' blocks{n} ']=block_coupling_counter_Jacobi_Rs(' blocks{n} ',' blocks{o} ', flowopt);']);

            else
                eval([blocks{n} '.Ja.Rs.coup_counter=sparse(size(' blocks{n} '.A,1),size(' blocks{n} '.A,2));'])

            end

        end
        
    end
    for o=1:length(blocks)
        if strcmp(blocks{n},blocks{o})==0
            eval(['[' blocks{n} ']=block_coupling_counter(' blocks{n} ',' blocks{o} ', flowopt);']);
            if flowopt.stability~=1
                Ar=[Ar eval([blocks{n} '.coup_counter'])];
                Jr=[Jr eval([blocks{n} '.Ja.coup_counter'])];
                
            elseif flowopt.stability==1
                Ar_stability=[Ar_stability eval([blocks{n} '.coup_counter'])];
                Br_stability=[Br_stability eval(['sparse(size(' blocks{n} '.coup_counter,1),size(' blocks{n} '.coup_counter,2))'])];
                
            end
                                    
        else
            if flowopt.stability~=1
                Ar=[Ar eval([blocks{n} '.A+' blocks{n} '.bc.A+' blocks{n} '.coup'])];
                Jr=[Jr eval([blocks{n} '.Jacobi+' blocks{n} '.bc.A+' blocks{n} '.Ja.coup+' blocks{n} '.Ja.Rs.coup_counter'])];
                
            elseif flowopt.stability==1
                Ar_stability=[Ar_stability eval([blocks{n} '.A_stability+' blocks{n} '.bc.A+' blocks{n} '.coup'])];
                Br_stability=[Br_stability eval([blocks{n} '.B_stability'])];
                
            end
            
        end
    end
    if flowopt.stability~=1
        A=[A; Ar];
        J=[J; Jr];
        b=[b; eval([blocks{n} '.rhs'])];
        if exist('x','var')==1
            eval(['res.u=[res.u; Ar(1:(' blocks{n} '.K+1)*(' blocks{n} '.J+2),:)*x - ' blocks{n} '.rhs(1:(' blocks{n} '.K+1)*(' blocks{n} '.J+2),1)];'])
            eval(['res.w=[res.w; Ar((' blocks{n} '.K+1)*(' blocks{n} '.J+2)+1:(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1),:)*x - ' blocks{n} '.rhs((' blocks{n} '.K+1)*(' blocks{n} '.J+2)+1:(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1),1)];'])
            if flowopt.energy==1
                eval(['res.p=[res.p; Ar((' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+(' blocks{n} '.K+2)*(' blocks{n} '.J+2)+1:(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+(' blocks{n} '.K+2)*(' blocks{n} '.J+2)+' blocks{n} '.K*' blocks{n} '.J,:)*x - ' blocks{n} '.rhs((' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+(' blocks{n} '.K+2)*(' blocks{n} '.J+2)+1:(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+(' blocks{n} '.K+2)*(' blocks{n} '.J+2)+' blocks{n} '.K*' blocks{n} '.J,1)];'])
                eval(['res.T=[res.T; Ar((' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+1:(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+(' blocks{n} '.K+2)*(' blocks{n} '.J+2),:)*x - ' blocks{n} '.rhs((' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+1:(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+(' blocks{n} '.K+2)*(' blocks{n} '.J+2),1)];'])

            else
                eval(['res.p=[res.p; Ar((' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+1:(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+' blocks{n} '.K*' blocks{n} '.J,:)*x - ' blocks{n} '.rhs((' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+1:(' blocks{n} '.K+1)*(' blocks{n} '.J+2)+(' blocks{n} '.K+2)*(' blocks{n} '.J+1)+' blocks{n} '.K*' blocks{n} '.J,1)];'])
            
            end
            
        end                
    
    else
        A_stability=[A_stability; Ar_stability];
        B_stability=[B_stability; Br_stability];
        
    end
        
end
clear Ar Jr Ar_stability Br_stability
if flowopt.stability==1
    for n=1:length(blocks)
        eval(['clear ' blocks{n}])

    end
    clear eigenvectors eigenvectors_cayley
    %{
    save(['intermediate_' file_name])
    clearvars -except file_name
    load(['intermediate_' file_name])
    %}
    
end
    
if flowopt.stability~=1
    if exist('x','var')==1
        %res_u=reshape(res.u,[b2.J+2,b2.K+1]).';
        %res_w=reshape(res.w,[b2.J+1,b2.K+2]).';
        %res_p=reshape(res.p,[b2.J,b2.K]).';
        res.u=norm(res.u);
        res.w=norm(res.w);
        res.p=norm(res.p);
        if flowopt.energy==1
            %res_T=reshape(res.T,[b2.J+2,b2.K+2]).';
            res.T=norm(res.T);
            if max(isnan([res.u, res.w, res.T, res.p]))
                error('residual = NaN')
            else
                residual=max([res.u, res.w, res.T, res.p]);
            end
            if exist('dx','var')
                fprintf('%2d u-residual = %3.2e; w-residual = %3.2e; T-residual = %3.2e; conti-residual = %3.2e; delta_x = %3.2e\n',[iteration, res.u, res.w, res.T, res.p, norm(dx)])
            else
                fprintf('%2d u-residual = %3.2e; w-residual = %3.2e; T-residual = %3.2e; conti-residual = %3.2e\n',[iteration, res.u, res.w, res.T, res.p])
            end

        else
            if max(isnan([res.u, res.w, res.p]))
                error('residual = NaN')
            else
                residual=max([res.u, res.w, res.p]);
            end
            if exist('dx','var')
                fprintf('%2d u-residual = %3.2e; w-residual = %3.2e; conti-residual = %3.2e; delta_x = %3.2e\n',[iteration, res.u, res.w, res.p, norm(dx)])
            else
                fprintf('%2d u-residual = %3.2e; w-residual = %3.2e; conti-residual = %3.2e\n',[iteration, res.u, res.w, res.p])
            end

        end

    else
        x=x0;

    end
    %%{
    scaling=max(max(abs(A*x-b),[],2), max(abs(J),[],2));
    A=spdiags(1./scaling,0,length(A(1,:)),length(A(:,1)))*A;
    J=spdiags(1./scaling,0,length(A(1,:)),length(A(:,1)))*J;
    b=b./scaling;
    clear scaling
    %}
    
else
    %%{
    scaling=max(max(abs(A_stability),[],2), max(abs(B_stability),[],2));
    A_stability=spdiags(1./scaling,0,length(A_stability(1,:)),length(A_stability(:,1)))*A_stability;
    B_stability=spdiags(1./scaling,0,length(A_stability(1,:)),length(A_stability(:,1)))*B_stability;
    clear scaling
    %}
    
end