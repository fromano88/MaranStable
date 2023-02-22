function real_pert = plus_cc(p,m,phi,gamma,time)

    % returns the sum with the complex conjugate
    
    if nargin == 3
        gamma = 0;
        time  = 0;
    end

    real_pert = p.*exp(-gamma(1,1)*time+1i*m*phi)+conj(p.*exp(-gamma(1,1)*time+1i*m*phi));