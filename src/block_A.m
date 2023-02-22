function[b]=block_A(b, flowopt)
    energy=flowopt.energy;
    
    [b]=r_momentum(b, flowopt);
    [b]=z_momentum(b, flowopt);
    if energy==1
        [b]=thermal_energy(b, flowopt);
        
    end
    [b]=continuity(b, flowopt);
    
    if energy==1
        b.A=[b.r_m.u                                   b.r_m.w          b.r_m.p                   sparse((b.K+1)*(b.J+2),(b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend);
             b.z_m.u                                   b.z_m.w          b.z_m.p                   sparse((b.K+2)*(b.J+1),(b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend);
             sparse((b.K+2)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J)                  b.th_e.T                                 sparse((b.K+2)*(b.J+2),b.sizeR1 + b.sizeRend);
             b.C.u                                     b.C.w            b.C.p                     sparse(b.K*b.J,(b.K+2)*(b.J+2)         + b.sizeR1 + b.sizeRend);
             sparse(b.sizeR1 + b.sizeRend, (b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) +b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend)];

    else
        b.A=[b.r_m.u                                   b.r_m.w          b.r_m.p                   sparse((b.K+1)*(b.J+2),b.sizeR1 + b.sizeRend);
             b.z_m.u                                   b.z_m.w          b.z_m.p                   sparse((b.K+2)*(b.J+1),b.sizeR1 + b.sizeRend);
             b.C.u                                     b.C.w            b.C.p                     sparse(b.K*b.J,        b.sizeR1 + b.sizeRend);
             sparse(b.sizeR1 + b.sizeRend, (b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) +b.K*b.J + b.sizeR1 + b.sizeRend)];

    end