function[b]=block_J_V3d2(b, flowopt)
    energy=flowopt.energy;
    
    [b]=r_momentum_Jacobian(b, flowopt);
    [b]=z_momentum_Jacobian(b, flowopt);
    if energy==1
        [b]=thermal_energy_Jacobian_V3d2(b, flowopt);
        
    end
    [b]=continuity_Jacobian(b, flowopt);

    if energy==1
        b.Jacobi=[b.Ja.r_m.u   b.Ja.r_m.w   b.r_m.p                          b.Ja.r_m.T   sparse((b.K+1)*(b.J+2),b.sizeR1 + b.sizeRend);
                  b.Ja.z_m.u   b.Ja.z_m.w   b.z_m.p                          b.Ja.z_m.T   sparse((b.K+2)*(b.J+1),b.sizeR1 + b.sizeRend);
                  b.Ja.th_e.u  b.Ja.th_e.w  sparse((b.K+2)*(b.J+2),b.K*b.J)  b.Ja.th_e.T  sparse((b.K+2)*(b.J+2),b.sizeR1+b.sizeRend);
                  b.C.u        b.C.w        b.C.p                            b.Ja.C.T     sparse(b.K*b.J,b.sizeR1 + b.sizeRend);
                  sparse(b.sizeR1 + b.sizeRend, (b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) +b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend)];
              
    else
        b.Jacobi=[b.Ja.r_m.u   b.Ja.r_m.w   b.r_m.p                          sparse((b.K+1)*(b.J+2),b.sizeR1 + b.sizeRend);
                  b.Ja.z_m.u   b.Ja.z_m.w   b.z_m.p                          sparse((b.K+2)*(b.J+1),b.sizeR1 + b.sizeRend);
                  b.C.u        b.C.w        b.C.p                            sparse(b.K*b.J,b.sizeR1 + b.sizeRend);
                  sparse(b.sizeR1 + b.sizeRend, (b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) +b.K*b.J + b.sizeR1 + b.sizeRend)];
        
    end