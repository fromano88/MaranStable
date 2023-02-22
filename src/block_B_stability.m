function[b]=block_B_stability(b, flowopt)
    energy=flowopt.energy;
    
    if energy==1
        b.B_stability=[b.st_t.r_m.u                             sparse((b.K+1)*(b.J+2),(b.K+2)*(b.J+2))  b.st_t.r_m.w        sparse((b.K+1)*(b.J+2),b.K*b.J)  b.st_t.r_m.T;
                       sparse((b.K+2)*(b.J+2),(b.K+1)*(b.J+2))  b.st_t.phi_m.v                           sparse((b.K+2)*(b.J+2),(b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2));
                       sparse((b.K+2)*(b.J+1),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2))                         b.st_t.z_m.w        sparse((b.K+2)*(b.J+1),b.K*b.J)  b.st_t.z_m.T;
                       sparse((b.K+2)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J)                                                  b.st_t.th_e.T;
                       sparse(b.K*b.J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J)                                                          b.st_t.C.T    ];
                   
    else
        b.B_stability=[b.st_t.r_m.u                             sparse((b.K+1)*(b.J+2),(b.K+2)*(b.J+2))  b.st_t.r_m.w        sparse((b.K+1)*(b.J+2),b.K*b.J);
                       sparse((b.K+2)*(b.J+2),(b.K+1)*(b.J+2))  b.st_t.phi_m.v                           sparse((b.K+2)*(b.J+2),(b.K+2)*(b.J+1) + b.K*b.J);
                       sparse((b.K+2)*(b.J+1),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2))                         b.st_t.z_m.w        sparse((b.K+2)*(b.J+1),b.K*b.J);
                       sparse(b.K*b.J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J)                                                        ];
                          
    end