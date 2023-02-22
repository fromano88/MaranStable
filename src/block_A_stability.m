function[b]=block_A_stability(b, flowopt)
    energy=flowopt.energy;
    
    [b]=r_momentum_stability(b, flowopt);
    [b]=phi_momentum_stability(b, flowopt);
    [b]=z_momentum_stability(b, flowopt);
    if energy==1
        [b]=thermal_energy_stability(b, flowopt);
        
    end
    [b]=continuity_stability(b, flowopt);
    
    if energy==1
        b.A_stability=[b.st.r_m.u    b.st.r_m.v    b.st.r_m.w    b.r_m.p                          b.Ja.r_m.T;
                       b.st.phi_m.u  b.st.phi_m.v  b.st.phi_m.w  b.st.phi_m.p                     b.st.phi_m.T;
                       b.Ja.z_m.u    b.st.z_m.v    b.st.z_m.w    b.z_m.p                          b.Ja.z_m.T;
                       b.Ja.th_e.u   b.st.th_e.v   b.Ja.th_e.w   sparse((b.K+2)*(b.J+2),b.K*b.J)  b.st.th_e.T;
                       b.C.u         b.st.C.v      b.C.w         b.C.p                            b.Ja.C.T     ];
                   
    else
        b.A_stability=[b.st.r_m.u    b.st.r_m.v    b.st.r_m.w    b.r_m.p;
                       b.st.phi_m.u  b.st.phi_m.v  b.st.phi_m.w  b.st.phi_m.p;
                       b.Ja.z_m.u    b.st.z_m.v    b.st.z_m.w    b.z_m.p;
                       b.C.u         b.st.C.v      b.C.w         b.C.p];
        
    end