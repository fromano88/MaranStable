function[b,Rs]=block_boundary(b, Rs, flowopt)
    ax=flowopt.ax;
    m=flowopt.m;
    g=flowopt.g;
    energy=flowopt.energy;
    stability=flowopt.stability;
    iv=flowopt.iv;
    rho=str2func(b.rho);
    
%%%%%%%%%%%%%%%%%%%%%%% lower/upper boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bc.rhs.C=sparse(b.K*b.J,1);
    if max(strncmp('ss', b.bc.z(:),2))==1 && (max(strncmp('o', b.bc.z(:),1)) || max(strncmp('o', b.bc.r(:),1)))~=1 && isfield(b, 'V_l')==1
        ri=find(abs(b.r.w(2:end-1)-b.rp)==min(min(abs(b.r.w(2:end-1)-b.rp))));
        zi=find(abs(b.z.u(2:end-1)-b.zp)==min(min(abs(b.z.u(2:end-1)-b.zp))));
        bc.rhs.C((ri(1)-1)*b.J+zi(1))=b.V_l;
        
    elseif (max(strncmp('o', b.bc.z(:),1)) || max(strncmp('o', b.bc.r(:),1)))~=1
        ri=find(abs(b.r.w(2:end-1)-b.rp)==min(min(abs(b.r.w(2:end-1)-b.rp))));
        zi=find(abs(b.z.u(2:end-1)-b.zp)==min(min(abs(b.z.u(2:end-1)-b.zp))));
        bc.rhs.C((ri(1)-1)*b.J+zi(1))=b.p_fluid;
    end    
    
    bc.rhs.u1=[];
    bc.rhs.w1=[];
    bc.rhs.T1=[];
    bc.rm.u1=[];
    bc.rm.w1=[];
    bc.rm.p1=[];
    bc.st.phim.v1=[];
    bc.zm.u1=[];
    bc.zm.w1=[];
    bc.th_e.T1=[];
    b.C.R1=[];
    b.counter.C.R1=[];
    b.Ja.C.R1=[];
    b.Ja.counter.C.R1=[];
    b.dbc.R1=[];
    
    bc.rhs.uend=[];
    bc.rhs.wend=[];
    bc.rhs.Tend=[];
    bc.rm.uend=[];
    bc.rm.wend=[];
    bc.rm.pend=[];
    bc.st.phim.vend=[];
    bc.zm.uend=[];
    bc.zm.wend=[];
    bc.th_e.Tend=[];
    b.C.Rend=[];
    b.counter.C.Rend=[];
    b.Ja.C.Rend=[];
    b.Ja.counter.C.Rend=[];
    b.dbc.Rend=[];
    
%%%%%%%%%%%%%%%%%%%%%%% left / right boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % r-momentum
        rm_u_Ssw=[];
        rm_u_Sse=[];
        rm_u_w=[];
        rm_u_Cl=[];
        rm_u_Cr=[];
        rm_u_e=[];
        rm_u_Nnw=[];
        rm_u_Nne=[];
        rm_u_W=[];
        rm_u_E=[];

        rm_w_swl=[];
        rm_w_sel=[];
        rm_w_nwl=[];
        rm_w_nel=[];
        rm_w_swr=[];
        rm_w_ser=[];
        rm_w_nwr=[];
        rm_w_ner=[];
        
        bc.rhs.ul=[];
        bc.rhs.ur=[];
        
    % phi-momentum
        phim_v_Ssw=[];
        phim_v_Sse=[];
        phim_v_w=[];
        phim_v_Cl=[];
        phim_v_Cr=[];
        phim_v_e=[];
        phim_v_Nnw=[];
        phim_v_Nne=[];
        phim_v_W=[];
        phim_v_E=[];
        
    % z-momentum
        zm_w_Ssw=[];
        zm_w_Sse=[];
        zm_w_w=[];
        zm_w_Cl=[];
        zm_w_Cr=[];
        zm_w_e=[];
        zm_w_Nnw=[];
        zm_w_Nne=[];
        zm_w_W=[];
        zm_w_E=[];
        zm_p_Cl=[];
        zm_p_Cr=[];
        
        bc.rhs.wl=[];
        bc.rhs.wr=[];

    % thermal energy
        th_e_T_Ssw=[];
        th_e_T_Sse=[];
        th_e_T_w=[];
        th_e_T_Cl=[];
        th_e_T_Cr=[];
        th_e_T_e=[];
        th_e_T_Nnw=[];
        th_e_T_Nne=[];
        th_e_T_W=[];
        th_e_T_E=[];
        
        bc.rhs.Tl=[];
        bc.rhs.Tr=[];

    N1.u=sqrt(1+b.R1_z.u(1,:).^2);
    N1.w=sqrt(1+b.R1_z.w(1,:).^2);
    Nend.u=sqrt(1+b.Rend_z.u(end,:).^2);
    Nend.w=sqrt(1+b.Rend_z.w(end,:).^2);
    
    b.sizeR1=0;
    b.sizeRend=0;
        
    for o=1:2
        for n=1:length(b.bc.z)/2
            J = length(b.Z.u(b.Z.u(1,:)>b.geom.z(n) & b.Z.u(1,:)<b.geom.z(n+1)));

            % r-momentum
                boundary_r_momentum_z

            if stability==1
                % phi-momentum
                    boundary_phi_momentum_z_stability

            end

            % z-momentum
                boundary_z_momentum_z

            % thermal energy
                if energy==1
                    boundary_thermal_energy_z % outer most boundaries in z direction + rhs for newton's cooling law

                end

        end
       
        for n=1:length(b.bc.r)/2
            K = length(b.R_cyl.w(b.R_cyl.w(:,1)>b.geom.r(n) & b.R_cyl.w(:,1)<b.geom.r(n+1)));

            % r-momentum
                boundary_r_momentum_r

            if stability==1
                % phi-momentum
                    boundary_phi_momentum_r_stability

            end

            % z-momentum
                boundary_z_momentum_r

            % thermal energy
                if energy==1
                    boundary_thermal_energy_r % outer most boundaries in r direction

                end

        end
        
    end

    % gravity
        if isfield(flowopt,'boussinesq')==1 && flowopt.boussinesq==1
            rho = @(theta) rho(b.T_0).*(1-b.beta*(theta-b.T_0));
        end
        
        gravity.rhs.w=[sparse(b.K,1) g*rho(b.T.w(2:end-1,2:end-1)).*b.R.w(2:end-1,2:end-1).^ax.*b.JA.w(2:end-1,2:end-1).*b.DXI.zm(:,2:end-1).*b.DETA.zm(:,2:end-1) sparse(b.K,1)];
        gravity.rhs.w=gravity.rhs.w.'; gravity.rhs.w=gravity.rhs.w(:);
        
        rho = str2func(b.rho);
    
%%%%%%%%%%%%%%%%%%%%%%% assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % r-momentum
        rm_u_s=rm_u_Ssw+rm_u_Sse;
        rm_u_n=rm_u_Nnw+rm_u_Nne;
        rm_u_C=rm_u_w+rm_u_e;

        rm_u_s=rm_u_s.'; rm_u_s=rm_u_s(:);
        rm_u_Cl=rm_u_Cl.'; rm_u_Cl=rm_u_Cl(:);
        rm_u_C=rm_u_C.'; rm_u_C=rm_u_C(:);
        rm_u_Cr=rm_u_Cr.'; rm_u_Cr=rm_u_Cr(:);
        rm_u_n=rm_u_n.'; rm_u_n=rm_u_n(:);
        rm_u_W=rm_u_W.'; rm_u_W=rm_u_W(:);
        rm_u_E=rm_u_E.'; rm_u_E=rm_u_E(:);
        bc.rm.u = [bc.rm.u1; spdiags([rm_u_s rm_u_W rm_u_Cr rm_u_C rm_u_Cl rm_u_E rm_u_n],[0 b.J b.J+1 b.J+2 b.J+3 b.J+4 2*(b.J+2)],(b.K-1)*(b.J+2),(b.K+1)*(b.J+2)); bc.rm.uend];
        
        rm_w_sw=rm_w_swl + rm_w_swr;
        rm_w_se=rm_w_sel + rm_w_ser;
        rm_w_nw=rm_w_nwl + rm_w_nwr;
        rm_w_ne=rm_w_nel + rm_w_ner;
        
        rm_w_sw=rm_w_sw.'; rm_w_sw=rm_w_sw(:);
        rm_w_se=rm_w_se.'; rm_w_se=rm_w_se(:);
        rm_w_nw=rm_w_nw.'; rm_w_nw=rm_w_nw(:);
        rm_w_ne=rm_w_ne.'; rm_w_ne=rm_w_ne(:);
        
        bc.rm.w = spdiags([rm_w_sw rm_w_se rm_w_nw rm_w_ne],[b.J+3 b.J+4 2*(b.J+2)+2 2*(b.J+2)+3],(b.K-1)*(b.J+3),(b.K+2)*(b.J+3));
        bc.rm.w(b.J+3:b.J+3:(b.K-1)*(b.J+3),:) = [];
        bc.rm.w(:,4:b.J+3:(b.K+2)*(b.J+3)) = [];
        bc.rm.w(:,4:b.J+2:(b.K+2)*(b.J+2)) = [];
        bc.rm.w = [bc.rm.w1; bc.rm.w; bc.rm.wend];
        
    % phi-momentum
        if stability==1
            phim_v_s=phim_v_Ssw+phim_v_Sse;
            phim_v_n=phim_v_Nnw+phim_v_Nne;
            phim_v_C=phim_v_w+phim_v_e;

            phim_v_s=phim_v_s.'; phim_v_s=phim_v_s(:);
            phim_v_Cl=phim_v_Cl.'; phim_v_Cl=phim_v_Cl(:);
            phim_v_C=phim_v_C.'; phim_v_C=phim_v_C(:);
            phim_v_Cr=phim_v_Cr.'; phim_v_Cr=phim_v_Cr(:);
            phim_v_n=phim_v_n.'; phim_v_n=phim_v_n(:);
            phim_v_W=phim_v_W.'; phim_v_W=phim_v_W(:);
            phim_v_E=phim_v_E.'; phim_v_E=phim_v_E(:);
            bc.st.phim.v = [bc.st.phim.v1; spdiags(iv*[phim_v_s phim_v_W phim_v_Cr phim_v_C phim_v_Cl phim_v_E phim_v_n],[0 b.J b.J+1 b.J+2 b.J+3 b.J+4 2*(b.J+2)],b.K*(b.J+2),(b.K+2)*(b.J+2)); bc.st.phim.vend];
            
        end

    % z-momentum
        zm_w_s=zm_w_Ssw+zm_w_Sse;
        zm_w_n=zm_w_Nnw+zm_w_Nne;
        zm_w_C=zm_w_w+zm_w_e;
        zm_p_C=zm_p_Cl+zm_p_Cr;
        
        zm_w_s=zm_w_s.'; zm_w_s=zm_w_s(:);
        zm_w_Cl=zm_w_Cl.'; zm_w_Cl=zm_w_Cl(:);
        zm_w_C=zm_w_C.'; zm_w_C=zm_w_C(:);
        zm_w_Cr=zm_w_Cr.'; zm_w_Cr=zm_w_Cr(:);
        zm_w_n=zm_w_n.'; zm_w_n=zm_w_n(:);
        zm_w_W=zm_w_W.'; zm_w_W=zm_w_W(:);
        zm_w_E=zm_w_E.'; zm_w_E=zm_w_E(:);
        zm_p_C=zm_p_C.'; zm_p_C=zm_p_C(:);
        bc.zm.u = [bc.zm.u1; sparse(b.K*(b.J+1),(b.K+1)*(b.J+2)); bc.zm.uend];
        bc.zm.w = [bc.zm.w1; spdiags([zm_w_s zm_w_W zm_w_Cr zm_w_C zm_w_Cl zm_w_E zm_w_n],[0 b.J-1 b.J b.J+1 b.J+2 b.J+3 2*(b.J+1)],b.K*(b.J+1),(b.K+2)*(b.J+1)); bc.zm.wend];
        bc.zm.p = [sparse(b.J+1,b.K*(b.J+1)); spdiags(zm_p_C,0,b.K*(b.J+1),b.K*(b.J+1)); sparse(b.J+1,b.K*(b.J+1))];
        bc.zm.p(:,b.J:b.J+1:end)=[];
        
    % thermal energy
        if energy==1
            th_e_T_s=th_e_T_Ssw+th_e_T_Sse;
            th_e_T_n=th_e_T_Nnw+th_e_T_Nne;
            th_e_T_C=th_e_T_w+th_e_T_e;

            th_e_T_s=th_e_T_s.'; th_e_T_s=th_e_T_s(:);
            th_e_T_Cl=th_e_T_Cl.'; th_e_T_Cl=th_e_T_Cl(:);
            th_e_T_C=th_e_T_C.'; th_e_T_C=th_e_T_C(:);
            th_e_T_Cr=th_e_T_Cr.'; th_e_T_Cr=th_e_T_Cr(:);
            th_e_T_n=th_e_T_n.'; th_e_T_n=th_e_T_n(:);
            th_e_T_W=th_e_T_W.'; th_e_T_W=th_e_T_W(:);
            th_e_T_E=th_e_T_E.'; th_e_T_E=th_e_T_E(:);
            bc.th_e.T = [bc.th_e.T1; spdiags([th_e_T_s th_e_T_W th_e_T_Cr th_e_T_C th_e_T_Cl th_e_T_E th_e_T_n],[0 b.J b.J+1 b.J+2 b.J+3 b.J+4 2*(b.J+2)],b.K*(b.J+2),(b.K+2)*(b.J+2)); bc.th_e.Tend];

        end
            
    if stability~=1
        if energy==1
            b.bc.A=[bc.rm.u  bc.rm.w  [bc.rm.p1; sparse((b.K-1)*(b.J+2),b.K*b.J); bc.rm.pend]  sparse((b.K+1)*(b.J+2),(b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend);
                    bc.zm.u  bc.zm.w  bc.zm.p                                                  sparse((b.K+2)*(b.J+1),(b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend);
                    sparse((b.K+2)*(b.J+2),(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1) + b.K*b.J)          bc.th_e.T  sparse((b.K+2)*(b.J+2),b.sizeR1 + b.sizeRend);
                    sparse(b.K*b.J, (b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend)
                    sparse(b.sizeR1 + b.sizeRend, (b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend)];

            b.rhs=[[bc.rhs.u1; bc.rhs.ul+bc.rhs.ur; bc.rhs.uend]; [bc.rhs.w1; bc.rhs.wl+bc.rhs.wr+gravity.rhs.w; bc.rhs.wend]; [bc.rhs.T1; bc.rhs.Tl+bc.rhs.Tr; bc.rhs.Tend]; bc.rhs.C; sparse(b.sizeR1 + b.sizeRend,1)];

        else
            b.bc.A=[bc.rm.u  bc.rm.w  [bc.rm.p1; sparse((b.K-1)*(b.J+2),b.K*b.J); bc.rm.pend]  sparse((b.K+1)*(b.J+2),b.sizeR1 + b.sizeRend);
                    bc.zm.u  bc.zm.w  bc.zm.p                                                  sparse((b.K+2)*(b.J+1),b.sizeR1 + b.sizeRend);
                    sparse(b.K*b.J, (b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + + b.sizeR1 + b.sizeRend)
                    sparse(b.sizeR1 + b.sizeRend, (b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + b.sizeRend)];
    
            b.rhs=[[bc.rhs.u1; bc.rhs.ul+bc.rhs.ur; bc.rhs.uend]; [bc.rhs.w1; bc.rhs.wl+bc.rhs.wr+gravity.rhs.w; bc.rhs.wend]; bc.rhs.C; sparse(b.sizeR1 + b.sizeRend,1)];

        end
        
    else
        if energy==1
            b.bc.A=[bc.rm.u                                  sparse((b.K+1)*(b.J+2),(b.K+2)*(b.J+2))  bc.rm.w  [bc.rm.p1; sparse((b.K-1)*(b.J+2),b.K*b.J); bc.rm.pend]  sparse((b.K+1)*(b.J+2),(b.K+2)*(b.J+2));
                    sparse((b.K+2)*(b.J+2),(b.K+1)*(b.J+2))  bc.st.phim.v                             sparse((b.K+2)*(b.J+2),(b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2));
                    bc.zm.u                                  sparse((b.K+2)*(b.J+1),(b.K+2)*(b.J+2))  bc.zm.w  bc.zm.p                                                  sparse((b.K+2)*(b.J+1),(b.K+2)*(b.J+2));
                    sparse((b.K+2)*(b.J+2),(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1) + b.K*b.J)                                                                   bc.th_e.T;
                    sparse(b.K*b.J, (b.K+1)*(b.J+2) + (b.K+2)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2))];
    
        else
            b.bc.A=[bc.rm.u                                  sparse((b.K+1)*(b.J+2),(b.K+2)*(b.J+2))  bc.rm.w  [bc.rm.p1; sparse((b.K-1)*(b.J+2),b.K*b.J); bc.rm.pend];
                    sparse((b.K+2)*(b.J+2),(b.K+1)*(b.J+2))  bc.st.phim.v                             sparse((b.K+2)*(b.J+2),(b.K+2)*(b.J+1) + b.K*b.J);
                    bc.zm.u                                  sparse((b.K+2)*(b.J+1),(b.K+2)*(b.J+2))  bc.zm.w  bc.zm.p;
                    sparse(b.K*b.J, (b.K+1)*(b.J+2) + (b.K+2)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J)];
    
        end
                
    end