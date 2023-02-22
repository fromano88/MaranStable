function[b]=block_coupling_counter_Jacobi_Rs(b, c, flowopt)
    energy=flowopt.energy;
    coupl_counter.rm.Rs1=sparse(b.J+2, b.sizeR1);
    coupl_counter.rm.Rsend=sparse(b.J+2, b.sizeRend);
    coupl_counter.zm.Rs1=sparse(b.J+1, b.sizeR1);
    coupl_counter.zm.Rsend=sparse(b.J+1, b.sizeRend);
    coupl_counter.the.Rs1=sparse(b.J+2, b.sizeR1);
    coupl_counter.the.Rsend=sparse(b.J+2, b.sizeRend);
    
    for n=1:length(b.bc.z)/2
        J=length(b.z.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)));
        
%%%%%%%%%%%%%%%%%%%%%%% lower boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % surface, interior or periodic
        if (strncmp('ss', {b.bc.z{n}},2) || strcmp(b.bc.z{n}(1),'c')==1 || strcmp(b.bc.z{n}(1),'p')==1) && eval(['b.' b.bc.z{n} '.o==1']) && max(strcmp(b.bc.z{n},c.bc.z))==1
            coupl_counter.rm.Rs1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) = -eval(['c.' b.bc.z{n} '.Ja.rmend.Rs']);
            coupl_counter.zm.Rs1(b.z.w>=b.geom.z(n) & b.z.w<=b.geom.z(n+1),:) = -eval(['c.' b.bc.z{n} '.Ja.zmend.Rs']);
            if energy==1
                coupl_counter.the.Rs1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) = -eval(['c.' b.bc.z{n} '.Ja.th_eend.Rs']);
                   
            end
                        
        end

%%%%%%%%%%%%%%%%%%%%%%% upper boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % surface, interior or periodic
        if (strncmp('ss', {b.bc.z{length(b.bc.z)/2+n}},2) || strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'c')==1 || strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'p')==1) && eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.o==1']) && max(strcmp(b.bc.z{length(b.bc.z)/2+n},c.bc.z))==1
            coupl_counter.rm.Rsend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) = -eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.Ja.rm1.Rs']);
            coupl_counter.zm.Rsend(b.z.w>=b.geom.z(n) & b.z.w<=b.geom.z(n+1),:) = -eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.Ja.zm1.Rs']);
            if energy==1
                coupl_counter.the.Rsend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) = -eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.Ja.th_e1.Rs']);
                        
            end
            
        end
        
    end
    
    if energy==1
        b.Ja.Rs.coup_counter=sparse((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1 + b.sizeRend,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1 + b.sizeRend);
        b.Ja.Rs.coup_counter(1:(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1) = coupl_counter.rm.Rs1;
        b.Ja.Rs.coup_counter(b.K*(b.J+2)+1:(b.K+1)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1 + b.sizeRend) = coupl_counter.rm.Rsend;
        b.Ja.Rs.coup_counter((b.K+1)*(b.J+2)+1:(b.K+1)*(b.J+2)+(b.J+1),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1) = coupl_counter.zm.Rs1;
        b.Ja.Rs.coup_counter((b.K+1)*(b.J+2)+(b.K+1)*(b.J+1)+1:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1 + b.sizeRend) = coupl_counter.zm.Rsend;b.Ja.Rs.coup_counter((b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+1:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1) = coupl_counter.the.Rs1;
        b.Ja.Rs.coup_counter((b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+(b.K+1)*(b.J+2)+1:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+(b.K+2)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1 + b.sizeRend) = coupl_counter.the.Rsend;
    
    else
        b.Ja.Rs.coup_counter=sparse((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + b.sizeRend,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + b.sizeRend);
        b.Ja.Rs.coup_counter(1:(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1) = coupl_counter.rm.Rs1;
        b.Ja.Rs.coup_counter(b.K*(b.J+2)+1:(b.K+1)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + b.sizeRend) = coupl_counter.rm.Rsend;
        b.Ja.Rs.coup_counter((b.K+1)*(b.J+2)+1:(b.K+1)*(b.J+2)+(b.J+1),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1) = coupl_counter.zm.Rs1;
        b.Ja.Rs.coup_counter((b.K+1)*(b.J+2)+(b.K+1)*(b.J+1)+1:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + b.sizeRend) = coupl_counter.zm.Rsend;
        
    end