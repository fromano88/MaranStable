function[b]=block_coupling_counter(b, c, flowopt)
    energy=flowopt.energy;
    stability=flowopt.stability;
    
    if stability~=1
        if energy==1
            b.coupl_counter.rm1   = sparse(b.J+2,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1 + c.sizeRend);
            b.coupl_counter.zm1   = sparse(b.J+1,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1 + c.sizeRend);
            b.coupl_counter.th_e1 = sparse(b.J+2,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1 + c.sizeRend);

            b.Ja.coupl_counter.rm1   = b.coupl_counter.rm1;
            b.Ja.coupl_counter.zm1   = b.coupl_counter.zm1;
            b.Ja.coupl_counter.th_e1 = b.coupl_counter.th_e1;

            b.rhs_counter.rm1 = sparse(b.J+2,1);
            b.rhs_counter.zm1 = sparse(b.J+1,1);

            b.coupl_counter.rmend   = sparse(b.J+2,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1 + c.sizeRend);
            b.coupl_counter.zmend   = sparse(b.J+1,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1 + c.sizeRend);
            b.coupl_counter.th_eend = sparse(b.J+2,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1 + c.sizeRend);

            b.Ja.coupl_counter.rmend   = b.coupl_counter.rmend;
            b.Ja.coupl_counter.zmend   = b.coupl_counter.zmend;
            b.Ja.coupl_counter.th_eend = b.coupl_counter.th_eend;

            b.rhs_counter.rmend = sparse(b.J+2,1);
            b.rhs_counter.zmend = sparse(b.J+1,1);
            
        else
            b.coupl_counter.rm1   = sparse(b.J+2,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + c.sizeR1 + c.sizeRend);
            b.coupl_counter.zm1   = sparse(b.J+1,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + c.sizeR1 + c.sizeRend);
            
            b.Ja.coupl_counter.rm1   = b.coupl_counter.rm1;
            b.Ja.coupl_counter.zm1   = b.coupl_counter.zm1;
            
            b.rhs_counter.rm1 = sparse(b.J+2,1);
            b.rhs_counter.zm1 = sparse(b.J+1,1);

            b.coupl_counter.rmend   = sparse(b.J+2,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + c.sizeR1 + c.sizeRend);
            b.coupl_counter.zmend   = sparse(b.J+1,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + c.sizeR1 + c.sizeRend);
            
            b.Ja.coupl_counter.rmend   = b.coupl_counter.rmend;
            b.Ja.coupl_counter.zmend   = b.coupl_counter.zmend;
            
            b.rhs_counter.rmend = sparse(b.J+2,1);
            b.rhs_counter.zmend = sparse(b.J+1,1);
            
        end
        
    else
        if energy==1
            b.st.coupl_counter.phim1 = sparse(b.J+2,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2));
            b.st.coupl_counter.zm1   = sparse(b.J+1,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2));
            b.st.coupl_counter.th_e1 = sparse(b.J+2,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2));

            b.st.coupl_counter.phimend = sparse(b.J+2,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2));
            b.st.coupl_counter.zmend   = sparse(b.J+1,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2));
            b.st.coupl_counter.th_eend = sparse(b.J+2,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2));
            
        else
            b.st.coupl_counter.phim1 = sparse(b.J+2,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J);
            b.st.coupl_counter.zm1   = sparse(b.J+1,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J);
            
            b.st.coupl_counter.phimend = sparse(b.J+2,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J);
            b.st.coupl_counter.zmend   = sparse(b.J+1,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J);
                        
        end
        
    end
    
    for n=1:length(b.bc.z)/2
        J=length(b.z.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)));
        
%%%%%%%%%%%%%%%%%%%%%%% lower boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % surface, interior or periodic
        if (strcmp(b.bc.z{n}(1),'s')==1 || strcmp(b.bc.z{n}(1),'c')==1 || strcmp(b.bc.z{n}(1),'p')==1) && max(strcmp(b.bc.z{n},c.bc.z))==1
            if stability~=1
                if strcmp(b.bc.z{n}(2),'i')~=1 && strcmp(b.bc.z{n}(2),'r')~=1
                    b.coupl_counter.rm1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)      = [eval(['c.' b.bc.z{n} '.rm']) sparse(J,c.sizeR1 + c.sizeRend)];
                    b.Ja.coupl_counter.rm1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)   = [eval(['c.' b.bc.z{n} '.Ja.rm']) sparse(J,c.sizeR1 + c.sizeRend)];
                    b.rhs_counter.rm1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),1)=eval(['c.' b.bc.z{n} '.rhs.rm']);
                    
                end
                b.coupl_counter.zm1(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)      = [eval(['c.' b.bc.z{n} '.zm']) sparse(J-1,c.sizeR1 + c.sizeRend)];
                b.Ja.coupl_counter.zm1(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)   = [eval(['c.' b.bc.z{n} '.Ja.zm']) sparse(J-1,c.sizeR1 + c.sizeRend)];
                b.rhs_counter.zm1(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),1)=eval(['c.' b.bc.z{n} '.rhs.zm']);
                
                if energy==1
                    b.coupl_counter.th_e1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)    = [eval(['c.' b.bc.z{n} '.th_e']) sparse(J,c.sizeR1 + c.sizeRend)];
                    b.Ja.coupl_counter.th_e1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) = [eval(['c.' b.bc.z{n} '.Ja.th_e']) sparse(J,c.sizeR1 + c.sizeRend)];

                end
                
            else
                b.st.coupl_counter.phim1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)    = eval(['c.' b.bc.z{n} '.st.phim']);
                b.st.coupl_counter.zm1(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)      = eval(['c.' b.bc.z{n} '.st.zm']);
                
                if energy==1
                    b.st.coupl_counter.th_e1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)    = eval(['c.' b.bc.z{n} '.st.th_e']);

                end

            end
        end

%%%%%%%%%%%%%%%%%%%%%%% upper boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % surface, interior or periodic
        if (strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'s')==1 || strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'c')==1 || strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'p')==1) && max(strcmp(b.bc.z{length(b.bc.z)/2+n},c.bc.z))==1
            if stability~=1
                if strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'i')~=1 && strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'r')~=1
                    b.coupl_counter.rmend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)      = [eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.rm']) sparse(J,c.sizeR1 + c.sizeRend)];
                    b.Ja.coupl_counter.rmend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)   = [eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.Ja.rm']) sparse(J,c.sizeR1 + c.sizeRend)];
                    b.rhs_counter.rmend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),1)        = eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.rhs.rm']);

                end
                b.coupl_counter.zmend(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)      = [eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.zm']) sparse(J-1,c.sizeR1 + c.sizeRend)];
                b.Ja.coupl_counter.zmend(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)   = [eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.Ja.zm']) sparse(J-1,c.sizeR1 + c.sizeRend)];
                b.rhs_counter.zmend(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),1)     = eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.rhs.zm']);
                
                if energy==1
                    b.coupl_counter.th_eend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)    = [eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.th_e']) sparse(J,c.sizeR1 + c.sizeRend)];
                    b.Ja.coupl_counter.th_eend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) = [eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.Ja.th_e']) sparse(J,c.sizeR1 + c.sizeRend)];
                    
                end
                
            else
                b.st.coupl_counter.phimend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)    = eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.st.phim']);
                b.st.coupl_counter.zmend(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)      = eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.st.zm']);
                
                if energy==1
                    b.st.coupl_counter.th_eend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)    = eval(['c.' b.bc.z{length(b.bc.z)/2+n} '.st.th_e']);
                    
                end
                
            end
            
        end
        
    end
    
    if stability~=1
        if energy==1
            b.coup_counter=sparse((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1 + b.sizeRend, (c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1 + c.sizeRend);
            b.Ja.coup_counter=b.coup_counter;

            b.coup_counter(2:(b.J+1),:)                                                                                           = b.coupl_counter.rm1(2:end-1,:);
            b.coup_counter(b.K*(b.J+2)+2:(b.K+1)*(b.J+2)-1,:)                                                                     = b.coupl_counter.rmend(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+2:(b.K+1)*(b.J+2)+b.J,:)                                                               = b.coupl_counter.zm1(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+(b.K+1)*(b.J+1)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)-1,:)                                 = b.coupl_counter.zmend(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+b.J+1,:)                             = b.coupl_counter.th_e1(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+(b.K+1)*(b.J+2)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+(b.K+2)*(b.J+2)-1,:) = b.coupl_counter.th_eend(2:end-1,:);

            b.Ja.coup_counter(2:(b.J+1),:)                                                                                           = b.Ja.coupl_counter.rm1(2:end-1,:);
            b.Ja.coup_counter(b.K*(b.J+2)+2:(b.K+1)*(b.J+2)-1,:)                                                                     = b.Ja.coupl_counter.rmend(2:end-1,:);
            b.Ja.coup_counter((b.K+1)*(b.J+2)+2:(b.K+1)*(b.J+2)+b.J,:)                                                               = b.Ja.coupl_counter.zm1(2:end-1,:);
            b.Ja.coup_counter((b.K+1)*(b.J+2)+(b.K+1)*(b.J+1)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)-1,:)                                 = b.Ja.coupl_counter.zmend(2:end-1,:);
            b.Ja.coup_counter((b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+b.J+1,:)                             = b.Ja.coupl_counter.th_e1(2:end-1,:);
            b.Ja.coup_counter((b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+(b.K+1)*(b.J+2)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)+(b.K+2)*(b.J+2)-1,:) = b.Ja.coupl_counter.th_eend(2:end-1,:);
            if isempty(b.Ja_counter.r_m.Rsend)==0
                b.coup_counter((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + 1:(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1)=b.counter.C.Rsend;
                
                b.Ja.coup_counter(1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + 1:(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1)=[b.Ja_counter.r_m.Rsend; b.Ja_counter.z_m.Rsend; b.Ja_counter.th_e.Rsend; b.Ja_counter.C.Rsend];

            end
            if isempty(b.Ja_counter.r_m.Rs1)==0
                b.coup_counter((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1 + 1:(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1 + c.sizeRend)=b.counter.C.Rs1;
                
                b.Ja.coup_counter(1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1 + 1:(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2) + c.sizeR1 + c.sizeRend)=[b.Ja_counter.r_m.Rs1; b.Ja_counter.z_m.Rs1; b.Ja_counter.th_e.Rs1; b.Ja_counter.C.Rs1];
                
            end

            b.rhs=b.rhs+[b.rhs_counter.rm1; sparse((b.K-1)*(b.J+2),1); b.rhs_counter.rmend; b.rhs_counter.zm1; sparse((b.K)*(b.J+1),1); b.rhs_counter.zmend; sparse((b.K+2)*(b.J+2),1); sparse(b.K*b.J,1); sparse(b.sizeR1 + b.sizeRend,1)];
            
        else
            b.coup_counter=sparse((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + b.sizeRend, (c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + c.sizeR1 + c.sizeRend);
            b.Ja.coup_counter=b.coup_counter;

            b.coup_counter(2:(b.J+1),:)                                                                                           = b.coupl_counter.rm1(2:end-1,:);
            b.coup_counter(b.K*(b.J+2)+2:(b.K+1)*(b.J+2)-1,:)                                                                     = b.coupl_counter.rmend(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+2:(b.K+1)*(b.J+2)+b.J,:)                                                               = b.coupl_counter.zm1(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+(b.K+1)*(b.J+1)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)-1,:)                                 = b.coupl_counter.zmend(2:end-1,:);
            
            b.Ja.coup_counter(2:(b.J+1),:)                                                                                           = b.Ja.coupl_counter.rm1(2:end-1,:);
            b.Ja.coup_counter(b.K*(b.J+2)+2:(b.K+1)*(b.J+2)-1,:)                                                                     = b.Ja.coupl_counter.rmend(2:end-1,:);
            b.Ja.coup_counter((b.K+1)*(b.J+2)+2:(b.K+1)*(b.J+2)+b.J,:)                                                               = b.Ja.coupl_counter.zm1(2:end-1,:);
            b.Ja.coup_counter((b.K+1)*(b.J+2)+(b.K+1)*(b.J+1)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+1)-1,:)                                 = b.Ja.coupl_counter.zmend(2:end-1,:);
            if isempty(b.Ja_counter.r_m.Rsend)==0
                b.coup_counter((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + 1:(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + c.sizeR1)=b.counter.C.Rsend;
                
                b.Ja.coup_counter(1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + 1:(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + c.sizeR1)=[b.Ja_counter.r_m.Rsend; b.Ja_counter.z_m.Rsend; b.Ja_counter.C.Rsend];

            end
            if isempty(b.Ja_counter.r_m.Rs1)==0
                b.coup_counter((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + c.sizeR1 + 1:(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + c.sizeR1 + c.sizeRend)=b.counter.C.Rs1;
                
                b.Ja.coup_counter(1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J,(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + c.sizeR1 + 1:(c.K+1)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + c.sizeR1 + c.sizeRend)=[b.Ja_counter.r_m.Rs1; b.Ja_counter.z_m.Rs1; b.Ja_counter.C.Rs1];

            end

            b.rhs=b.rhs+[b.rhs_counter.rm1; sparse((b.K-1)*(b.J+2),1); b.rhs_counter.rmend; b.rhs_counter.zm1; sparse((b.K)*(b.J+1),1); b.rhs_counter.zmend; sparse(b.K*b.J,1); sparse(b.sizeR1 + b.sizeRend,1)];
            
        end
        
    else
        if energy==1
            b.coup_counter=sparse((b.K+1)*(b.J+2) + (b.K+2)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J, (c.K+1)*(c.J+2) + (c.K+2)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J + (c.K+2)*(c.J+2));
            
            b.coup_counter((b.K+1)*(b.J+2)+2:(b.K+1)*(b.J+2)+b.J+1,:)                                                                                             = b.st.coupl_counter.phim1(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+(b.K+1)*(b.J+2)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)-1,:)                                                                 = b.st.coupl_counter.phimend(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+b.J,:)                                                               = b.st.coupl_counter.zm1(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+1)*(b.J+1)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)-1,:)                                 = b.st.coupl_counter.zmend(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)+b.J+1,:)                             = b.st.coupl_counter.th_e1(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)+(b.K+1)*(b.J+2)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)+(b.K+2)*(b.J+2)-1,:) = b.st.coupl_counter.th_eend(2:end-1,:);

        else
            b.coup_counter=sparse((b.K+1)*(b.J+2) + (b.K+2)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J, (c.K+1)*(c.J+2) + (c.K+2)*(c.J+2) + (c.K+2)*(c.J+1) + c.K*c.J);
            
            b.coup_counter((b.K+1)*(b.J+2)+2:(b.K+1)*(b.J+2)+b.J+1,:)                                                                                             = b.st.coupl_counter.phim1(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+(b.K+1)*(b.J+2)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)-1,:)                                                                 = b.st.coupl_counter.phimend(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+b.J,:)                                                               = b.st.coupl_counter.zm1(2:end-1,:);
            b.coup_counter((b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+1)*(b.J+1)+2:(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)-1,:)                                 = b.st.coupl_counter.zmend(2:end-1,:);
            
        end
        
    end