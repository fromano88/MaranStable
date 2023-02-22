function[b]=block_coupling(b, T_fluid, flowopt)
    ax=flowopt.ax;
    g=flowopt.g;
    m=flowopt.m;
    energy=flowopt.energy;
    stability=flowopt.stability;
    thermcapcon=flowopt.thermcapcon;
    iv=flowopt.iv;
    
    rho=str2func(b.rho);
    mu=str2func(b.mu);
    lambda=str2func(b.lambda);
    drho=str2func(b.drho);
    dmu=str2func(b.dmu);
    dlambda=str2func(b.dlambda);
    if isfield(b,'sigma')
        sigma=str2func(b.sigma);
        dsigma=str2func(b.dsigma);
        
    end
    
    for o=1:2
        % normal stresses
            block_coupling_normal_stresses

        % axial shear stresses
            block_coupling_axial_shear_stresses

        % normal heat transfer
            if energy==1
                block_coupling_normal_heat_transfer

            end

        % azimuthal shear stresses
            if stability==1
                block_coupling_azimuthal_shear_stresses

            end
            
    end
            
    coup.rm1.u=sparse(b.J+2,(b.K+1)*(b.J+2));
    coup.rm1.w=sparse(b.J+2,(b.K+2)*(b.J+1));
    coup.rm1.p=sparse(b.J+2,b.K*b.J);
    coup.zm1.u=sparse(b.J+1,(b.K+1)*(b.J+2));
    coup.zm1.w=sparse(b.J+1,(b.K+2)*(b.J+1));
    
    rhs.rm1=sparse(b.J+2,1);
    rhs.zm1=sparse(b.J+1,1);
    
    coup.rmend.u=sparse(b.J+2,(b.K+1)*(b.J+2));
    coup.rmend.w=sparse(b.J+2,(b.K+2)*(b.J+1));
    coup.rmend.p=sparse(b.J+2,b.K*b.J);
    coup.zmend.u=sparse(b.J+1,(b.K+1)*(b.J+2));
    coup.zmend.w=sparse(b.J+1,(b.K+2)*(b.J+1));
    
    rhs.rmend=sparse(b.J+2,1);
    rhs.zmend=sparse(b.J+1,1);
    
    rhs.p=sparse(b.K*b.J,1);
    
    if energy==1
        coup.th_e1.T=sparse(b.J+2,(b.K+2)*(b.J+2));
        coup.th_eend.T=sparse(b.J+2,(b.K+2)*(b.J+2));
        
        coup.Ja.rm1.T=sparse(b.J+2,(b.K+2)*(b.J+2));
        coup.Ja.zm1.T=sparse(b.J+1,(b.K+2)*(b.J+2));
        coup.Ja.th_e1.T=sparse(b.J+2,(b.K+2)*(b.J+2));

        coup.Ja.rmend.T=sparse(b.J+2,(b.K+2)*(b.J+2));
        coup.Ja.zmend.T=sparse(b.J+1,(b.K+2)*(b.J+2));
        coup.Ja.th_eend.T=sparse(b.J+2,(b.K+2)*(b.J+2));
        
    end
    
    if stability~=1
        coup.Ja.rm1.Rs=[];
        coup_counter.Ja.rm1.Rs=[];
        coup.Ja.zm1.Rs=[];
        coup_counter.Ja.zm1.Rs=[];
        coup.Ja.rmend.Rs=[];
        coup_counter.Ja.rmend.Rs=[];
        coup.Ja.zmend.Rs=[];
        coup_counter.Ja.zmend.Rs=[];
        
        if energy==1
            coup.Ja.th_e1.Rs=[];
            coup_counter.Ja.th_e1.Rs=[];
            coup.Ja.th_eend.Rs=[];
            coup_counter.Ja.th_eend.Rs=[];

        end
        
    else
        coup.st.phim1.u=sparse(b.J+2,(b.K+1)*(b.J+2));
        coup.st.phim1.v=sparse(b.J+2,(b.K+2)*(b.J+2));
        coup.st.phimend.u=sparse(b.J+2,(b.K+1)*(b.J+2));
        coup.st.phimend.v=sparse(b.J+2,(b.K+2)*(b.J+2));
        
        if energy==1            
            coup.st.phim1.T=sparse(b.J+2,(b.K+2)*(b.J+2));
            coup.st.phimend.T=sparse(b.J+2,(b.K+2)*(b.J+2));
    
        end
        
    end

    for n=1:length(b.bc.z)/2
        J=length(b.z.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)));
        
%%%%%%%%%%%%%%%%%%%%%%% lower boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % surface, interior or periodic
        if strcmp(b.bc.z{n}(1),'s')==1 || strcmp(b.bc.z{n}(1),'c')==1 || strcmp(b.bc.z{n}(1),'p')==1
            if eval(['b.' b.bc.z{n} '.o==1'])
                % r-momentum
                    % isothermal, no influence of shear stresses on static surface deformation or rigid cylindrical surface
                        if strcmp(b.bc.z{n}(2),'i')==1 || strcmp(b.bc.z{n}(2),'r')==1 || stability==1
                            coup.rm1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.e1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                    % static surface deformation, interior or periodic
                        else
                            coup.rm1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.ns1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                            coup.rm1.w(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.ns1.w(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                            coup.rm1.p(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.ns1.p(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                            ns.Rs1=Ja.ns.Rs1(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                            if n>1
                                ns.Rs1(length(b.z.u(b.z.u<b.geom.z(n))),1)=0;
                                ns.Rs1(length(b.z.u(b.z.u<b.geom.z(n)))+1,1)=ns.Rs1(length(b.z.u(b.z.u<b.geom.z(n)))+1,1)-Ja.ns.Rs1_ad_W(length(b.z.u(b.z.u<b.geom.z(n))));
                                
                            end
                            if n<length(b.bc.z)/2
                                ns.Rs1(length(b.z.u(b.z.u<b.geom.z(n+1)))+1,end)=0;
                                ns.Rs1(length(b.z.u(b.z.u<b.geom.z(n+1))),end)=ns.Rs1(length(b.z.u(b.z.u<b.geom.z(n+1))),end)-Ja.ns.Rs1_ad_E(length(b.z.u(b.z.u<b.geom.z(n+1)))-1);
                                
                            end
                                                        
                            coup.Ja.rm1.Rs=[coup.Ja.rm1.Rs ns.Rs1];
                            
                            if energy==1
                                coup.Ja.rm1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=Ja.ns.T1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                                
                                eval(['b.' b.bc.z{n} '.rm=-[b.e1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2))];'])
                            
                            else
                                eval(['b.' b.bc.z{n} '.rm=-[b.e1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+1) + b.K*b.J)];'])
                            
                            end

                            eval(['b.' b.bc.z{n} '.Ja.rm=b.' b.bc.z{n} '.rm;'])
                            eval(['b.' b.bc.z{n} '.rhs.rm=sparse(length(b.z.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1))),1);'])

                            % surface
                                if strcmp(b.bc.z{n}(1),'s')==1 && isfield(b, 'sigma')==1
                                    rhs.rm1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),1)=-HP1(1,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)).'-LP1(1,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)).';

                            % interior or periodic
                                elseif strcmp(b.bc.z{n}(1),'s')==1 || strcmp(b.bc.z{n}(1),'c')==1 || strcmp(b.bc.z{n}(1),'p')==1
                                    rhs.rm1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),1)=-HP1(1,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)).';

                                end

                        end
                        
                % phi-momentum
                    if stability==1
                        coup.st.phim1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=st.ass1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                        coup.st.phim1.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=st.ass1.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                        
                        if energy==1
                            coup.st.phim1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=st.ass1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                            eval(['b.' b.bc.z{n} '.st.phim=-[sparse(J,(b.K+1)*(b.J+2)) b.e1.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2))];'])
                            
                        else
                            eval(['b.' b.bc.z{n} '.st.phim=-[sparse(J,(b.K+1)*(b.J+2)) b.e1.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+1) + b.K*b.J)];'])

                        end
                        
                    end                

                % z-momentum
                    coup.zm1.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)=b.ss1.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:);
                    coup.zm1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)=b.ss1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:);
                    
                    if stability~=1
                        if energy==1
                            eval(['b.' b.bc.z{n} '.zm=-[sparse(J-1,(b.K+1)*(b.J+2)) b.e1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J + (b.K+2)*(b.J+2))];'])

                            if strcmp(b.bc.z{n}(2),'s')==1
                                ss.Rs1=Ja.ss.Rs1(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                                if n>1
                                    ss.Rs1(find(b.z.w==b.geom.z(n))-1,1)=0;
                                    ss.Rs1(b.z.w==b.geom.z(n),1:2)=0;
                                    ss.Rs1(find(b.z.w==b.geom.z(n))+1,1)=ss.Rs1(find(b.z.w==b.geom.z(n))+1,1)-Ja.ss.Rs1_ad_w(b.z.w==b.geom.z(n));

                                end
                                if n<length(b.bc.z)/2
                                    ss.Rs1(b.z.w==b.geom.z(n+1),end-1:end)=0;
                                    ss.Rs1(find(b.z.w==b.geom.z(n+1))+1,end)=0;
                                    ss.Rs1(find(b.z.w==b.geom.z(n+1))-1,end)=ss.Rs1(find(b.z.w==b.geom.z(n+1))-1,end)-Ja.ss.Rs1_ad_e(find(b.z.w==b.geom.z(n+1))-2);

                                end
                                
                                coup.Ja.zm1.Rs=[coup.Ja.zm1.Rs ss.Rs1];
                                                            
                            end
                            coup.Ja.zm1.T(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)=Ja.ss.T1(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:);
                            rhs.zm1(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),1)=-THS1(1,b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1)).';

                % thermal energy
                            coup.th_e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.nht1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                            if strcmp(b.bc.z{n}(2),'s')==1
                                nht.Rs1=Ja.nht.Rs1(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                                if n>1
                                    nht.Rs1(b.z.w==b.geom.z(n),1)=0;
                                    nht.Rs1(find(b.z.w==b.geom.z(n))+1,1)=nht.Rs1(find(b.z.w==b.geom.z(n))+1,1)-Ja.nht.Rs1_ad_w(b.z.w==b.geom.z(n));
                                    if  ((strcmp(b.bc.z{n-1}(1),'w')==1 || strcmp(b.bc.z{n-1}(1),'i')==1 || strcmp(b.bc.z{n-1}(1),'o')==1) && strcmp(b.bc.z{n-1}(3),'a')==1)
                                        nht.Rs1(b.z.w==b.geom.z(n),1)=Ja.ad.Rs1_E(b.z.w==b.geom.z(n));

                                    end

                                end
                                if n<length(b.bc.z)/2
                                    nht.Rs1(find(b.z.w==b.geom.z(n+1))+1,end)=0;
                                    nht.Rs1(b.z.w==b.geom.z(n+1),end)=nht.Rs1(b.z.w==b.geom.z(n+1),end)-Ja.nht.Rs1_ad_e(find(b.z.w==b.geom.z(n+1))-1);
                                    if  ((strcmp(b.bc.z{n+1}(1),'w')==1 || strcmp(b.bc.z{n+1}(1),'i')==1 || strcmp(b.bc.z{n+1}(1),'o')==1) && strcmp(b.bc.z{n+1}(3),'a')==1)
                                        nht.Rs1(find(b.z.w==b.geom.z(n+1))+1,end)=Ja.ad.Rs1_W(find(b.z.w==b.geom.z(n+1))+1);

                                    end

                                end

                                coup.Ja.th_e1.Rs=[coup.Ja.th_e1.Rs nht.Rs1];

                            end

                            coup.Ja.th_e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=Ja.nht.T1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                            eval(['b.' b.bc.z{n} '.th_e=-[sparse(J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J) b.e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])
                            eval(['b.' b.bc.z{n} '.Ja.th_e=b.' b.bc.z{n} '.th_e;'])

                        else
                            if strcmp(b.bc.z{n}(2),'s')==1
                                ss.Rs1=Ja.ss.Rs1(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                                if n>1
                                    ss.Rs1(find(b.z.w==b.geom.z(n))-1,1)=0;
                                    ss.Rs1(b.z.w==b.geom.z(n),1:2)=0;
                                    ss.Rs1(find(b.z.w==b.geom.z(n))+1,1)=ss.Rs1(find(b.z.w==b.geom.z(n))+1,1)-Ja.ss.Rs1_ad_w(b.z.w==b.geom.z(n));

                                end
                                if n<length(b.bc.z)/2
                                    ss.Rs1(b.z.w==b.geom.z(n+1),end-1:end)=0;
                                    ss.Rs1(find(b.z.w==b.geom.z(n+1))+1,end)=0;
                                    ss.Rs1(find(b.z.w==b.geom.z(n+1))-1,end)=ss.Rs1(find(b.z.w==b.geom.z(n+1))-1,end)-Ja.ss.Rs1_ad_e(find(b.z.w==b.geom.z(n+1))-2);

                                end
                                
                                coup.Ja.zm1.Rs=[coup.Ja.zm1.Rs ss.Rs1];
                        
                            end
                            eval(['b.' b.bc.z{n} '.zm=-[sparse(J-1,(b.K+1)*(b.J+2)) b.e1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J)];'])

                        end

                        eval(['b.' b.bc.z{n} '.Ja.zm=b.' b.bc.z{n} '.zm;'])
                        eval(['b.' b.bc.z{n} '.rhs.zm=sparse(length(b.z.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1))),1);'])

                    else
                        if energy==1
                            coup.Ja.zm1.T(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)=Ja.ss.T1(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:);
                                                                
                            eval(['b.' b.bc.z{n} '.st.zm=-[sparse(J-1,(b.K+1)*(b.J+2)+(b.K+2)*(b.J+2)) b.e1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J + (b.K+2)*(b.J+2))];'])

                % thermal energy
                                coup.th_e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.nht1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                                coup.Ja.th_e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=Ja.nht.T1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                                eval(['b.' b.bc.z{n} '.st.th_e=-[sparse(J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J) b.e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])

                        else
                            eval(['b.' b.bc.z{n} '.st.zm=-[sparse(J-1,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)) b.e1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J)];'])

                        end

                    end
                
            else
                
                % r-momentum
                    % isothermal, no influence of shear stresses on static surface deformation or rigid cylindrical surface
                        if strcmp(b.bc.z{n}(2),'i')==1 || strcmp(b.bc.z{n}(2),'r')==1 || stability==1
                            coup.rm1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.e1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                    % static surface deformation, interior or periodic
                        else
                            coup.rm1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.e1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                            if energy==1
                                eval(['b.' b.bc.z{n} '.rm=-[b.ns1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) b.ns1.w(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) b.ns1.p(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+2))];'])

                                eval(['b.' b.bc.z{n} '.Ja.rm=-[b.ns1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) b.ns1.w(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) b.ns1.p(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) Ja.ns.T1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])
                                
                            else
                                eval(['b.' b.bc.z{n} '.rm=-[b.ns1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) b.ns1.w(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) b.ns1.p(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])

                                eval(['b.' b.bc.z{n} '.Ja.rm=b.' b.bc.z{n} '.rm;'])
                                
                            end
                            
                            ns.Rs1=Ja.ns.Rs1(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                            if n>1
                                ns.Rs1(length(b.z.u(b.z.u<b.geom.z(n))),1)=0;
                                ns.Rs1(length(b.z.u(b.z.u<b.geom.z(n)))+1,1)=ns.Rs1(length(b.z.u(b.z.u<b.geom.z(n)))+1,1)-Ja.ns.Rs1_ad_W(length(b.z.u(b.z.u<b.geom.z(n))));
                                
                            end
                            if n<length(b.bc.z)/2
                                ns.Rs1(length(b.z.u(b.z.u<b.geom.z(n+1)))+1,end)=0;
                                ns.Rs1(length(b.z.u(b.z.u<b.geom.z(n+1))),end)=ns.Rs1(length(b.z.u(b.z.u<b.geom.z(n+1))),end)-Ja.ns.Rs1_ad_E(length(b.z.u(b.z.u<b.geom.z(n+1)))-1);
                                
                            end
                               
                            eval(['b.' b.bc.z{n} '.Ja.rm1.Rs=ns.Rs1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);'])
                            
                            % surface
                                if strcmp(b.bc.z{n}(1),'s')==1 && isfield(b, 'sigma')==1
                                    eval(['b.' b.bc.z{n} '.rhs.rm=HP1(1,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)).''+LP1(1,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)).'';'])

                            % surface, interior or periodic
                                elseif strcmp(b.bc.z{n}(1),'s')==1 || strcmp(b.bc.z{n}(1),'c')==1 || strcmp(b.bc.z{n}(1),'p')==1
                                    eval(['b.' b.bc.z{n} '.rhs.rm=HP1(1,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)).'';'])

                                end

                        end
                        
                % phi-momentum
                    if stability==1
                        coup.st.phim1.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.e1.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                    
                        if energy==1
                            eval(['b.' b.bc.z{n} '.st.phim=-[st.ass1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) st.ass1.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+1)+b.K*b.J) st.ass1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])

                        else
                            eval(['b.' b.bc.z{n} '.st.phim=-[st.ass1.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) st.ass1.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+1)+b.K*b.J)];'])

                        end
                        
                    end                       
                        
                % z-momentum
                    coup.zm1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)=b.e1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:);
                    
                    if stability~=1
                        if energy==1
                            eval(['b.' b.bc.z{n} '.zm=-[b.ss1.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) b.ss1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J + (b.K+2)*(b.J+2))];'])

                            eval(['b.' b.bc.z{n} '.Ja.zm=-[b.ss1.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) b.ss1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J) Ja.ss.T1(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)];'])
                            if strcmp(b.bc.z{n}(2),'s')==1
                                ss.Rs1=Ja.ss.Rs1(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                                if n>1
                                    ss.Rs1(find(b.z.w==b.geom.z(n))-1,1)=0;
                                    ss.Rs1(b.z.w==b.geom.z(n),1:2)=0;
                                    ss.Rs1(find(b.z.w==b.geom.z(n))+1,1)=ss.Rs1(find(b.z.w==b.geom.z(n))+1,1)-Ja.ss.Rs1_ad_w(b.z.w==b.geom.z(n));

                                end
                                if n<length(b.bc.z)/2
                                    ss.Rs1(b.z.w==b.geom.z(n+1),end-1:end)=0;
                                    ss.Rs1(find(b.z.w==b.geom.z(n+1))+1,end)=0;
                                    ss.Rs1(find(b.z.w==b.geom.z(n+1))-1,end)=ss.Rs1(find(b.z.w==b.geom.z(n+1))-1,end)-Ja.ss.Rs1_ad_e(find(b.z.w==b.geom.z(n+1))-2);

                                end
                                
                                eval(['b.' b.bc.z{n} '.Ja.zm1.Rs=ss.Rs1(b.z.w>=b.geom.z(n) & b.z.w<=b.geom.z(n+1),:);'])

                            end
                            eval(['b.' b.bc.z{n} '.rhs.zm=THS1(1,b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1)).'';'])

                % thermal energy
                            coup.th_e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                            coup.Ja.th_e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=coup.th_e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                            
                            eval(['b.' b.bc.z{n} '.th_e=-[sparse(J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J) b.nht1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])

                            if strcmp(b.bc.z{n}(2),'s')==1
                                nht.Rs1=Ja.nht.Rs1(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                                if n>1
                                    nht.Rs1(b.z.w==b.geom.z(n),1)=0;
                                    nht.Rs1(find(b.z.w==b.geom.z(n))+1,1)=nht.Rs1(find(b.z.w==b.geom.z(n))+1,1)-Ja.nht.Rs1_ad_w(b.z.w==b.geom.z(n));
                                    if  ((strcmp(b.bc.z{n-1}(1),'w')==1 || strcmp(b.bc.z{n-1}(1),'i')==1 || strcmp(b.bc.z{n-1}(1),'o')==1) && strcmp(b.bc.z{n-1}(3),'a')==1)
                                        nht.Rs1(b.z.w==b.geom.z(n),1)=Ja.ad.Rs1_E(b.z.w==b.geom.z(n));

                                    end

                                end
                                if n<length(b.bc.z)/2
                                    nht.Rs1(find(b.z.w==b.geom.z(n+1))+1,end)=0;
                                    nht.Rs1(b.z.w==b.geom.z(n+1),end)=nht.Rs1(b.z.w==b.geom.z(n+1),end)-Ja.nht.Rs1_ad_e(find(b.z.w==b.geom.z(n+1))-1);
                                    if  ((strcmp(b.bc.z{n+1}(1),'w')==1 || strcmp(b.bc.z{n+1}(1),'i')==1 || strcmp(b.bc.z{n+1}(1),'o')==1) && strcmp(b.bc.z{n+1}(3),'a')==1)
                                        nht.Rs1(find(b.z.w==b.geom.z(n+1))+1,end)=Ja.ad.Rs1_W(find(b.z.w==b.geom.z(n+1))+1);

                                    end

                                end

                                eval(['b.' b.bc.z{n} '.Ja.th_e1.Rs=nht.Rs1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);'])
                                
                            end

                            eval(['b.' b.bc.z{n} '.Ja.th_e=-[sparse(J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J) Ja.nht.T1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])

                        else
                            eval(['b.' b.bc.z{n} '.zm=-[b.ss1.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) b.ss1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J)];'])

                            eval(['b.' b.bc.z{n} '.Ja.zm=b.' b.bc.z{n} '.zm;'])

                            if strcmp(b.bc.z{n}(2),'s')==1
                                ss.Rs1=Ja.ss.Rs1(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                                if n>1
                                    ss.Rs1(find(b.z.w==b.geom.z(n))-1,1)=0;
                                    ss.Rs1(b.z.w==b.geom.z(n),1:2)=0;
                                    ss.Rs1(find(b.z.w==b.geom.z(n))+1,1)=ss.Rs1(find(b.z.w==b.geom.z(n))+1,1)-Ja.ss.Rs1_ad_w(b.z.w==b.geom.z(n));

                                end
                                if n<length(b.bc.z)/2
                                    ss.Rs1(b.z.w==b.geom.z(n+1),end-1:end)=0;
                                    ss.Rs1(find(b.z.w==b.geom.z(n+1))+1,end)=0;
                                    ss.Rs1(find(b.z.w==b.geom.z(n+1))-1,end)=ss.Rs1(find(b.z.w==b.geom.z(n+1))-1,end)-Ja.ss.Rs1_ad_e(find(b.z.w==b.geom.z(n+1))-2);

                                end
                                
                                eval(['b.' b.bc.z{n} '.Ja.zm1.Rs=ss.Rs1(b.z.w>=b.geom.z(n) & b.z.w<=b.geom.z(n+1),:);'])
                                
                            end
                                
                            eval(['b.' b.bc.z{n} '.rhs.zm=sparse(length(b.z.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1))),1);'])

                        end
                        
                    else
                        if energy==1
                            eval(['b.' b.bc.z{n} '.st.zm=-[b.ss1.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,(b.K+2)*(b.J+2)) b.ss1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J) Ja.ss.T1(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)];'])
                                
                % thermal energy
                            coup.th_e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                            coup.Ja.th_e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=coup.th_e1.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                            
                            eval(['b.' b.bc.z{n} '.st.th_e=-[sparse(J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J) Ja.nht.T1(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])

                        else
                            eval(['b.' b.bc.z{n} '.st.zm=-[b.ss1.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,(b.K+2)*(b.J+2)) b.ss1.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J)];'])

                        end
                        
                    end
                                        
            end

        end

%%%%%%%%%%%%%%%%%%%%%%% upper boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % surface, interior or periodic
        if strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'s')==1 || strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'c')==1 || strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'p')==1
            if eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.o==1'])
                % r-momentum
                    % isothermal, no influence of shear stresses on static surface deformation or rigid cylindrical surface
                        if strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'i')==1 || strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'r')==1 || stability==1
                            coup.rmend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.eend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                    % static surface deformation, interior or periodic
                        else
                            coup.rmend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.nsend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                            coup.rmend.w(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.nsend.w(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                            coup.rmend.p(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.nsend.p(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                            ns.Rsend=Ja.ns.Rsend(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                            if n>1
                                ns.Rsend(length(b.z.u(b.z.u<b.geom.z(n))),1)=0;
                                ns.Rsend(length(b.z.u(b.z.u<b.geom.z(n)))+1,1)=ns.Rsend(length(b.z.u(b.z.u<b.geom.z(n)))+1,1)-Ja.ns.Rsend_ad_W(length(b.z.u(b.z.u<b.geom.z(n))));
                                
                            end
                            if n<length(b.bc.z)/2
                                ns.Rsend(length(b.z.u(b.z.u<b.geom.z(n+1)))+1,end)=0;
                                ns.Rsend(length(b.z.u(b.z.u<b.geom.z(n+1))),end)=ns.Rsend(length(b.z.u(b.z.u<b.geom.z(n+1))),end)-Ja.ns.Rsend_ad_E(length(b.z.u(b.z.u<b.geom.z(n+1)))-1);
                                
                            end
                                                        
                            coup.Ja.rmend.Rs=[coup.Ja.rmend.Rs ns.Rsend];
                            
                            if energy==1
                                coup.Ja.rmend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=Ja.ns.Tend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                                eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.rm=-[b.eend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2))];'])
                            
                            else
                                eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.rm=-[b.eend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+1) + b.K*b.J)];'])
                            
                            end
                            
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.Ja.rm=b.' b.bc.z{length(b.bc.z)/2+n} '.rm;'])
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.rhs.rm=sparse(length(b.z.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1))),1);'])

                            % surface
                                if strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'s')==1
                                    rhs.rmend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),1)=-HPend(1,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)).'-LPend(1,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)).';

                            % interior or periodic
                                elseif strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'c')==1 || strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'p')==1
                                    rhs.rmend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),1)=-HPend(1,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)).';

                                end

                        end

                % phi-momentum
                    if stability==1
                        coup.st.phimend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=st.assend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                        coup.st.phimend.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=st.assend.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                        if energy==1
                            coup.st.phimend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=st.assend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.st.phim=-[sparse(J,(b.K+1)*(b.J+2)) b.eend.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2))];'])
                            
                        else
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.st.phim=-[sparse(J,(b.K+1)*(b.J+2)) b.eend.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+1) + b.K*b.J)];'])

                        end
                        
                    end                

                % z-momentum
                    coup.zmend.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)=b.ssend.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:);
                    coup.zmend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)=b.ssend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:);

                    if stability~=1
                        if energy==1
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.zm=-[sparse(J-1,(b.K+1)*(b.J+2)) b.eend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J + (b.K+2)*(b.J+2))];'])
                    
                            if strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'s')==1
                                ss.Rsend=Ja.ss.Rsend(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                                if n>1
                                    ss.Rsend(find(b.z.w==b.geom.z(n))-1,1)=0;
                                    ss.Rsend(b.z.w==b.geom.z(n),1:2)=0;
                                    ss.Rsend(find(b.z.w==b.geom.z(n))+1,1)=ss.Rsend(find(b.z.w==b.geom.z(n))+1,1)-Ja.ss.Rsend_ad_w(b.z.w==b.geom.z(n));

                                end
                                if n<length(b.bc.z)/2
                                    ss.Rsend(b.z.w==b.geom.z(n+1),end-1:end)=0;
                                    ss.Rsend(find(b.z.w==b.geom.z(n+1))+1,end)=0;
                                    ss.Rsend(find(b.z.w==b.geom.z(n+1))-1,end)=ss.Rsend(find(b.z.w==b.geom.z(n+1))-1,end)-Ja.ss.Rsend_ad_e(find(b.z.w==b.geom.z(n+1))-2);

                                end
                                
                                coup.Ja.zmend.Rs=[coup.Ja.zmend.Rs ss.Rsend];

                            end
                            coup.Ja.zmend.T(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)=Ja.ss.Tend(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:);
                            rhs.zmend(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),1)=-THSend(1,b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1)).';

                % thermal energy
                            coup.th_eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.nhtend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                            if strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'s')==1
                                nht.Rsend=Ja.nht.Rsend(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                                if n>1
                                    nht.Rsend(b.z.w==b.geom.z(n),1)=0;
                                    nht.Rsend(find(b.z.w==b.geom.z(n))+1,1)=nht.Rsend(find(b.z.w==b.geom.z(n))+1,1)-Ja.nht.Rsend_ad_w(b.z.w==b.geom.z(n));
                                    if  ((strcmp(b.bc.z{n-1}(1),'w')==1 || strcmp(b.bc.z{n-1}(1),'i')==1 || strcmp(b.bc.z{n-1}(1),'o')==1) && strcmp(b.bc.z{n-1}(3),'a')==1)
                                        nht.Rsend(b.z.w==b.geom.z(n),1)=Ja.ad.Rsend_E(b.z.w==b.geom.z(n));

                                    end

                                end
                                if n<length(b.bc.z)/2
                                    nht.Rsend(find(b.z.w==b.geom.z(n+1))+1,end)=0;
                                    nht.Rsend(b.z.w==b.geom.z(n+1),end)=nht.Rsend(b.z.w==b.geom.z(n+1),end)-Ja.nht.Rsend_ad_e(find(b.z.w==b.geom.z(n+1))-1);
                                    if  ((strcmp(b.bc.z{n+1}(1),'w')==1 || strcmp(b.bc.z{n+1}(1),'i')==1 || strcmp(b.bc.z{n+1}(1),'o')==1) && strcmp(b.bc.z{n+1}(3),'a')==1)
                                        nht.Rsend(find(b.z.w==b.geom.z(n+1))+1,end)=Ja.ad.Rsend_W(find(b.z.w==b.geom.z(n+1))+1);

                                    end

                                end

                                coup.Ja.th_eend.Rs=[coup.Ja.th_eend.Rs nht.Rsend];

                            end

                            coup.Ja.th_eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=Ja.nht.Tend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.th_e=-[sparse(J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J) b.eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.Ja.th_e=b.' b.bc.z{length(b.bc.z)/2+n} '.th_e;'])

                        else
                            if strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'s')==1
                                ss.Rsend=Ja.ss.Rsend(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                                if n>1
                                    ss.Rsend(find(b.z.w==b.geom.z(n))-1,1)=0;
                                    ss.Rsend(b.z.w==b.geom.z(n),1:2)=0;
                                    ss.Rsend(find(b.z.w==b.geom.z(n))+1,1)=ss.Rsend(find(b.z.w==b.geom.z(n))+1,1)-Ja.ss.Rsend_ad_w(b.z.w==b.geom.z(n));

                                end
                                if n<length(b.bc.z)/2
                                    ss.Rsend(b.z.w==b.geom.z(n+1),end-1:end)=0;
                                    ss.Rsend(find(b.z.w==b.geom.z(n+1))+1,end)=0;
                                    ss.Rsend(find(b.z.w==b.geom.z(n+1))-1,end)=ss.Rsend(find(b.z.w==b.geom.z(n+1))-1,end)-Ja.ss.Rsend_ad_e(find(b.z.w==b.geom.z(n+1))-2);

                                end
                                
                                coup.Ja.zmend.Rs=[coup.Ja.zmend.Rs ss.Rsend];
                                
                            end
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.zm=-[sparse(J-1,(b.K+1)*(b.J+2)) b.eend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J)];'])

                        end

                        eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.Ja.zm=b.' b.bc.z{length(b.bc.z)/2+n} '.zm;'])
                        eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.rhs.zm=sparse(length(b.z.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1))),1);'])

                    else
                        if energy==1
                            coup.Ja.zmend.T(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)=Ja.ss.Tend(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:);
                                
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.st.zm=-[sparse(J-1,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)) b.eend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J + (b.K+2)*(b.J+2))];'])

                % thermal energy
                                coup.th_eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.nhtend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                                coup.Ja.th_eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=Ja.nht.Tend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                                eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.st.th_e=-[sparse(J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J) b.eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])
                                
                        else
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.st.zm=-[sparse(J-1,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)) b.eend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J)];'])

                        end
                        
                    end
                    
            else
                
                % r-momentum
                    % isothermal, no influence of shear stresses on static surface deformation or rigid cylindrical surface
                        if strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'i')==1 || strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'r')==1 || stability==1
                            coup.rmend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.eend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                    % static surface deformation, interior or periodic
                        else
                            coup.rmend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.eend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);

                            if energy==1
                                eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.rm=-[b.nsend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) b.nsend.w(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) b.nsend.p(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+2))];'])

                                eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.Ja.rm=-[b.nsend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) b.nsend.w(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) b.nsend.p(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) Ja.ns.Tend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])

                            else
                                eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.rm=-[b.nsend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) b.nsend.w(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) b.nsend.p(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])

                                eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.Ja.rm=b.' b.bc.z{length(b.bc.z)/2+n} '.rm;'])

                            end
                            
                            ns.Rsend=Ja.ns.Rsend(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                            if n>1
                                ns.Rsend(length(b.z.u(b.z.u<b.geom.z(n))),1)=0;
                                ns.Rsend(length(b.z.u(b.z.u<b.geom.z(n)))+1,1)=ns.Rsend(length(b.z.u(b.z.u<b.geom.z(n)))+1,1)-Ja.ns.Rsend_ad_W(length(b.z.u(b.z.u<b.geom.z(n))));
                                
                            end
                            if n<length(b.bc.z)/2
                                ns.Rsend(length(b.z.u(b.z.u<b.geom.z(n+1)))+1,end)=0;
                                ns.Rsend(length(b.z.u(b.z.u<b.geom.z(n+1))),end)=ns.Rsend(length(b.z.u(b.z.u<b.geom.z(n+1))),end)-Ja.ns.Rsend_ad_E(length(b.z.u(b.z.u<b.geom.z(n+1)))-1);
                                
                            end
                                                        
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.Ja.rmend.Rs=ns.Rsend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);'])
                            
                            % surface (liquid)
                                if strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'s')==1 && isfield(b, 'sigma')==1
                                    eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.rhs.rm=HPend(1,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)).''+LPend(1,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)).'';'])

                            % surface (gas), interior or periodic
                                elseif strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'s')==1 || strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'c')==1 || strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'p')==1
                                    eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.rhs.rm=HPend(1,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1)).'';'])

                                end

                        end

                % phi-momentum
                    if stability==1
                        coup.st.phimend.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.eend.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                    
                        if energy==1
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.st.phim=-[st.assend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) st.assend.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+1)+b.K*b.J) st.assend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])

                        else
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.st.phim=-[st.assend.u(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) st.assend.v(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:) sparse(J,(b.K+2)*(b.J+1)+b.K*b.J)];'])

                        end
                        
                    end
                        
                % z-momentum
                    coup.zmend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)=b.eend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:);

                    if stability~=1
                        if energy==1
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.zm=-[b.ssend.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) b.ssend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J + (b.K+2)*(b.J+2))];'])

                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.Ja.zm=-[b.ssend.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) b.ssend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J) Ja.ss.Tend(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)];'])
                            if strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'s')==1
                                ss.Rsend=Ja.ss.Rsend(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                                if n>1
                                    ss.Rsend(find(b.z.w==b.geom.z(n))-1,1)=0;
                                    ss.Rsend(b.z.w==b.geom.z(n),1:2)=0;
                                    ss.Rsend(find(b.z.w==b.geom.z(n))+1,1)=ss.Rsend(find(b.z.w==b.geom.z(n))+1,1)-Ja.ss.Rsend_ad_w(b.z.w==b.geom.z(n));

                                end
                                if n<length(b.bc.z)/2
                                    ss.Rsend(b.z.w==b.geom.z(n+1),end-1:end)=0;
                                    ss.Rsend(find(b.z.w==b.geom.z(n+1))+1,end)=0;
                                    ss.Rsend(find(b.z.w==b.geom.z(n+1))-1,end)=ss.Rsend(find(b.z.w==b.geom.z(n+1))-1,end)-Ja.ss.Rsend_ad_e(find(b.z.w==b.geom.z(n+1))-2);

                                end
                                
                                eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.Ja.zmend.Rs=ss.Rsend(b.z.w>=b.geom.z(n) & b.z.w<=b.geom.z(n+1),:);'])

                            end
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.rhs.zm=THSend(1,b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1)).'';'])

                % thermal energy
                            coup.th_eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                            coup.Ja.th_eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=coup.th_eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                            
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.th_e=-[sparse(J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J) b.nhtend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])

                            if strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'s')==1
                                nht.Rsend=Ja.nht.Rsend(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                                if n>1
                                    nht.Rsend(b.z.w==b.geom.z(n),1)=0;
                                    nht.Rsend(find(b.z.w==b.geom.z(n))+1,1)=nht.Rsend(find(b.z.w==b.geom.z(n))+1,1)-Ja.nht.Rsend_ad_w(b.z.w==b.geom.z(n));
                                    if  ((strcmp(b.bc.z{n-1}(1),'w')==1 || strcmp(b.bc.z{n-1}(1),'i')==1 || strcmp(b.bc.z{n-1}(1),'o')==1) && strcmp(b.bc.z{n-1}(3),'a')==1)
                                        nht.Rsend(b.z.w==b.geom.z(n),1)=Ja.ad.Rsend_E(b.z.w==b.geom.z(n));

                                    end

                                end
                                if n<length(b.bc.z)/2
                                    nht.Rsend(find(b.z.w==b.geom.z(n+1))+1,end)=0;
                                    nht.Rsend(b.z.w==b.geom.z(n+1),end)=nht.Rsend(b.z.w==b.geom.z(n+1),end)-Ja.nht.Rsend_ad_e(find(b.z.w==b.geom.z(n+1))-1);
                                    if  ((strcmp(b.bc.z{n+1}(1),'w')==1 || strcmp(b.bc.z{n+1}(1),'i')==1 || strcmp(b.bc.z{n+1}(1),'o')==1) && strcmp(b.bc.z{n+1}(3),'a')==1)
                                        nht.Rsend(find(b.z.w==b.geom.z(n+1))+1,end)=Ja.ad.Rsend_W(find(b.z.w==b.geom.z(n+1))+1);

                                    end

                                end

                                eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.Ja.th_eend.Rs=nht.Rsend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);'])

                            end
                            
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.Ja.th_e=-[sparse(J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J) Ja.nht.Tend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])

                        else
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.zm=-[b.ssend.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) b.ssend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J)];'])

                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.Ja.zm=b.' b.bc.z{length(b.bc.z)/2+n} '.zm;'])

                            if strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'s')==1
                                ss.Rsend=Ja.ss.Rsend(:,b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1));
                                if n>1
                                    ss.Rsend(find(b.z.w==b.geom.z(n))-1,1)=0;
                                    ss.Rsend(b.z.w==b.geom.z(n),1:2)=0;
                                    ss.Rsend(find(b.z.w==b.geom.z(n))+1,1)=ss.Rsend(find(b.z.w==b.geom.z(n))+1,1)-Ja.ss.Rsend_ad_w(b.z.w==b.geom.z(n));

                                end
                                if n<length(b.bc.z)/2
                                    ss.Rsend(b.z.w==b.geom.z(n+1),end-1:end)=0;
                                    ss.Rsend(find(b.z.w==b.geom.z(n+1))+1,end)=0;
                                    ss.Rsend(find(b.z.w==b.geom.z(n+1))-1,end)=ss.Rsend(find(b.z.w==b.geom.z(n+1))-1,end)-Ja.ss.Rsend_ad_e(find(b.z.w==b.geom.z(n+1))-2);

                                end
                                
                                eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.Ja.zmend.Rs=ss.Rsend(b.z.w>=b.geom.z(n) & b.z.w<=b.geom.z(n+1),:);'])

                            end 
                            
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.rhs.zm=sparse(length(b.z.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1))),1);'])

                        end
                        
                    else
                        if energy==1
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.st.zm=-[b.ssend.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,(b.K+2)*(b.J+2)) b.ssend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J) Ja.ss.Tend(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:)];'])
                                    
                % thermal energy
                            coup.th_eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=b.eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                            coup.Ja.th_eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)=coup.th_eend.T(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:);
                            
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.st.th_e=-[sparse(J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J) Ja.nht.Tend(b.z.u>b.geom.z(n) & b.z.u<b.geom.z(n+1),:)];'])

                        else
                            eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.st.zm=-[b.ssend.u(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,(b.K+2)*(b.J+2)) b.ssend.w(b.z.w>b.geom.z(n) & b.z.w<b.geom.z(n+1),:) sparse(J-1,b.K*b.J)];'])

                        end
                    end

            end

        end
        
    end
    
    if stability~=1
        if energy==1
            b.coup=[coup.rm1.u   coup.rm1.w    coup.rm1.p    sparse((b.J+2),(b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend);
                   sparse((b.K-1)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) +b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend);
                   coup.rmend.u  coup.rmend.w  coup.rmend.p  sparse((b.J+2),(b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend);
                   coup.zm1.u    coup.zm1.w    sparse((b.J+1),b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend);
                   sparse((b.K)*(b.J+1),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) +b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend);
                   coup.zmend.u  coup.zmend.w  sparse((b.J+1),b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend);
                   sparse((b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J)  coup.th_e1.T  sparse((b.J+2),b.sizeR1 + b.sizeRend);
                   sparse((b.K)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) +b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend);
                   sparse((b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J)  coup.th_eend.T  sparse((b.J+2),b.sizeR1 + b.sizeRend);
                   sparse(b.K*b.J + b.sizeR1 + b.sizeRend, (b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) +b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend)];
            
            b.Ja.coup=b.coup;
            b.Ja.coup(1:(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2))                                                                                                  = coup.Ja.rm1.T;
            b.Ja.coup(b.K*(b.J+2) + 1:(b.K+1)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2))                                                                            = coup.Ja.rmend.T;
            b.Ja.coup((b.K+1)*(b.J+2) + 1:(b.K+1)*(b.J+2) + b.J+1,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2))                                                                = coup.Ja.zm1.T;
            b.Ja.coup((b.K+1)*(b.J+2) + (b.K+1)*(b.J+1) + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2))                                    = coup.Ja.zmend.T;
            b.Ja.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.J+2,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2))                            = coup.Ja.th_e1.T;
            b.Ja.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+1)*(b.J+2) +1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2)) = coup.Ja.th_eend.T;

            if b.sizeR1>0
                b.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1,1:(b.K+1)*(b.J+2))                                                  = b.dbc.R1;
                b.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2) + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1) = b.C.Rs1;

                b.Ja.r_m.Rs1(1:b.J+2,:)=coup.Ja.rm1.Rs;
                b.Ja.z_m.Rs1(1:b.J+1,:)=coup.Ja.zm1.Rs;
                b.Ja.th_e.Rs1(1:b.J+2,:)=coup.Ja.th_e1.Rs;

                b.Ja.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1,1:(b.K+1)*(b.J+2))                                               = b.dbc.R1;
                b.Ja.coup(1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2) + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1) = [b.Ja.r_m.Rs1; b.Ja.z_m.Rs1; b.Ja.th_e.Rs1; b.Ja.C.Rs1];

            end
            if b.sizeRend>0
                b.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J+ b.sizeR1 + b.sizeRend,1:(b.K+1)*(b.J+2))                                                   = b.dbc.Rend;
                b.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend) = b.C.Rsend;

                b.Ja.r_m.Rsend(b.K*(b.J+2)+1:(b.K+1)*(b.J+2),:)      = coup.Ja.rmend.Rs;
                b.Ja.z_m.Rsend((b.K+1)*(b.J+1)+1:(b.K+2)*(b.J+1),:)  = coup.Ja.zmend.Rs;
                b.Ja.th_e.Rsend((b.K+1)*(b.J+2)+1:(b.K+2)*(b.J+2),:) = coup.Ja.th_eend.Rs;

                b.Ja.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J+ b.sizeR1 + b.sizeRend,1:(b.K+1)*(b.J+2))                                                = b.dbc.Rend;
                b.Ja.coup(1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + (b.K+2)*(b.J+2) + b.K*b.J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + (b.K+2)*(b.J+2) + b.sizeR1 + b.sizeRend) = [b.Ja.r_m.Rsend; b.Ja.z_m.Rsend; b.Ja.th_e.Rsend; b.Ja.C.Rsend];

            end

            b.rhs=b.rhs + [rhs.rm1; sparse((b.K-1)*(b.J+2),1); rhs.rmend; rhs.zm1; sparse((b.K)*(b.J+1),1); rhs.zmend; sparse((b.K+2)*(b.J+2),1); rhs.p; sparse(b.sizeR1 + b.sizeRend,1)];
            
        else
            b.coup=[coup.rm1.u   coup.rm1.w    coup.rm1.p    sparse((b.J+2),b.sizeR1 + b.sizeRend);
                   sparse((b.K-1)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) +b.K*b.J + b.sizeR1 + b.sizeRend);
                   coup.rmend.u  coup.rmend.w  coup.rmend.p  sparse((b.J+2),b.sizeR1 + b.sizeRend);
                   coup.zm1.u    coup.zm1.w    sparse((b.J+1),b.K*b.J + b.sizeR1 + b.sizeRend);
                   sparse((b.K)*(b.J+1),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) +b.K*b.J + b.sizeR1 + b.sizeRend);
                   coup.zmend.u  coup.zmend.w  sparse((b.J+1),b.K*b.J + b.sizeR1 + b.sizeRend);
                   sparse(b.K*b.J + b.sizeR1 + b.sizeRend, (b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) +b.K*b.J + b.sizeR1 + b.sizeRend)];

            b.Ja.coup=b.coup;
            
            if b.sizeR1>0
                b.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1,1:(b.K+1)*(b.J+2))                                = b.dbc.R1;
                b.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1) = b.C.Rs1;

                b.Ja.r_m.Rs1(1:b.J+2,:)  = coup.Ja.rm1.Rs;
                b.Ja.z_m.Rs1(1:b.J+1,:)  = coup.Ja.zm1.Rs;
                
                b.Ja.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1,1:(b.K+1)*(b.J+2))                             = b.dbc.R1;
                b.Ja.coup(1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1) = [b.Ja.r_m.Rs1; b.Ja.z_m.Rs1; b.Ja.C.Rs1];

            end
            if b.sizeRend>0
                b.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J+ b.sizeR1 + b.sizeRend,1:(b.K+1)*(b.J+2))                                                   = b.dbc.Rend;
                b.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + b.sizeRend) = b.C.Rsend;

                b.Ja.r_m.Rsend(b.K*(b.J+2)+1:(b.K+1)*(b.J+2),:)      = coup.Ja.rmend.Rs;
                b.Ja.z_m.Rsend((b.K+1)*(b.J+1)+1:(b.K+2)*(b.J+1),:)  = coup.Ja.zmend.Rs;
                
                b.Ja.coup((b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J+ b.sizeR1 + b.sizeRend,1:(b.K+1)*(b.J+2))                                                = b.dbc.Rend;
                b.Ja.coup(1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + 1:(b.K+1)*(b.J+2) + (b.K+2)*(b.J+1) + b.K*b.J + b.sizeR1 + b.sizeRend) = [b.Ja.r_m.Rsend; b.Ja.z_m.Rsend; b.Ja.C.Rsend];

            end

            b.rhs=b.rhs + [rhs.rm1; sparse((b.K-1)*(b.J+2),1); rhs.rmend; rhs.zm1; sparse((b.K)*(b.J+1),1); rhs.zmend; rhs.p; sparse(b.sizeR1 + b.sizeRend,1)];
          
        end
        
    else
        if energy==1
            b.coup=[coup.rm1.u   sparse(b.J+2,(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J+(b.K+2)*(b.J+2));
                   sparse((b.K-1)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)+(b.K+2)*(b.J+1) +b.K*b.J + (b.K+2)*(b.J+2));
                   coup.rmend.u  sparse(b.J+2,(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J+(b.K+2)*(b.J+2));
                   coup.st.phim1.u coup.st.phim1.v sparse(b.J+2,(b.K+2)*(b.J+1)+b.K*b.J) coup.st.phim1.T;
                   sparse((b.K)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)+(b.K+2)*(b.J+1) +b.K*b.J + (b.K+2)*(b.J+2));
                   coup.st.phimend.u coup.st.phimend.v sparse(b.J+2,(b.K+2)*(b.J+1)+b.K*b.J) coup.st.phimend.T;
                   coup.zm1.u sparse(b.J+1,(b.K+2)*(b.J+2)) coup.zm1.w sparse(b.J+1,b.K*b.J) coup.Ja.zm1.T;
                   sparse((b.K)*(b.J+1),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)+(b.K+2)*(b.J+1) +b.K*b.J + (b.K+2)*(b.J+2));
                   coup.zmend.u sparse(b.J+1,(b.K+2)*(b.J+2)) coup.zmend.w sparse(b.J+1,b.K*b.J) coup.Ja.zmend.T;
                   sparse((b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)+(b.K+2)*(b.J+1) +b.K*b.J) coup.Ja.th_e1.T;
                   sparse((b.K)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)+(b.K+2)*(b.J+1) +b.K*b.J + (b.K+2)*(b.J+2));
                   sparse((b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)+(b.K+2)*(b.J+1) +b.K*b.J) coup.Ja.th_eend.T;
                   sparse(b.K*b.J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)+(b.K+2)*(b.J+1) +b.K*b.J + (b.K+2)*(b.J+2))];
            
        else
            b.coup=[coup.rm1.u   sparse(b.J+2,(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J);
                   sparse((b.K-1)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)+(b.K+2)*(b.J+1) +b.K*b.J);
                   coup.rmend.u  sparse(b.J+2,(b.K+2)*(b.J+2)+(b.K+2)*(b.J+1)+b.K*b.J);
                   coup.st.phim1.u coup.st.phim1.v sparse(b.J+2,(b.K+2)*(b.J+1)+b.K*b.J);
                   sparse((b.K)*(b.J+2),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)+(b.K+2)*(b.J+1) +b.K*b.J);
                   coup.st.phimend.u coup.st.phimend.v sparse(b.J+2,(b.K+2)*(b.J+1)+b.K*b.J);
                   coup.zm1.u sparse(b.J+1,(b.K+2)*(b.J+2)) coup.zm1.w sparse(b.J+1,b.K*b.J);
                   sparse((b.K)*(b.J+1),(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)+(b.K+2)*(b.J+1) +b.K*b.J);
                   coup.zmend.u sparse(b.J+1,(b.K+2)*(b.J+2)) coup.zmend.w sparse(b.J+1,b.K*b.J);
                   sparse(b.K*b.J,(b.K+1)*(b.J+2) + (b.K+2)*(b.J+2)+(b.K+2)*(b.J+1) +b.K*b.J)];
                        
        end        
        
    end