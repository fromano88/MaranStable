if o==1
    % lower boundary
        index_bc_z=0;
        index_bc_r=1;
        K=1;
        R_z.w=b.R1_z.w;

elseif o==2
    % upper boundary
        index_bc_z=length(b.bc.z)/2;
        index_bc_r=length(b.bc.r)/2;
        K=b.K+2;
        R_z.w=b.Rend_z.w;

end

if max(strncmp('s', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},1))
    N.w=sqrt(1+R_z.w(1,:).^2);

    % axial shear stresses
        F.w=-mu(b.T.w(K,2:end-1)).*b.XI_z.w(K,2:end-1).*b.dz.zm(2:end-1).*N.w(1,2:end-1);

        u_sw=zeros(1,b.J+1); u_se=u_sw; u_nw=u_sw; u_ne=u_sw;

        u_sw(1,2:end-1) = (-1)^o*F.w./b.XI_r.v(K-(o-1),2:end-1)./b.DXI.zm(1,2:end-1)/2 - mu(b.T.w(K,2:end-1))./N.w(1,2:end-1).*(1-R_z.w(1,2:end-1).^2).*b.ETA_z.w(K,2:end-1)./b.XI_r.u(K-(o-1),2:end-2)./b.DETA.zm(1,2:end-1).*b.dz.zm(2:end-1);
        u_se(1,2:end-1) = (-1)^o*F.w./b.XI_r.v(K-(o-1),2:end-1)./b.DXI.zm(1,2:end-1)/2 + mu(b.T.w(K,2:end-1))./N.w(1,2:end-1).*(1-R_z.w(1,2:end-1).^2).*b.ETA_z.w(K,2:end-1)./b.XI_r.u(K-(o-1),3:end-1)./b.DETA.zm(1,2:end-1).*b.dz.zm(2:end-1);
        u_nw(1,2:end-1) = - (-1)^o*F.w./b.XI_r.v(K-o*(-1)^o,2:end-1)./b.DXI.zm(1,2:end-1)/2;
        u_ne(1,2:end-1) = - (-1)^o*F.w./b.XI_r.v(K-o*(-1)^o,2:end-1)./b.DXI.zm(1,2:end-1)/2;

        F.w=mu(b.T.w(K,2:end-1)).*b.XI_z.w(K,2:end-1).*b.dz.zm(2:end-1).*N.w(1,2:end-1);
        F2.w=-mu(b.T.w(K,2:end-1))./N.w(1,2:end-1).*(1-R_z.w(1,2:end-1).^2).*b.ETA_z.w(K,2:end-1).*b.dz.zm(2:end-1);
        F4.w=mu(b.T.w(K,2:end-1)).*b.XI_r.w(K,2:end-1).*b.dz.zm(2:end-1).*N.w(1,2:end-1);
        F5.w=-mu(b.T.w(K,2:end-1))./N.w(1,2:end-1).*2.*R_z.w(1,2:end-1).*b.ETA_z.w(K,2:end-1).*b.dz.zm(2:end-1);

        w_Wsw=zeros(1,b.J+1); w_s=w_Wsw; w_Ese=w_Wsw; w_n=w_Wsw;

        w_Wsw(1,2:end-1) = - F2.w.*b.JA.u(K-(o-1),2:end-2).*b.XI_z.u(K-(o-1),2:end-2)./b.DETA.zm(1,2:end-1)/2 - F5.w./b.ETA_z.u(K-(o-1),2:end-2)./b.DETA.zm(1,2:end-1)/2;
        w_Ese(1,2:end-1) = + F2.w.*b.JA.u(K-(o-1),3:end-1).*b.XI_z.u(K-(o-1),3:end-1)./b.DETA.zm(1,2:end-1)/2 + F5.w./b.ETA_z.u(K-(o-1),3:end-1)./b.DETA.zm(1,2:end-1)/2;
        w_s(1,2:end-1)   = w_Wsw(1,2:end-1) + w_Ese(1,2:end-1) + (-1)^o*(F.w.*b.JA.w(K,2:end-1).*b.XI_z.w(K,2:end-1)./b.DXI.zm(1,2:end-1)*2 + F4.w./b.ETA_z.w(K,2:end-1)./b.DXI.zm(1,2:end-1)*2);
        w_n(1,2:end-1)   = - (-1)^o*(F.w.*b.JA.w(K-(-1)^o,2:end-1).*b.XI_z.w(K-(-1)^o,2:end-1)./b.DXI.zm(1,2:end-1)*2 +F4.w./b.ETA_z.w(K-(-1)^o,2:end-1)./b.DXI.zm(1,2:end-1)*2);

    % equal velocity
        equal_w=N.w./b.ETA_z.w(K,:);
    
    % thermocapillary stresses
        if energy==1
            if thermcapcon==1 && isfield(b, 'sigma')==1
                if size(T_fluid)==0
                    THS=-(-1)^o*b.ETA_z.w(K,:).*(sigma(b.T.u(K-(o-1),2:end))-sigma(b.T.u(K-(o-1),1:end-1)))./b.DETA.zm(1,:).*b.dz.zm;

                else
                    THS=-(-1)^o*b.ETA_z.w(K,:).*dsigma(b.T.w(K,:)).*(b.T.u(K-(o-1),2:end)-b.T.u(K-(o-1),1:end-1))./b.DETA.zm(1,:).*b.dz.zm;

                end
                if isfield(flowopt,'regularization') == 1
                    THS = THS.*regularizationFunction(flowopt.regularization.function,b.J+1,flowopt.regularization.percentage);
                end
            else
                THS=sparse(1,b.J+2);

            end

        % derivative of axial shear stresses with respect to T
            % -dmu/dT*xi_z*d(u0/xi_r)/delta.xi
                FC.w=-dmu(b.T.w(K,:)).*b.XI_z.w(K,:).*(b.u.v(K-(-1)^o,:)./b.XI_r.v(K-(-1)^o,:)-b.u.v(K-(o-1)*2,:)./b.XI_r.v(K-(o-1)*2,:))./b.DXI.zm(1,:);

            % dmu/dT*xi_z*d(J*xi_z*w0)/delta.xi
                FC.w=FC.w+dmu(b.T.w(K,:)).*b.XI_z.w(K,:).*(b.JA.w(K+(2-o),:).*b.XI_z.w(K+(2-o),:).*b.w.w(K+(2-o),:)-b.JA.w(K+(1-o),:).*b.XI_z.w(K+(1-o),:).*b.w.w(K+(1-o),:))./b.DXI.zm(1,:)*2;

            % dmu/dT*xi_r*d(w0/eta_z)/delta.xi
                FC.w=FC.w+dmu(b.T.w(K,:)).*b.XI_r.w(K,:).*(b.w.w(K+(2-o),:)./b.ETA_z.w(K+(2-o),:)-b.w.w(K+(1-o),:)./b.ETA_z.w(K+(1-o),:))./b.DXI.zm(1,:)*2;

            % dmu/dT/N^2*(1-r_z^2)*eta_z*d(u0/xi_r)/delta.eta
                FC.w=FC.w+dmu(b.T.w(K,:))./N.w.^2.*(1-R_z.w(1,:).^2).*b.ETA_z.w(K,:).*(b.u.u(K-(o-1),2:end)./b.XI_r.u(K-(o-1),2:end)-b.u.u(K-(o-1),1:end-1)./b.XI_r.u(K-(o-1),1:end-1))./b.DETA.zm(1,:);

            % -dmu/dT/N^2*(1-r_z^2)*eta_z*d(J*xi_z*w0)/delta.eta
                FC.w=FC.w-dmu(b.T.w(K,:))./N.w.^2.*(1-R_z.w(1,:).^2).*b.ETA_z.w(K,:).*(b.JA.u(K-(o-1),2:end).*b.XI_z.u(K-(o-1),2:end).*b.w.u(K-(o-1),2:end)-b.JA.u(K-(o-1),1:end-1).*b.XI_z.u(K-(o-1),1:end-1).*b.w.u(K-(o-1),1:end-1))./b.DETA.zm(1,:);

            % -dmu/dT/N^2*2*r_z*eta_z*d(w0/ETA_z)/delta.eta
                FC.w=FC.w-dmu(b.T.w(K,:))./N.w.^2.*2.*R_z.w(1,:).*b.ETA_z.w(K,:).*(b.w.u(K-(o-1),2:end)./b.ETA_z.u(K-(o-1),2:end)-b.w.u(K-(o-1),1:end-1)./b.ETA_z.u(K-(o-1),1:end-1))./b.DETA.zm(1,:);

            % derivative of thermocapillary stresses with respect to T
                if thermcapcon==1 && isfield(b, 'sigma')==1
                    F_w.w=-(-1)^o*1./N.w.*b.ETA_z.w(K,:).*(dsigma(b.T.u(K-(o-1),1:end-1)))./b.DETA.zm(1,:);
                    F_e.w=-(-1)^o*1./N.w.*b.ETA_z.w(K,:).*(dsigma(b.T.u(K-(o-1),2:end)))./b.DETA.zm(1,:);
                    
                else
                    F_w.w=sparse(1,b.J+1);
                    F_e.w=sparse(1,b.J+1);
                    
                end

            T_w=(b.W.w(1,2:end).*FC.w(2:end-1) - F_w.w(2:end-1)).*b.dz.zm(2:end-1).*N.w(1,2:end-1);
            T_e=(b.E.w(1,1:end-1).*FC.w(2:end-1) + F_e.w(2:end-1)).*b.dz.zm(2:end-1).*N.w(1,2:end-1);        

        else
            THS=sparse(1,b.J+2);

        end

   if stability~=1 && max(strncmp('ss', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},2))
        for p=1:2
            if max(strncmp('ss', {b.bc.z{(p-1)*length(b.bc.z)/2+1:(p-1)*length(b.bc.z)/2+length(b.bc.z)/2}},2))
                if p==1
                    % derivative with respect to lower boundary
                        R_z_rlb.w=b.R1_z_rlb1.w;
                        XI_r_rlb.u=b.XI_r_rlb1.u;
                        XI_r_rlb.w=b.XI_r_rlb1.w;
                        XI_r_inv_rlb.u=b.XI_r_inv_rlb1.u;
                        XI_r_inv_rlb.v=b.XI_r_inv_rlb1.v;
                        XI_z_rlb_C.u=b.XI_z_rlb1_C.u;
                        XI_z_rlb_C.w=b.XI_z_rlb1_C.w;
                        XI_z_rlb_E.u=b.XI_z_rlb1_E.u;
                        XI_z_rlb_E.w=b.XI_z_rlb1_E.w;
                        JA_rlb.u=b.JA_rlb1.u;
                        JA_rlb.w=b.JA_rlb1.w;
                        JA_rlb.v=b.JA_rlb1.v;

                elseif p==2
                    % derivative with respect to upper boundary
                        R_z_rlb.w=b.Rend_z_rlbend.w;
                        XI_r_rlb.u=b.XI_r_rlbend.u;
                        XI_r_rlb.w=b.XI_r_rlbend.w;
                        XI_r_inv_rlb.u=b.XI_r_inv_rlbend.u;
                        XI_r_inv_rlb.v=b.XI_r_inv_rlbend.v;
                        XI_z_rlb_C.u=b.XI_z_rlbend_C.u;
                        XI_z_rlb_C.w=b.XI_z_rlbend_C.w;
                        XI_z_rlb_E.u=b.XI_z_rlbend_E.u;
                        XI_z_rlb_E.w=b.XI_z_rlbend_E.w;
                        JA_rlb.u=b.JA_rlbend.u;
                        JA_rlb.w=b.JA_rlbend.w;
                        JA_rlb.v=b.JA_rlbend.v;

                end

                % derivative of axial shear stresses with respect to Rs
                    if o==p
                        N_inv_sq_rlb.w=-1./N.w.^4.*2.*R_z.w(1,:).*R_z_rlb.w(1,:);
                        N_inv_rlb.w=-0.5./N.w.^3.*2.*R_z.w(1,:).*R_z_rlb.w(1,:);
                        N_rlb.w=0.5./N.w.*2.*R_z.w(1,:).*R_z_rlb.w(1,:);
                        
                        % -mu*xi_z*N*d(u0/xi_r)/delta.xi*delta.z
                            % -mu*d(xi_z*N)/dRs*d(u0/xi_r)/delta.xi*delta.z
                                F.w=(b.u.v(K-(-1)^o,:)./b.XI_r.v(K-(-1)^o,:)-b.u.v(K-(o-1)*2,:)./b.XI_r.v(K-(o-1)*2,:))./b.DXI.zm(1,:).*b.dz.zm;
                                FC.w=-mu(b.T.w(K,:)).*XI_z_rlb_C.w(K,:).*N.w.*F.w;
                                FE.w=-mu(b.T.w(K,:)).*(XI_z_rlb_E.w(K,:).*N.w+N_rlb.w.*b.XI_z.w(K,:)).*F.w;

                            % -mu*xi_z*N*d(u0*(1/xi_r)/dRs)/delta.xi*delta.z
                                FC.w=FC.w-mu(b.T.w(K,:)).*b.XI_z.w(K,:).*N.w.*(b.u.v(K-(-1)^o,:).*XI_r_inv_rlb.v(K-(-1)^o,:)-b.u.v(K-(o-1)*2,:).*XI_r_inv_rlb.v(K-(o-1)*2,:))./b.DXI.zm(1,:).*b.dz.zm;

                        % mu*xi_z*N*d(J*xi_z*w0)/delta.xi*delta.z
                            % mu*d(xi_z*N)/dRs*d(J*xi_z*w0)/delta.xi*delta.z
                                F.w=(b.JA.w(K+(2-o),:).*b.XI_z.w(K+(2-o),:).*b.w.w(K+(2-o),:)-b.JA.w(K+(1-o),:).*b.XI_z.w(K+(1-o),:).*b.w.w(K+(1-o),:))./b.DXI.zm(1,:)*2.*b.dz.zm;
                                FC.w=FC.w+mu(b.T.w(K,:)).*XI_z_rlb_C.w(K,:).*N.w.*F.w;
                                FE.w=FE.w+mu(b.T.w(K,:)).*(XI_z_rlb_E.w(K,:).*N.w+N_rlb.w.*b.XI_z.w(K,:)).*F.w;

                            % mu*xi_z*N*d((dJ/dRs*xi_z+J*dxi_z/dRs)*w0)/delta.xi*delta.z
                                FC.w=FC.w+mu(b.T.w(K,:)).*b.XI_z.w(K,:).*N.w.*((JA_rlb.w(K+(2-o),:).*b.XI_z.w(K+(2-o),:)+b.JA.w(K+(2-o),:).*XI_z_rlb_C.w(K+(2-o),:)).*b.w.w(K+(2-o),:)-(JA_rlb.w(K+(1-o),:).*b.XI_z.w(K+(1-o),:)+b.JA.w(K+(1-o),:).*XI_z_rlb_C.w(K+(1-o),:)).*b.w.w(K+(1-o),:))./b.DXI.zm(1,:)*2.*b.dz.zm;

                            % mu*xi_z*N*d((J*dxi_z/dRs)*w0)/delta.xi*delta.z
                                FE.w=FE.w+mu(b.T.w(K,:)).*b.XI_z.w(K,:).*N.w.*(b.JA.w(K+(2-o),:).*XI_z_rlb_E.w(K+(2-o),:).*b.w.w(K+(2-o),:)-b.JA.w(K+(1-o),:).*XI_z_rlb_E.w(K+(1-o),:).*b.w.w(K+(1-o),:))./b.DXI.zm(1,:)*2.*b.dz.zm;

                        % mu*xi_r*N*d(w0/eta_z)/delta.xi*delta.z
                            % mu*d(xi_r*N)/dRs*d(w0/eta_z)/delta.xi*delta.z
                                F.w=(b.w.w(K+(2-o),:)./b.ETA_z.w(K+(2-o),:)-b.w.w(K+(1-o),:)./b.ETA_z.w(K+(1-o),:))./b.DXI.zm(1,:)*2.*b.dz.zm;
                                FC.w=FC.w+mu(b.T.w(K,:)).*XI_r_rlb.w(K,:).*N.w.*F.w;
                                FE.w=FE.w+mu(b.T.w(K,:)).*b.XI_r.w(K,:).*N_rlb.w.*F.w;

                        % mu/N*(1-r_z^2)*eta_z*d(u0/xi_r)/delta.eta*delta.z
                            % mu/N*(1-r_z^2)*eta_z*d(u0*(1/xi_r)/dRs)/delta.eta*delta.z
                                F.w=mu(b.T.w(K,:))./N.w.*(1-R_z.w(1,:).^2).*b.ETA_z.w(K,:).*b.dz.zm;
                                F_W.w=F.w.*b.u.u(K-(o-1),1:end-1).*XI_r_inv_rlb.u(K-(o-1),1:end-1)./b.DETA.zm(1,:);
                                F_E.w=F.w.*b.u.u(K-(o-1),2:end).*XI_r_inv_rlb.u(K-(o-1),2:end)./b.DETA.zm(1,:);

                            % mu/N*(-2*r_z*dr_z/dRs)*eta_z*d(u0/xi_r)/delta.eta*delta.z
                                F.w=(b.u.u(K-(o-1),2:end)./b.XI_r.u(K-(o-1),2:end) - b.u.u(K-(o-1),1:end-1)./b.XI_r.u(K-(o-1),1:end-1))./b.DETA.zm(1,:).*b.dz.zm;
                                FE.w=FE.w+mu(b.T.w(K,:))./N.w.*(-2*R_z.w(1,:).*R_z_rlb.w(1,:)).*b.ETA_z.w(K,:).*F.w;

                            % mu*d(1/N)/dRs*(1-r_z^2)*eta_z*d(u0/xi_r)/delta.eta*delta.z
                                FE.w=FE.w+mu(b.T.w(K,:)).*N_inv_rlb.w.*(1-R_z.w(1,:).^2).*b.ETA_z.w(K,:).*F.w;

                        % -mu/N*(1-r_z^2)*eta_z*d(J*xi_z*w0)/delta.eta*delta.z
                            % -mu/N*(1-r_z^2)*eta_z*d((dJ/dRs*xi_z+J*dxi_z/dRs)*w0)/delta.eta*delta.z
                                F.w=-mu(b.T.w(K,:))./N.w.*(1-R_z.w(1,:).^2).*b.ETA_z.w(K,:).*b.dz.zm;
                                F_W.w=F_W.w+F.w.*((JA_rlb.u(K-(o-1),1:end-1).*b.XI_z.u(K-(o-1),1:end-1)+b.JA.u(K-(o-1),1:end-1).*XI_z_rlb_C.u(K-(o-1),1:end-1)).*b.w.u(K-(o-1),1:end-1))./b.DETA.zm(1,:);
                                F_E.w=F_E.w+F.w.*((JA_rlb.u(K-(o-1),2:end).*b.XI_z.u(K-(o-1),2:end)+b.JA.u(K-(o-1),2:end).*XI_z_rlb_C.u(K-(o-1),2:end)).*b.w.u(K-(o-1),2:end))./b.DETA.zm(1,:);

                            % -mu/N*(1-r_z^2)*eta_z*d(J*dxi_z/dRs*w0)/delta.eta*delta.z
                                FE_W.w=F.w.*(b.JA.u(K-(o-1),1:end-1).*XI_z_rlb_E.u(K-(o-1),1:end-1).*b.w.u(K-(o-1),1:end-1))./b.DETA.zm(1,:);
                                FE_E.w=F.w.*(b.JA.u(K-(o-1),2:end).*XI_z_rlb_E.u(K-(o-1),2:end).*b.w.u(K-(o-1),2:end))./b.DETA.zm(1,:);

                            % -mu/N*(-2*r_z*dr_z/dRs)*eta_z*d(J*xi_z*w0)/delta.eta*delta.z
                                F.w=(b.JA.u(K-(o-1),2:end).*b.XI_z.u(K-(o-1),2:end).*b.w.u(K-(o-1),2:end) - b.JA.u(K-(o-1),1:end-1).*b.XI_z.u(K-(o-1),1:end-1).*b.w.u(K-(o-1),1:end-1))./b.DETA.zm(1,:).*b.dz.zm;
                                FE.w=FE.w-mu(b.T.w(K,:))./N.w.*(-2*R_z.w(1,:).*R_z_rlb.w(1,:)).*b.ETA_z.w(K,:).*F.w;

                            % -mu*d(1/N)/dRs*(1-r_z^2)*eta_z*d(J*xi_z*w0)/delta.eta*delta.z
                                FE.w=FE.w-mu(b.T.w(K,:)).*N_inv_rlb.w.*(1-R_z.w(1,:).^2).*b.ETA_z.w(K,:).*F.w;

                       % -mu/N^2*2*r_z*eta_z*d(w0/eta_z)/delta.eta*delta.z
                            % -mu/N*2*dr_z/dRs*eta_z*d(w0/eta_z)/delta.eta*delta.z
                                F.w=(b.w.u(K-(o-1),2:end)./b.ETA_z.u(K-(o-1),2:end) - b.w.u(K-(o-1),1:end-1)./b.ETA_z.u(K-(o-1),1:end-1))./b.DETA.zm(1,:).*b.dz.zm;
                                FE.w=FE.w-mu(b.T.w(K,:))./N.w.*2.*R_z_rlb.w(1,:).*b.ETA_z.w(K,:).*F.w;

                            % -mu*d(1/N)/dRs*2*r_z*eta_z*d(J*xi_r*w0)/delta.eta*delta.z
                                FE.w=FE.w-mu(b.T.w(K,:)).*N_inv_rlb.w.*2.*R_z.w(1,:).*b.ETA_z.w(K,:).*F.w;

                        % thermocapillary stresses
                            %{
                            if energy==1 && thermcapcon==1 && isfield(b, 'sigma')==1
                                if size(T_fluid)==0
                                    FE.w=FE.w-(-1)^o*N_inv_rlb.w.*b.ETA_z.w(K,:).*(sigma(b.T.u(K-(o-1),2:end))-sigma(b.T.u(K-(o-1),1:end-1)))./b.DETA.zm(1,:).*b.dz.zm.*N.w;

                                else
                                    FE.w=FE.w-(-1)^o*N_inv_rlb.w.*b.ETA_z.w(K,:).*dsigma(b.T.w(K,:)).*(b.T.u(K-(o-1),2:end)-b.T.u(K-(o-1),1:end-1))./b.DETA.zm(1,:).*b.dz.zm.*N.w;

                                end

                            end
                            %}

                    else
                        FC.w=sparse(1,b.J+1);
                        FE.w=sparse(1,b.J+1);
                        F_W.w=sparse(1,b.J+1);
                        F_E.w=sparse(1,b.J+1);
                        FE_W.w=sparse(1,b.J+1);
                        FE_E.w=sparse(1,b.J+1);

                    end

                % Jacobi of axial shear stresses with respect to Rs
                    Rs_ad_w=[0 b.E.w(1,1:end-2).*FE_W.w(1,3:end-1)];
                    Rs_ad_e=[b.W.w(1,3:end).*FE_E.w(1,2:end-2) 0];
                    Rs_w=b.W.w(1,2:end).*(FC.w(1,2:end-1) - FE_E.w(1,2:end-1) - FE_W.w(1,2:end-1)) - FE.w(1,2:end-1) - F_W.w(1,2:end-1) + Rs_ad_w;
                    Rs_e=b.E.w(1,1:end-1).*(FC.w(1,2:end-1) - FE_E.w(1,2:end-1) - FE_W.w(1,2:end-1)) + FE.w(1,2:end-1) + F_E.w(1,2:end-1) + Rs_ad_e;
                    Rs_W=[0 b.W.w(1,2:end-1).*FE_W.w(1,3:end-1)];
                    Rs_E=[b.E.w(1,2:end-1).*FE_E.w(1,2:end-2) 0];

                    Ja.ss.Rs=sparse(b.J+1,b.J+2);
                    Ja.ss.Rs(2:end-1,:)=spdiags([Rs_W.' Rs_w.' Rs_e.' Rs_E.'],[0 1 2 3],b.J-1,b.J+2);

                if o==1 && o==p
                    % Jacobi of axial shear stresses with respect to Rs
                        Ja.ss.Rs1=Ja.ss.Rs;
                        Ja.ss.Rs1_ad_w=Rs_ad_w;
                        Ja.ss.Rs1_ad_e=Rs_ad_e;

                elseif o==2 && o==p
                    % Jacobi of axial shear stresses with respect to Rs
                        Ja.ss.Rsend=Ja.ss.Rs;
                        Ja.ss.Rsend_ad_w=Rs_ad_w;
                        Ja.ss.Rsend_ad_e=Rs_ad_e;

                end

            end

        end

   end

end
if o==1 && max(strncmp('s', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},1))
    b.ss1.u=spdiags([u_sw.' u_se.' u_nw.' u_ne.'],[0 1 b.J+2 b.J+3],(b.J+1),(b.K+1)*(b.J+2));
    b.ss1.w=spdiags([w_Wsw.' w_s.' w_Ese.' w_n.'],[-1 0 1 b.J+1],(b.J+1),(b.K+2)*(b.J+1));
    if energy==1
        THS1=THS;
    % Jacobi of axial shear stresses with respect to T
        Ja.ss.T1=sparse(b.J+1,(b.K+2)*(b.J+2));
        Ja.ss.T1(2:end-1,:)=spdiags([T_w.' T_e.'],[1 2],b.J-1,(b.K+2)*(b.J+2));

    end
    % equal velocity
        b.e1.w=spdiags(equal_w.',0,b.J+1,(b.K+2)*(b.J+1));

elseif o==2 && max(strncmp('s', {b.bc.z{index_bc_z+1:index_bc_z+length(b.bc.z)/2}},1))
    b.ssend.u=spdiags([u_nw.' u_ne.' u_sw.' u_se.'],(b.K-1)*(b.J+2)+[0 1 b.J+2 b.J+3],(b.J+1),(b.K+1)*(b.J+2));
    b.ssend.w=spdiags([w_n.' w_Wsw.' w_s.' w_Ese.'],(b.K)*(b.J+1)+[0 b.J b.J+1 b.J+2],(b.J+1),(b.K+2)*(b.J+1));
    if energy==1
        THSend=THS;
    % Jacobi of axial shear stresses with respect to T
        Ja.ss.Tend=sparse(b.J+1,(b.K+2)*(b.J+2));
        Ja.ss.Tend(2:end-1,:)=spdiags([T_w.' T_e.'],(b.K)*(b.J+2)+[b.J+3 b.J+4],b.J-1,(b.K+2)*(b.J+2));

    end
    % equal velocity
        b.eend.w=[sparse(b.J+1, (b.K+1)*(b.J+1)) spdiags(equal_w.',0,b.J+1,b.J+1)];

end