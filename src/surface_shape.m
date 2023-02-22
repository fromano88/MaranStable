function[r_lb_l, dp, convergence] = surface_shape(r_lb_l, dp, V_r, g, rho_l, rho_g, sigma, r_c, r_i, l_lb, l_d1, delta_lb)

V_l=V_r*(r_i^2-r_c^2)*pi*l_lb;

%nf=1;
J_s=round(l_lb/delta_lb*2);%11*2^nf-2^(nf-1)+1;
dz=ones(1,J_s)*l_lb/J_s;
z=(0:J_s)*l_lb/J_s;

if isempty(r_lb_l)==0
    % disp('Solution vector from previous calculations!!!') 
    r=r_lb_l(z+l_d1);
    rr=[r,dp].';
    
else
    r=zeros(1,J_s+1);
    r(1,1:J_s+1)=sign(r_i)*sqrt(V_l/pi/l_lb);
    dp=sign(r_i)*sigma/r_i;
    rr=[r,dp].';
    
end
tol=10^-6;

res=1;    

while res>tol && res < 1/tol
    r1=zeros(1,J_s+1);
    r2=zeros(1,J_s+1);
    r3=zeros(1,J_s+1);
    delta_p=zeros(1,J_s+1);
    vol1=zeros(1,J_s+1);
    vol2=zeros(1,J_s+1);
    
    r2(1,1)=1;
    r2(1,J_s+1)=1;


    r1(1,1:J_s-1)=r(1,2:J_s)./dz(1,2:J_s).^2+(r(1,3:J_s+1)-r(1,1:J_s-1))./(4*dz(1,2:J_s).^2);
    r2(1,2:J_s)=-2*r(1,2:J_s)./dz(1,2:J_s).^2-(rho_l-rho_g)*g/sigma.*z(1,2:J_s).*(1+((r(1,3:J_s+1)-r(1,1:J_s-1))./(2*dz(1,2:J_s))).^2).^(3/2);
    r3(1,3:J_s+1)=r(1,2:J_s)./dz(1,2:J_s).^2-(r(1,3:J_s+1)-r(1,1:J_s-1))./(4*dz(1,2:J_s).^2);
    rad=spdiags([r1.' r2.' r3.'],[-1 0 1],J_s+1,J_s+1);

    delta_p(1,2:J_s)=1/sigma*r(1,2:J_s).*(1+((r(1,3:J_s+1)-r(1,1:J_s-1))./(2*dz(1,2:J_s))).^2).^(3/2);
    delta_p=delta_p.';

    vol1(1,1:J_s)=(r(1,1:J_s)+r(1,2:J_s+1))/4.*dz(1,1:J_s);
    vol2(1,2:J_s+1)  =(r(1,1:J_s)+r(1,2:J_s+1))/4.*dz(1,1:J_s);
    vol=[vol1+vol2,0];
    
    A=[[rad delta_p];vol];
    
    b=[r_i;ones(J_s-1,1);r_i;V_l/pi+r_c^2*l_lb];

    Jr1=zeros(1,J_s+1);
    Jr2=zeros(1,J_s+1);
    Jr3=zeros(1,J_s+1);
    Jvol1=zeros(1,J_s+1);
    Jvol2=zeros(1,J_s+1);
    
    Jr2(1,1)=1;
    Jr2(1,J_s+1)=1;

    Jr1(1,1:J_s-1)=r(1,2:J_s)./dz(1,2:J_s).^2+(r(1,3:J_s+1)-r(1,1:J_s-1))./(2*dz(1,2:J_s).^2)-(dp/sigma-(rho_l-rho_g)*g/sigma*z(1,2:J_s))*3/2.*r(1,2:J_s).*(1+((r(1,3:J_s+1)-r(1,1:J_s-1))./(2*dz(1,2:J_s))).^2).^0.5.*(r(1,3:J_s+1)-r(1,1:J_s-1))./(2*dz(1,2:J_s).^2);
    Jr2(1,2:J_s)=(r(1,3:J_s+1)-4*r(1,2:J_s)+r(1,1:J_s-1))./dz(1,2:J_s).^2+(dp/sigma-(rho_l-rho_g)*g/sigma*z(1,2:J_s)).*(1+((r(1,3:J_s+1)-r(1,1:J_s-1))./(2*dz(1,2:J_s))).^2).^(3/2);
    Jr3(1,3:J_s+1)=r(1,2:J_s)./dz(1,2:J_s).^2-(r(1,3:J_s+1)-r(1,1:J_s-1))./(2*dz(1,2:J_s).^2)+(dp/sigma-(rho_l-rho_g)*g/sigma*z(1,2:J_s))*3/2.*r(1,2:J_s).*(1+((r(1,3:J_s+1)-r(1,1:J_s-1))./(2*dz(1,2:J_s))).^2).^0.5.*(r(1,3:J_s+1)-r(1,1:J_s-1))./(2*dz(1,2:J_s).^2);
    Jrad=spdiags([Jr1.' Jr2.' Jr3.'],[-1 0 1],J_s+1,J_s+1);

    Jvol1(1,1:J_s)=(r(1,1:J_s)+r(1,2:J_s+1)).*dz(1,1:J_s)/2;
    Jvol2(1,2:J_s+1)  =(r(1,1:J_s)+r(1,2:J_s+1)).*dz(1,1:J_s)/2;
    Jvol=[Jvol1+Jvol2,0];
    
    J=[[Jrad delta_p];Jvol];
    
    warning off
    drr=-J\(A*rr-b);
    warning on
    rr=rr+drr;
    %rr=(A\b);
    res=norm(A*rr-b,inf);
    dp=rr(J_s+2,1);
    r=rr(1:J_s+1,1).';
end

% display(['res = ', num2str(res)])

%mean_curvature=1./r(1,2:J_s)./(1+((r(1,3:J_s+1)-r(1,1:J_s-1))./2./dz(1,1:J_s-1)).^2).^0.5 - (r(1,3:J_s+1)-2.*r(1,2:J_s)+r(1,1:J_s-1))./dz(1,1:J_s-1).^2./(1+((r(1,3:J_s+1)-r(1,1:J_s-1))./2./dz(1,1:J_s-1)).^2).^1.5;
%Bond1=[z(1,1:J_s+1).'-0.5 r(1,1:J_s+1).' [1./r(1,1)./(1+((r(1,2)-r(1,1))./dz(1,1)).^2).^0.5; mean_curvature'; 1./r(1,J_s+1)./(1+((r(1,J_s+1)-r(1,J_s))./dz(1,J_s)).^2).^0.5]];
%{
figcount = 1;
plot(z+l_d1,r)
axis([l_d1 l_lb+l_d1 0 max(r)*1.05])
daspect([1 1 1])
%}
r_lb_l=spline(z+l_d1,r);
%r_lb_l=spapi(3,z+l_d1,r);
%dr_lb_l= fnder(r_lb_l);
r_lb_l=@(z_lb) fnval(r_lb_l,z_lb);
%dr_lb_l=@(z_lb) fnval(dr_lb_l,z_lb);
%r_lb_l=@(z_lb) interp1(z,r,z_lb-l_d1,'spline');
%pp=cscvn([z;r]);
%[breaks,coefs,l,k,d] = unmkpp(pp);
%ppval(pp,breaks(end));
%dr_lb_l=@(z_lb) diff(interp1(z,r,z_lb-l_d1,'spline'),z_lb);

% checking if while loop ended because solution converged or because endless loop was prevented
if res > 1/tol
    convergence = 0;
else
    convergence = 1;
end