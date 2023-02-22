function ray = optical_ray(b,r_0,z_p,N_coeff,alpha_0,reflection)
% Integrate with a Dormand-Prince method (RK4-RK5 schemes)
% Note: n and grad(n) depend on space, not on time
% Note: the grid can also be deformed, the deformation of the liquid bridge is handled by
%       the mapping back and forth from a square prototypical domain. This strategy 
%       allows to interpolate on a tensor-product grid. This is much faster than 
%       interpolanting by scattered interpolants such as griddata. 
%       -> griddata requires a searching algorithm on Nr*Nz nodes
%       -> interp2 is based on a tensor product grid, so, it searches over Nr+Nz nodes

% Jiang He, Wei Liu, Yao-Xiong Huang: Simultaneous Determination of Glass, Transition Temperatures of Several Polymers (p.6)
% Attention: there is a mistake in the paper: replace '+' with '-' (see figure 4 of He16)
% N_coeff = [1.4094 3.9*1E-4 20]; --> N(T) = N_coeff(1)-N_coeff(2)*(T-N_coeff(3))

if nargin < 2
    error('Not enough input parameter.')
    
elseif nargin == 2
    z_p = 0; N_coeff = [1.4094 3.9*1E-4 20]; alpha_0 = 0; reflection = 0;
    
elseif nargin == 3
    N_coeff = [1.4094 3.9*1E-4 20]; alpha_0 = 0; reflection = 0;
    
elseif nargin == 4
    alpha_0 = 0; reflection = 0;

elseif nargin == 5
    reflection = 0;
    
end

% Initial conditions of the ray
t_max   = 5.0;            % maximum intergration time for the light-ray tracing
alpha_0 = alpha_0*pi/180; % input of alpha_0 in °
z_in    = 1.0;            % axial  position of the ray (starts at top rod)
T_d1l   = b.T.T(1,1);     % Temperature of the top wall
u_in    = (N_coeff(1) - N_coeff(2)*(T_d1l-N_coeff(3)))*round(sin(alpha_0),14);   % radial velocity of the ray
w_in    = -(N_coeff(1) - N_coeff(2)*(T_d1l-N_coeff(3)))*round(cos(alpha_0),14);   % axial  velocity of the ray


l_lb   = max(b.Z.T(:))-min(b.Z.T(:));
Z      = (b.Z.T-b.geom.z(1))/l_lb;
R      = b.R.T/l_lb;
T      = b.T.T;
dTdXi  = (T(2:end,:)-T(1:end-1,:))./[b.DXI.rm(:,1) b.DXI.rm b.DXI.rm(:,end)];
dTdEta = (T(:,2:end)-T(:,1:end-1))./[b.DETA.zm(1,:); b.DETA.zm; b.DETA.zm(end,:)];
dTdXi  = interp2(b.ETA.u,b.XI.u,dTdXi,b.ETA.T,b.XI.T);
dTdEta = interp2(b.ETA.w,b.XI.w,dTdEta,b.ETA.T,b.XI.T);
dTdR = b.XI_r.T.*dTdXi*l_lb;
dTdZ = (b.XI_z.T.*dTdXi + b.ETA_z.T.*dTdEta)*l_lb;

N    = N_coeff(1) - N_coeff(2)*(T-N_coeff(3));
dNdR = -N_coeff(2)*dTdR;
dNdZ = -N_coeff(2)*dTdZ;

% Interpolate the location of the interface to define the domain of validity of n
r_surf  = @(zEval) interp1(Z(end,:)',R(end,:)',zEval,'spline',b.R.T(end,1)/l_lb);
gradN2r = scatteredInterpolant(R(:)./r_surf(Z(:)),Z(:),N(:).*dNdR(:),'linear','linear');
gradN2z = scatteredInterpolant(R(:)./r_surf(Z(:)),Z(:),N(:).*dNdZ(:),'linear','linear');

time    = [0 t_max];
IC      = [u_in w_in r_0 z_in];
if z_p == 1
    rzuw    = [r_0 z_in u_in w_in; r_0 z_in u_in w_in];
else
    opts    = odeset('RelTol',1e-8,'AbsTol',1e-8,'Events',@(t,y)Wall(t,y,z_p,r_surf,'b'));
    [~,y]   = ode45(@(t,y) OpticalRay2D(t,y,gradN2r,gradN2z,r_surf), time, IC, opts);
    rzuw    = [y(:,3),y(:,4),y(:,1),y(:,2)];
end



ray.r_0       = r_0;
ray.z_p       = rzuw(end,2);
ray.r_p       = rzuw(end,1);
ray.alpha_0   = alpha_0;
ray.N         = N;
ray.path_real = rzuw(:,1:2);
ray.path_fict = [r_0 1; r_0+u_in/w_in*(ray.z_p-1) ray.z_p];

if ray.r_p > r_surf(ray.z_p)
    disp('The ray reached the free surface')
    ray.r_p   = NaN;
    return
end

% Compute the reflected ray
if reflection == 1
    IC    = [rzuw(end,3) -rzuw(end,4) rzuw(end,1) z_p];
    opts  = odeset('RelTol',1e-8,'AbsTol',1e-8,'Events',@(t,y)Wall(t,y,z_p,r_surf,'t'));
    [~,y] = ode45(@(t,y) OpticalRay2D(t,y,gradN2r,gradN2z,r_surf), time, IC, opts);
    rzuw  = [y(:,3),y(:,4),y(:,1),y(:,2)];

    ray.alpha_out = atan(rzuw(end,4)/rzuw(end,3));
    ray.path_real = [ray.path_real; rzuw(:,1:2)];
    ray.path_fict = [ray.path_fict; r_0-2*u_in/w_in*(1-z_p) 1];
    
end

% nested functions
    function dydt = OpticalRay2D(~,y,gradN2r,gradN2z,r_surf)
        dydt = [gradN2r(y(3)/r_surf(y(4)),y(4)); gradN2z(y(3)/r_surf(y(4)),y(4)); y(1); y(2)];

    end

    function [value,isterminal,direction] = Wall(~,y,z_p,r_surf,wall)
        % Locate the time when the value is zero
        if strcmp(wall,'b') % bottom wall
            value = y(4)-z_p;
        elseif strcmp(wall,'t') % top wall
            value = y(4)-1;
        else
            error('The parameter ''wall'' can either be ''b''(ottom) or ''t''(op).')
        end
        value = [value; 1.001*r_surf(y(4))-y(3)];
        isterminal = [1; 1];
        direction  = [0; 0];
    end

end