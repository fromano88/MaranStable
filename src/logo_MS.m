function logo_MS(flag)

% flag = 1 --> with border lines
% flag = 0 --> without border lines

if nargin == 0
    flag = 1;
end


color.M = [0.8500 0.3250 0.0980];
color.S = [0 0.4470 0.7410];

C = [0 0];

% define line via normal distance from c and steepness
    nd = 0.6;
    k  = -0.4;
    
    g  = @(x) k*x+nd*sqrt(1+k^2);

% define upper circle
    r = 0.8;                 % radius
    xs = -0.9; ys = g(xs); % starting point
    alpha = 0.7*pi;        % apex angle

    B = [xs ys];
    Cs = B+r/sqrt(1+k^2)*[-k 1];
    beta = 3*pi/2-abs(atan(k));
    phi = beta-linspace(0,alpha);
    circ_s = Cs+r*[cos(phi)' sin(phi)'];

% mirror the end point of upper circle
    E = 2*C-circ_s(end,:);

% define lower circle
    R = 1.5;

    % center M of big lower circle is obtained from the quadratic equation
    % a*M(1)^2 + b*M(1) + c = 0;
        a = 1+k^2;
        b = 2*(k*((nd-R)*sqrt(1+k^2)-E(2))-E(1));
        c = E(1)^2+((nd-R)*sqrt(1+k^2)-E(2))^2-R^2;
        M(1) = (-b-sqrt(b^2-4*a*c))/(2*a);

    M(2) = M(1)*k+nd*sqrt(1+k^2)-R;
    B(1) = M(1)-R*k/sqrt(1+k^2); B(2) = g(B(1));

    beta = pi/2-abs(atan(k));
    M = B+R/sqrt(1+k^2)*[k -1];
    alpha = acos(dot(B-M,E-M)/R^2);
    phi = beta-linspace(0,alpha);
    circ_b = M+R*[cos(phi)' sin(phi)'];

% compute middle part
    x_lin = linspace(xs,B(1));
    lin = [x_lin(2:end-1)' g(x_lin(2:end-1)')];


% merge letter S
    s = [flipud(circ_s); lin; circ_b]; % merge
    s_m = 2*C-s;                       % mirror at center
    s = [s; s_m(2:end-1,:)];           % merge again

% shift letter S such that its minimum x and y value is at 0
    shift = [-min(s(:,1)) -min(s(:,2))];
    s = s+shift;


% define circles
    rc = 0.55;
    D = Cs+[1.4 0];

    phi = linspace(0,2*pi);
    c1 =  D+rc*[cos(phi)' sin(phi)']+shift;
    c2 = -D+rc*[cos(phi)' sin(phi)']+shift;


% define letter M
    H = max(s(:,2))-min(s(:,2)); % total height of letter M
    b = 1.5*rc;                  % width of rectangle
    LL(1) = -b-0.5;              % left lower point of vertical rectangle
    w = 3+b;                     % total width of letter M
    h = 0.2*(w-2*b);             % horizontal distance from LL to middle point
    x = LL(1)-h/2;               % horizontal distance from LL to right most point of inclined quadrangle
    v = 0.3*b;                   % small vertical length of inclined rectangle in left upper corner

    LL(2) = 0;
    RU    = [LL(1)-w+2*b H];     % right upper point of inclined quadrangle
    f     = @(x) (H-v)/(w-b-h)*(RU(1)-x)+RU(2); % linear function parallel to left inlined line which goes through RU

    c3 = [RU(1)-b+rc LL(2)+rc] + rc*[cos(phi)' sin(phi)']; % circle


hold on

lw = 1.2; % linewidth
if flag == 1
    fill(s(:,1),s(:,2),color.S,'LineWidth',lw)
    fill(c1(:,1),c1(:,2),color.S,'LineWidth',lw)
    fill(c2(:,1),c2(:,2),color.S,'LineWidth',lw)
    fill([LL(1)-h x RU(1) RU(1)-b RU(1)-b],[LL(2) f(x) RU(2) RU(2) RU(2)-v],color.M,'LineWidth',lw)
    fill([LL(1) LL(1)+b LL(1)+b LL(1)],[LL(2) LL(2) LL(2)+H LL(2)+H],color.M,'LineWidth',lw)
    fill(c3(:,1),c3(:,2),color.M,'LineWidth',lw)
else
    fill(s(:,1),s(:,2),color.S,'LineStyle','none')
    fill(c1(:,1),c1(:,2),color.S,'LineStyle','none')
    fill(c2(:,1),c2(:,2),color.S,'LineStyle','none')
    fill([LL(1)-h x RU(1) RU(1)-b RU(1)-b],[LL(2) f(x) RU(2) RU(2) RU(2)-v],color.M,'LineStyle','none')
    fill([LL(1) LL(1)+b LL(1)+b LL(1)],[LL(2) LL(2) LL(2)+H LL(2)+H],color.M,'LineStyle','none')
    fill(c3(:,1),c3(:,2),color.M,'LineStyle','none')
end

axis equal
set(gca,'visible','off')
set(gcf,'color','w')