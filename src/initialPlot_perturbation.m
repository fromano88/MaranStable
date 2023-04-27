%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% set default parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define mycolors if it does not exist
if evalin('base','~exist(''mycolor'',''var'')')
    evalin('base','mycolor.blue = [0 0.4470 0.7410];')
    evalin('base','mycolor.red = [0.8500 0.3250 0.0980];')
    evalin('base','mycolor.ocher = [0.9290 0.6940 0.1250];')
    evalin('base','mycolor.violet = [0.4940 0.1840 0.5560];')
    evalin('base','mycolor.green = [0.4660 0.6740 0.1880];')
end

h.l_lb = evalin('base','l_lb'); % needed for non dimensionalization
try
    h.l_d1 = evalin('base','l_d1');
catch
    h.l_d1 = 0;
end

h.dim = 0;
if h.dim == 0
    h.f = 1/h.l_lb;
else
    h.f = 1000; % mm
end

% default values for contour plot
    h.blocks = evalin('base','blocks');
    if length(h.blocks) == 1
        h.ct.domain = [1 0];
    else
        h.ct.domain = [1 1];
    end
    h.ct.draw = 1; h.ct.equal = 1; h.ct.colormap = evalin('base','myColormap(''BlueWhiteRed'');');
    h.ct.r = 1; h.ct.phi = 0; h.z_T_max = evalin('base','l_lb/2+b1.geom.z(1)-z_T_max'); h.ct.z = h.z_T_max*h.f;
    if evalin('base','flowopt.energy') == 1
        h.ct.quantity = 'T';
    else
        h.ct.quantity = 'u';
    end

% default values for vector plot
    h.vp.draw = 0; h.vp.color = '[0 0 0]';
    h.vp.numPhi = 20;
    if evalin('base','length(blocks)') == 2
        h.vp.numR = 30; h.vp.numZ = 40;
    else
        h.vp.numR = 12; h.vp.numZ = 15;
    end
    
% show streamlines or temperature isolines
    h.sf_bS = 0; h.temp_bS = 0; h.prod_pF = 0; h.temp_pF = 0;
    h.iso_sf_bS   = 8; h.iso_temp_bS = 12;
    h.iso_prod_pF = 8; h.iso_temp_pF = 12;
    h.ct.color.sf_bS   = '[0 0 0]'; h.ct.color.temp_bS = '[0 0 0]';
    h.ct.color.prod_pF = '[0 0 0]'; h.ct.color.temp_pF = '[0 0 0]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% update default values %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    assignin('base','tmp',h.dim);
    evalin('base','dim=tmp;');
    
% contour plot
    assignin('base','tmp',h.ct.draw);
    evalin('base','ct.draw=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.ct.equal);
    evalin('base','ct.equal=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.ct.domain);
    evalin('base','ct.domain=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.ct.quantity);
    evalin('base','ct.quantity=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.ct.colormap);
    evalin('base','ct.colormap=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.ct.phi);
    evalin('base','ct.phi=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.ct.r);
    evalin('base','ct.r=tmp;');
    % ------------------------------------------------------------------ %
    h.ct.z_shift = h.l_lb/2+h.l_d1-h.ct.z/h.f;
    assignin('base','tmp',h.ct.z_shift);
    evalin('base','ct.z=tmp;');

% vector plot
    assignin('base','tmp',h.vp.draw);
    evalin('base','vp.draw=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.vp.numR);
    evalin('base','vp.numR=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.vp.numZ);
    evalin('base','vp.numZ=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.vp.numPhi);
    evalin('base','vp.numPhi=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.vp.color);
    evalin('base','vp.color=eval(tmp);');
    
% streamlines and isolines
    assignin('base','tmp',h.sf_bS);
    evalin('base','ct.sf_bS=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.temp_bS);
    evalin('base','ct.temp_bS=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.prod_pF);
    evalin('base','ct.prod_pF=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.temp_pF);
    evalin('base','ct.temp_pF=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.iso_sf_bS);
    evalin('base','ct.iso_sf_bS=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.iso_temp_bS);
    evalin('base','ct.iso_temp_bS=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.iso_prod_pF);
    evalin('base','ct.iso_prod_pF=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.iso_temp_pF);
    evalin('base','ct.iso_temp_pF=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.ct.color.sf_bS);
    evalin('base','ct.color.sf_bS=eval(tmp);');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.ct.color.temp_bS);
    evalin('base','ct.color.T_bS=eval(tmp);');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.ct.color.prod_pF);
    evalin('base','ct.color.prod_pF=eval(tmp);');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.ct.color.temp_pF);
    evalin('base','ct.color.T_pF=eval(tmp);');

evalin('base','clear tmp');