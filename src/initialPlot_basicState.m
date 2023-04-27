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

h.dim = 1; h.unity = 'mm';

% default values for contour plot
    h.ct.draw = 1; h.ct.domain = [1 1]; h.ct.equal = 1; h.ct.colormap = evalin('base','myColormap(''Cool2Warm'');');
    if evalin('base','flowopt.energy') == 1
        h.ct.quantity = 'T';
    else
        h.ct.quantity = 'u';
    end

% default values for vector plot
    h.vp.draw = 0; h.vp.color = '[0 0 0]';
    if evalin('base','length(blocks)') == 2
        h.vp.numR = 30; h.vp.numZ = 40;
    else
        h.vp.numR = 12; h.vp.numZ = 15;
    end
    
% show streamlines or temperature isolines
    h.streamlines = 0; h.temperatureIsolines = 0; h.ct.color.sf = '[0 0 0]';
    h.numStreamlines = 8; h.numTemperatureIsolines = 10; h.ct.color.T = '[0 0 0]';

% default values for line plot
    h.lp.n = 2; h.lp.axis = 'z'; h.lp.domain = [1 0]; h.lp.coord = 'curved';
    if evalin('base','flowopt.energy') == 1
        h.lp.quantity = ["T", "u"];
    else
        h.lp.quantity = ["u", "w"];
    end
    h.lp.e_r = [1 0];
    if evalin('base','length(blocks)') == 2
        h.lp.e_z = [1 0.5 0];
    else
        h.lp.e_z = [0.5 0 0];
    end
    h.lp.sum = sum(h.lp.e_r);
    h.lp.zoom = [0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% update default values %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    assignin('base','tmp',h.unity);
    evalin('base','unity=tmp;');
    % ------------------------------------------------------------------ %
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
    assignin('base','tmp',h.ct.color.sf);
    evalin('base','ct.color.sf=eval(tmp);');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.ct.color.T);
    evalin('base','ct.color.T=eval(tmp);');

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
    assignin('base','tmp',h.vp.color);
    evalin('base','vp.color=eval(tmp);');
    
% streamlines and isolines
    assignin('base','tmp',h.streamlines);
    evalin('base','ct.streamlines=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.temperatureIsolines);
    evalin('base','ct.isolines=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.numStreamlines);
    evalin('base','ct.numStreamlines=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.numTemperatureIsolines);
    evalin('base','ct.numIsolines=tmp;');

% line plot
    assignin('base','tmp',h.lp.n);
    evalin('base','lp.n=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.lp.quantity);
    evalin('base','lp.quantity=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.lp.axis);
    evalin('base','lp.axis=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.lp.domain);
    evalin('base','lp.domain=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.lp.coord);
    evalin('base','lp.coord=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.lp.e_r);
    evalin('base','lp.e_r=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.lp.e_z);
    evalin('base','lp.e_z=tmp;');
    % ------------------------------------------------------------------ %
    assignin('base','tmp',h.lp.zoom);
    evalin('base','lp.zoom=tmp;');

evalin('base','clear tmp');