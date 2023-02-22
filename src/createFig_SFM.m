%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%                  ____    ____     ____                                  %
%                 |       |    |   |    |    ___                          %
%                 |       |    |   |  __|   |___                          %
%                 |____   |____|   |  \     |___                          %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% windows or linux
    % w=1 --> windows
    % w=0 --> linux
    w = ispc;
    
% save writings
    % tabs
    TB = {'style','push','units','pix','ForegroundColor',...
        h.unselectedTabColor,'position'};
    % core panels
    CP = {'unit','pix','pos',[4 4 720 727-80*w],'visible','off'};
    % edit texts
    ET = {'style','edit','units','pix','parent'};
    % static texts
    ST = {'units','pix','style','text','HorizontalAlignment','left','pos'};
    % axes
    AS = {'unit','pix','visible','off','parent'};
    
% create tabs
    tab_width = [72 83 57 85 156 95 98 75 62];
    h.pb.general = uicontrol(TB{:},   [3                     733 tab_width(1) 27],'string','General','ForegroundColor',h.selectedTabColor);
    h.pb.geometry = uicontrol(TB{:},  [3+sum(tab_width(1:1)) 733 tab_width(2) 27],'string','Geometry');
    h.pb.mesh = uicontrol(TB{:},      [3+sum(tab_width(1:2)) 733 tab_width(3) 27],'string','Mesh');
    h.pb.equations = uicontrol(TB{:}, [3+sum(tab_width(1:3)) 733 tab_width(4) 27],'string','Equations');
    h.pb.bcs = uicontrol(TB{:},       [3+sum(tab_width(1:4)) 733 tab_width(5) 27],'string','Boundary Conditions');
    h.pb.simulation = uicontrol(TB{:},[3+sum(tab_width(1:5)) 733 tab_width(6) 27],'string','Basic State');
    h.pb.ray = uicontrol(TB{:},       [3+sum(tab_width(1:6)) 733 tab_width(7) 27],'string','Optical Ray');
    h.pb.lsa = uicontrol(TB{:},       [3+sum(tab_width(1:7)) 733 tab_width(8) 27],'string','Stability');
    h.pb.about = uicontrol(TB{:},     [3+sum(tab_width(1:8)) 733 tab_width(9) 27],'string','About');
    if ispc %windows
        tab_width = [60 68 50 70 120 75 76 60 52];
        set(h.pb.general,'position',   [3                     653 tab_width(1) 25])
        set(h.pb.geometry,'position',  [3+sum(tab_width(1:1)) 653 tab_width(2) 25])
        set(h.pb.mesh,'position',      [3+sum(tab_width(1:2)) 653 tab_width(3) 25])
        set(h.pb.equations,'position', [3+sum(tab_width(1:3)) 653 tab_width(4) 25])
        set(h.pb.bcs,'position',       [3+sum(tab_width(1:4)) 653 tab_width(5) 25])
        set(h.pb.simulation,'position',[3+sum(tab_width(1:5)) 653 tab_width(6) 25])
        set(h.pb.ray,'position',       [3+sum(tab_width(1:6)) 653 tab_width(7) 25])
        set(h.pb.lsa,'position',       [3+sum(tab_width(1:7)) 653 tab_width(8) 25])
        set(h.pb.about,'position',     [3+sum(tab_width(1:8)) 653 tab_width(9) 25])
    end

% summarize all to save writings (later)
    h.tabs = {h.pb.general h.pb.geometry h.pb.mesh h.pb.equations ...
        h.pb.bcs h.pb.simulation h.pb.ray h.pb.lsa h.pb.about};

% sketch of the liquid bridge & logo
    h.pl.sketch = uipanel('unit','pix','pos',[726 4 403 727],...
        'BackgroundColor',[1 1 1]);
    h.as.sketch = axes(AS{:},h.pl.sketch,'pos',[19 8 385 709]);
    h.as.logo   = axes('unit','pix','visible','off','pos',[1087 733 40 25]);
    if ispc
        set(h.pl.sketch,'pos',[726 4 373 647])
        set(h.as.sketch,'pos',[19 8 355 629])
        set(h.as.logo,'pos',[1057 653 40 25])
    end
    
% create main panels
    h.pl.general = uipanel(CP{:},'visible','on');
    h.pl.geometry = uipanel(CP{:});
    h.pl.mesh = uipanel(CP{:});
    h.pl.equations = uipanel(CP{:});
    h.pl.bcs = uipanel(CP{:});
    h.pl.simulation = uipanel(CP{:});
    h.pl.ray = uipanel(CP{:});
    h.pl.lsa = uipanel(CP{:});
    h.pl.about = uipanel(CP{:});
    
% summarize all to save writings (later)
    h.mainPanels = {h.pl.general h.pl.geometry h.pl.mesh h.pl.equations ...
        h.pl.bcs h.pl.simulation h.pl.ray h.pl.lsa h.pl.about};
    
% creating info panel, when user comes across the logo with the mouse
    if ispc
        font_name = 'Times';
    else
        font_name = 'Likhan';
    end 
    SA = {'units','pix','style','text','HorizontalAlignment','left','BackgroundColor',[1 1 1],'FontName',font_name,'FontSize',10,'pos'};
    h.pl.info(1) = uipanel(CP{:},'pos',[727 350 401 410],...
        'BackgroundColor',[1 1 1]);
    h.pl.info(2) = uipanel('parent',h.pl.info(1),'units','pix',...
        'BackgroundColor',h.unselectedTabColor,'pos',[0 379 401 48]);
    if ispc
        set(h.pl.info(1),'pos',[727 230 370 420])
        set(h.pl.info(2),'pos',[0 385 368 34])
    end
    uicontrol('parent',h.pl.info(2),'unit','pix','style','text','pos',...
        [1 6+2*w 398-15*w 22],'BackgroundColor',h.unselectedTabColor,'ForegroundColor',[1 1 1],...
        'string','MaranStable','FontSize',16,'FontWeight','bold','FontName','Likhan');
    uicontrol('parent',h.pl.info(1),SA{:},[5 355 70 17],'string','University:');
    uicontrol('parent',h.pl.info(1),SA{:},[15 335 60 17],'string','TU Wien');
    uicontrol('parent',h.pl.info(1),SA{:},[5 295 70 17],'string','Department:');
    uicontrol('parent',h.pl.info(1),SA{:},[15 275 270-10*w 17],'string','Institute of Fluid Mechanics and Heat Transfer');
    uicontrol('parent',h.pl.info(1),SA{:},[5 235 140 17],...
        'BackgroundColor',[1 1 1],'string','Developed by:');
    uicontrol('parent',h.pl.info(1),SA{:},[15 215 140 17],...
        'BackgroundColor',[1 1 1],'string','Mario Stojanovic');
    uicontrol('parent',h.pl.info(1),SA{:},[15 195 140 17],...
        'BackgroundColor',[1 1 1],'string','Michael Lukasser');
    uicontrol('parent',h.pl.info(1),SA{:},[15 175 140 17],...
        'BackgroundColor',[1 1 1],'string','Francesco Romano');
    uicontrol('parent',h.pl.info(1),SA{:},[15 155 140 17],...
        'BackgroundColor',[1 1 1],'string','Hendrik C. Kuhlmann');
    uicontrol('parent',h.pl.info(1),SA{:},[5 115 140 17],'string','Sponsored by:');
    uicontrol('parent',h.pl.info(1),SA{:},[15 95 220-20*w 17],'string','Austrian Research Promotion Agency');
    uicontrol('parent',h.pl.info(1),SA{:},[15 75 140-10*w 17],'string','European Space Agency');
    uicontrol('parent',h.pl.info(1),SA{:},[5 35 80 17],'string','License:');
    uicontrol('parent',h.pl.info(1),SA{:},[15 15 100-15*w 17],'string','CC BY-NC-SA');
    axes(AS{:},h.pl.info(1),'pos',[80 328+1*w 30 30]);
        imshow(imread('logo_TU.png'));
    axes(AS{:},h.pl.info(1),'pos',[292-20*w 262+2*w 45 45]);
        imshow(imread('logo_isw.png'));
    axes(AS{:},h.pl.info(1),'pos',[238-12*w 77+1*w 50 50]);
        imshow(imread('logo_ffg.png'));
    axes(AS{:},h.pl.info(1),'pos',[162-5*w 58 50 50]);
        imshow(imread('logo_esa.png'));
    axes(AS{:},h.pl.info(1),'pos',[115-2*w -7+2*w 100 60]);
        imshow(imread('license_white.png'));
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%           ____                           ____                           %
%          |         ___    |\  |   ___   |    |    /\    |               %
%          |  __    |___    | \ |  |___   |  __|   /__\   |               %
%          |____|   |___    |  \|  |___   |  \    /    \  |___            %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subpanels
    h.gl.pl.g = uipanel('parent',h.pl.general,'units','pix','pos',...
        [40 470-80*w 300 90],'title','Gravitational Acceleration');
    h.gl.pl.Vr = uipanel('parent',h.pl.general,'units','pix','pos',...
        [40 360-80*w 300 90],'title','Volume Ratio');
    h.gl.pl.fluidProperties = uipanel('parent',h.pl.general,'units',...
        'pix','pos',[380 360-80*w 300 330],'title','Fluid Properties');

% button groups
    h.gl.bg.formulation = uibuttongroup('parent',h.pl.general,'units',...
        'pix','pos',[40 580-80*w 300 110],'title','2D Space');
    h.gl.bg.model = uibuttongroup('parent',...
        h.gl.pl.fluidProperties,'units','pix','pos',[-5 189 310 155]);
    
% radio buttons
    h.gl.rb.planar = uicontrol(h.gl.bg.formulation,'style','rad','unit',...
        'pix','pos',[26 60 150 20],'string','Planar');
    h.gl.rb.axisymmetric = uicontrol(h.gl.bg.formulation,'style','rad',...
        'unit','pix','pos',[26 20 150 20],'value',1,'string','Axisymmetric');
    h.gl.rb.NS1 = uicontrol(h.gl.bg.model,'style','rad',...
        'unit','pix','pos',[16+15*w 90 260 20],'string','Oberbeck-Boussinesq Approximation');
    h.gl.rb.NS2 = uicontrol(h.gl.bg.model,'style','rad','unit','pix',...
        'pos',[16+15*w 55 260 20],'string','Linearly Temperature Dependent');
    h.gl.rb.NS3 = uicontrol(h.gl.bg.model,'style','rad','unit','pix',...
        'pos',[16+15*w 20 260 20],'string','Fully Temperature Dependent');
    switch h.NS
        case 1
            set(h.gl.rb.NS1,'value',1)
        case 2
            set(h.gl.rb.NS2,'value',1)
        case 3
            set(h.gl.rb.NS3,'value',1)
    end
    
% edit texts
    h.gl.et.g = uicontrol('parent',h.gl.pl.g,'style','edit','units',...
        'pix','pos',[62 27 60 23],'string',num2str(h.flowopt.g));
    h.gl.et.Vr = uicontrol('parent',h.gl.pl.Vr,'style','edit','units',...
        'pix','pos',[120 27 60 23],'string',num2str(h.V_r));

% pushbuttons
    h.gl.pb.create = uicontrol('parent',h.gl.pl.fluidProperties,'units','pix',...
        'pos',[198+12*w 6 93-11*w 25],'string','Create Fluid');
    
% popup menus
    str = {'KF-96L-5cs','KF-96L-2cs','KF-96L-20cs','H2O','Br2','HF',...
        'HCN','Others ...'};
    h.gl.pm.liquid = uicontrol('parent',h.gl.pl.fluidProperties,'style',...
        'pop','units','pix','pos',[35 126 130 24],'string',str);
    if max(contains(str,h.selectedLiquid))
        idx = find(ismember(str,h.selectedLiquid));
        set(h.gl.pm.liquid,'Value',idx)
    else
        str=[str(1:end-1), h.selectedLiquid, str(end)];
        set(h.gl.pm.liquid,'String',str,'value',length(str)-1);
    end
    
% axes (for latex texts)
    h.gl.as.g = axes(AS{:},h.gl.pl.g,'pos',h.gl.et.g.Position);
    h.gl.as.Vr = axes(AS{:},h.gl.pl.Vr,'pos',h.gl.et.Vr.Position);
    
% static texts
    h.gl.st.liquid = uicontrol(ST{:},[26 155 100 17],'parent',...
        h.gl.pl.fluidProperties,'string','Select Liquid');
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%        ____             ___                 _____    ____               %
%       |         ___    |   |   |\  /|   ___   |     |    |   \  /       %
%       |  __    |___    |   |   | \/ |  |___   |     |  __|    \/        %
%       |____|   |___    |___|   |    |  |___   |     |  \       \        %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% subpanels
    h.gy.pl.lengths = uipanel('parent',h.pl.geometry,'units','pix',...
        'pos',[20 600-80*w 680 100],'title','Length');
    h.gy.pl.radii = uipanel('parent',h.pl.geometry,'units','pix','pos',...
        [20 434-80*w 680 136],'title','Radii');
    
% edit texts
    h.gy.et.llb = uicontrol('parent',h.gy.pl.lengths,'style','edit',...
        'units','pix','pos',[90 33 60 23],'string',num2str(h.l_lb*1000));
    h.gy.et.rc = uicontrol('parent',h.gy.pl.radii,'style','edit',...
        'units','pix','pos',[90 76 60 23],'string',num2str(h.r_c*1000));
    h.gy.et.ri = uicontrol('parent',h.gy.pl.radii,'style','edit',...
        'units','pix','pos',[90 30 60 23],'string',num2str(h.r_i*1000));
    
% axes (for latex texts)
    h.gy.as.llb = axes(AS{:},h.gy.pl.lengths,'pos',h.gy.et.llb.Position);
    h.gy.as.rc = axes(AS{:},h.gy.pl.radii,'pos',h.gy.et.rc.Position);
    h.gy.as.ri = axes(AS{:},h.gy.pl.radii,'pos',h.gy.et.ri.Position);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%                                                                         %
%                      |\  /|   ___   ____                                %
%                      | \/ |  |___  |____   |__|                         %
%                      |    |  |___   ____|  |  |                         %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save writings
    % panels
    P = {'units','pix','pos',[0 260-65*w 720 415],'BorderType','none',...
        'visible','off','parent'};
    % panels: spacing, coeffiecient, function
    PS = {'units','pix','pos',[10 280 700 132],'title','Spacing','parent'};
    PC = {'units','pix','pos',[10 143 700 132],'title','Stretching Factor','parent'};
    PF = {'units','pix','pos',[10 6 700 132],'title','Stretching Function','parent'};
    % subpanels: coefficient
    SC = {'units','pix','pos',[234 -5 233 150],'parent'};
    % button (groups): start, fit, end
    BS = {'units','pix','BorderType','none','pos',[20 2 190 90],'parent'};
    BF = {'units','pix','BorderType','none','pos',[20 6 190 90],'parent'};
    BE = {'units','pix','BorderType','none','pos',[487 2 190 90],'parent'};
    BG = {BS{:}; BF{:}; BE{:}};
    % radio (button): tangent hyp., geometric progression, linear
    RT = {'style','rad','unit','pix','pos',[5 65 180 20],'string',...
        'Hyperbolic Tangent'};
    RG = {'style','rad','unit','pix','pos',[5 37 180 20],'string',...
        'Geometric Progression'};
    RL = {'style','rad','unit','pix','pos',[5 9 180 20],'string','Linear'};
    % start - fit - end
    SFT = {'start', 'fit', 'end'};
    
% subpanels
    h.mh.pl.radialLB.main = uipanel(P{:},h.pl.mesh,'visible','on');
    h.mh.pl.axialLB.main  = uipanel(P{:},h.pl.mesh);
    
% subsubpanels
    h.mh.pl.radialLB.spacing   = uipanel(PS{:},h.mh.pl.radialLB.main);
    h.mh.pl.axialLB.spacing    = uipanel(PS{:},h.mh.pl.axialLB.main);
    % --------------------------------------------------------------------
    h.mh.pl.radialLB.stretchingC    = uipanel(PC{:},h.mh.pl.radialLB.main);
    h.mh.pl.axialLB.stretchingC     = uipanel(PC{:},h.mh.pl.axialLB.main);
    % --------------------------------------------------------------------
    h.mh.pl.radialLB.stretchingF(1) = uipanel(PF{:},h.mh.pl.radialLB.main);
    h.mh.pl.axialLB.stretchingF(1)  = uipanel(PF{:},h.mh.pl.axialLB.main);
    
% subsubsubpanels
    h.mh.pl.radialLB.stretchingF(2) = uipanel(SC{:},h.mh.pl.radialLB.stretchingF(1));
    h.mh.pl.axialLB.stretchingF(2)  = uipanel(SC{:},h.mh.pl.axialLB.stretchingF(1));
    
% button groups
    for i = 1:3
        try
            h.mh.bg.radialLB(i) = uibuttongroup(BG{i,:},h.mh.pl.radialLB.stretchingF(i));
            h.mh.bg.axialLB(i) = uibuttongroup(BG{i,:},h.mh.pl.axialLB.stretchingF(i));
        catch
            h.mh.bg.radialLB(i) = uibuttongroup(BG{i,:},h.mh.pl.radialLB.stretchingF(i-2));
            h.mh.bg.axialLB(i) = uibuttongroup(BG{i,:},h.mh.pl.axialLB.stretchingF(i-2));
        end
    end
    % --------------------------------------------------------------------
    h.mh.buttonGroups.radial = {h.mh.bg.radialLB(1) h.mh.bg.radialLB(3)};
    
    h.mh.buttonGroups.axial = {h.mh.bg.axialLB(1) h.mh.bg.axialLB(3)};
    
% radio buttons
    for i = 1:3
        h.mh.rb.radialLB(i).tanh = uicontrol(h.mh.bg.radialLB(i),RT{:});
        h.mh.rb.radialLB(i).gp = uicontrol(h.mh.bg.radialLB(i),RG{:});
        h.mh.rb.radialLB(i).lin = uicontrol(h.mh.bg.radialLB(i),RL{:});
        % ----------------------------------------------------------------
        h.mh.rb.axialLB(i).tanh = uicontrol(h.mh.bg.axialLB(i),RT{:});
        h.mh.rb.axialLB(i).gp = uicontrol(h.mh.bg.axialLB(i),RG{:});
        h.mh.rb.axialLB(i).lin = uicontrol(h.mh.bg.axialLB(i),RL{:});
        % ----------------------------------------------------------------
        if strcmp(h.mesh.r.sf,'tanh')
            set(h.mh.rb.radialLB(i).tanh,'value',1)
        else
            set(h.mh.rb.radialLB(i).gp,'value',1)
        end
        if strcmp(h.mesh.z.sf,'tanh')
            set(h.mh.rb.axialLB(i).tanh,'value',1)
        else
            set(h.mh.rb.axialLB(i).gp,'value',1)
        end
    end
    % --------------------------------------------------------------------
    set([h.mh.rb.radialLB(2).tanh h.mh.rb.axialLB(2).tanh],'enable','off','value',0)
    set([h.mh.rb.radialLB(2).gp h.mh.rb.axialLB(2).gp],'enable','off')
    set([h.mh.rb.radialLB(2).lin h.mh.rb.axialLB(2).lin],'value',1)
    
% edit texts
    for i = 1:3
        h.mh.et.radialLB.spacing(i) = uicontrol('pos',...
            [325 (120-35*i) 60 23],ET{:},h.mh.pl.radialLB.spacing,...
            'string',eval(['num2str(h.mesh.r.delta.' SFT{i} ')']));
        h.mh.et.axialLB.spacing(i) = uicontrol('pos',...
            [325 (120-35*i) 60 23],ET{:},h.mh.pl.axialLB.spacing,...
            'string',eval(['num2str(h.mesh.z.delta.' SFT{i} ')']));
    end
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    h.mh.et.radialLB.stretchingC(1) = uicontrol('pos',[325 85 60 23],ET{:},...
        h.mh.pl.radialLB.stretchingC,'string',num2str(h.mesh.r.f(1)));
    h.mh.et.radialLB.stretchingC(2) = uicontrol('pos',[325 50 60 23],ET{:},...
        h.mh.pl.radialLB.stretchingC,'string',num2str(1),'enable','off');
    h.mh.et.radialLB.stretchingC(3) = uicontrol('pos',[325 15 60 23],ET{:},...
        h.mh.pl.radialLB.stretchingC,'string',num2str(h.mesh.r.f(2)));
    % --------------------------------------------------------------------
    h.mh.et.axialLB.stretchingC(1) = uicontrol('pos',[325 85 60 23],ET{:},...
        h.mh.pl.axialLB.stretchingC,'string',num2str(h.mesh.z.f(1)));
    h.mh.et.axialLB.stretchingC(2) = uicontrol('pos',[325 50 60 23],ET{:},...
        h.mh.pl.axialLB.stretchingC,'string',num2str(1),'enable','off');
    h.mh.et.axialLB.stretchingC(3) = uicontrol('pos',[325 15 60 23],ET{:},...
        h.mh.pl.axialLB.stretchingC,'string',num2str(h.mesh.z.f(2)));    
    
% pushbuttons
    h.mh.pb.checkMesh = uicontrol('parent',h.pl.mesh,'units','pix',...
        'pos',[515 227-63*w 93 25],'string','Check Mesh');
    h.mh.pb.drawMesh = uicontrol('parent',h.pl.mesh,'units','pix',...
        'pos',[618 227-63*w 90 25],'string','Draw Mesh');
    
% popup menus
    h.mh.pm.headline = uicontrol('parent',h.pl.mesh,'style','pop',...
        'pos',[373 684-75*w 88 24],'string',{'radial','axial'});
    
% axes (for latex texts)
    h.mh.as.spacing(1) = axes(AS{:},h.mh.pl.radialLB.spacing,'pos',...
        [62 27 60 23]);
    for i = 1:3
        h.mh.as.radialLB.spacing(i) = axes(AS{:},...
            h.mh.pl.radialLB.spacing,'pos',h.mh.et.radialLB.spacing(i).Position);
        h.mh.as.axialLB.spacing(i) = axes(AS{:},...
            h.mh.pl.axialLB.spacing,'pos',h.mh.et.axialLB.spacing(i).Position);
        % ----------------------------------------------------------------
        h.mh.as.radialLB.stretchingC(i) = axes(AS{:},...
            h.mh.pl.radialLB.stretchingC,'pos',h.mh.et.radialLB.stretchingC(i).Position);
        h.mh.as.axialLB.stretchingC(i) = axes(AS{:},...
            h.mh.pl.axialLB.stretchingC,'pos',h.mh.et.axialLB.stretchingC(i).Position);
    end
    h.mh.as.radialLB.stretchingF(1) = axes(AS{:},...
        h.mh.pl.radialLB.stretchingF(1),'pos',[12 92 20 17]);
    h.mh.as.radialLB.stretchingF(2) = axes(AS{:},...
        h.mh.pl.radialLB.stretchingF(2),'pos',[12 96 20 17]);
    h.mh.as.radialLB.stretchingF(3) = axes(AS{:},...
        h.mh.pl.radialLB.stretchingF(1),'pos',[477 92 20 17]);
    % --------------------------------------------------------------------
    h.mh.as.axialLB.stretchingF(1) = axes(AS{:},...
        h.mh.pl.axialLB.stretchingF(1),'pos',[12 92 20 17]);
    h.mh.as.axialLB.stretchingF(2) = axes(AS{:},...
        h.mh.pl.axialLB.stretchingF(2),'pos',[12 96 20 17]);
    h.mh.as.axialLB.stretchingF(3) = axes(AS{:},...
        h.mh.pl.axialLB.stretchingF(1),'pos',[477 92 20 17]);
    
% static texts
    h.mh.st.headline1 = uicontrol(ST{:},[60 687-75*w 310 17],'parent',...
        h.pl.mesh,'HorizontalAlignment','right','string',...
        'Liquid Phase: Mesh Parameters in');
    h.mh.st.headline2 = uicontrol(ST{:},[464 687-75*w 100 17],'parent',...
        h.pl.mesh,'string','Direction');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%            ____                    _____        ____                    %
%     ___   |    |   |   |     /\      |     |   |    |   |\  |    ___    %
%    |___   |    |   |   |    /__\     |     |   |    |   | \ |   |___    %
%    |___   |___\|   |___|   /    \    |     |   |____|   |  \|    ___|   %
%                      \                                                  %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save writings
    % panels
    P = {'units','pix','parent',h.pl.equations,'pos'};
    
% subpanels
    h.es.pl.continuity = uipanel(P{:},[20 640-80*w 680 60],'title','Continuity Equation');
    h.es.pl.momentum   = uipanel(P{:},[20 560-80*w 680 60],'title','Momentum Equation');
    h.es.pl.energy     = uipanel(P{:},[20 480-80*w 680 60],'title','Energy Equation');
    h.es.pl.additional = uipanel(P{:},[20 365-80*w 680 95],'title','Additional Options');
        
% check boxes
    h.es.cb.energy = uicontrol('parent',h.es.pl.energy,'style',...
        'check','pos',[30 14 85 20],'String','On','Value',1);
    h.es.cb.marangoni = uicontrol('parent',h.es.pl.additional,'style',...
        'check','pos',[30 45 210 20],'value',1,'String','Thermocapillary Convection');
    h.es.cb.stokes = uicontrol('parent',h.es.pl.additional,'style',...
        'check','pos',[30 15 210 20],'String','Stokes Flow (Creeping Flow)');

% axes
    h.es.as.continuity = axes(AS{:},h.es.pl.continuity,'pos',[320 23 40 25]);
    h.es.as.momentum   = axes(AS{:},h.es.pl.momentum,  'pos',[320 23 40 25]);
    h.es.as.energy     = axes(AS{:},h.es.pl.energy,    'pos',[320 23 40 25]);
    h.es.as.properties = axes(AS{:},h.pl.equations,    'pos',[22 50 40 25]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%                  ____           _____                                   %
%                 |    |         |                                        %
%                 |  __|         |               ___                      %
%                 |    |    _    |         _    |___     _                %
%                 |____|   |_|   |_____   |_|    ___|   |_|               %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save writings
    % panels, subpanels, subsubpanels
    P = {'units','pix','pos',[0 0 720 670],'BorderType','none',...
        'visible','off','parent'};
    SP1 = {'units','pix','pos',[20 480-40*w 680 160-20*w],'parent'};
    SP2 = {'units','pix','pos',[20 260-20*w 680 160-20*w],'parent'};
    SSP1 = {'units','pix','pos',[227 -5 226 170-20*w],'parent'};
    SSP2 = {'units','pix','pos',[-5 -5 690 55],'parent'};
    % button (groups):
    BV = {'units','pix','BorderType','none','pos',[15 70-10*w 195 80],'parent'};
    BH = {'units','pix','BorderType','none','pos',[467 66-10*w 195 80],'parent'};
    % radio (button): wall, noSlip, slip, adiabatic, conductive
    RW = {'style','rad','unit','pix','pos',[25 95-10*w 50 20],'string',...
        'Wall','value',1,'enable','inactive','parent'};
    RN = {'style','rad','unit','pix','pos',[10 50 110 20],'string',...
        'No Slip','value',1};
    RS = {'style','rad','unit','pix','pos',[10 10 110 20],'string',...
        'Slip'};
    RA = {'style','rad','unit','pix','pos',[10 50 110 20],'string','Adiabatic'};
    RC = {'style','rad','unit','pix','pos',[10 10 110 20],'string',...
        'Conductive','value',1};
    % edit texts
    E = {'style','edit','units','pix','pos',[90 18 142 23],'parent'};
    % popup menu
    PM = {'style','pop','units','pix','pos',[540 17 88 24],...
        'string',{'constant','variable'},'parent'};
    % (statix) texts
    TT = {'units','pix','style','text','HorizontalAlignment','left',...
        'pos',[400+25*w 20 135-25*w 17],'string','Temperature Profile:','parent'};
    
% subpanels
    h.bc.pl.radialLB.main = uipanel(P{:},h.pl.bcs,'visible','on');
    h.bc.pl.axialLB.main  = uipanel(P{:},h.pl.bcs);
    
% subsubpanels
    h.bc.pl.axialLB.d1(1) = uipanel(SP1{:},h.bc.pl.axialLB.main);
    h.bc.pl.axialLB.d2(1) = uipanel(SP2{:},h.bc.pl.axialLB.main);
    h.bc.pl.radialLB.rc(1) = uipanel(SP1{:},h.bc.pl.radialLB.main);
    h.bc.pl.radialLB.ri(1) = uipanel(SP2{:},h.bc.pl.radialLB.main);
    h.bc.pl.radialLB.ri(3) = uipanel(SP2{:},h.bc.pl.radialLB.main,'pos',[20 211-20*w 680 51]);
    
% subsubsubpanels
    h.bc.pl.axialLB.d1(2) = uipanel(SSP1{:},h.bc.pl.axialLB.d1(1));
    h.bc.pl.axialLB.d1(3) = uipanel(SSP2{:},h.bc.pl.axialLB.d1(1));
    h.bc.pl.axialLB.d2(2) = uipanel(SSP1{:},h.bc.pl.axialLB.d2(1));
    h.bc.pl.axialLB.d2(3) = uipanel(SSP2{:},h.bc.pl.axialLB.d2(1));
    % --------------------------------------------------------------------
    h.bc.pl.radialLB.rc(2) = uipanel(SSP1{:},h.bc.pl.radialLB.rc(1));
    h.bc.pl.radialLB.rc(3) = uipanel(SSP2{:},h.bc.pl.radialLB.rc(1));
    h.bc.pl.radialLB.ri(2) = uipanel(SSP1{:},h.bc.pl.radialLB.ri(1),'pos',[451 -5 232 170]);
    
% popup menus
    h.bc.pm.headline = uicontrol('parent',h.pl.bcs,'style','pop',...
        'pos',[373 683-75*w 88 24],'string',{'radial','axial'});
    
% button groups
    h.bc.bg.axialLB.d1(2) = uibuttongroup(BV{:},h.bc.pl.axialLB.d1(2));
    h.bc.bg.axialLB.d1(3) = uibuttongroup(BH{:},h.bc.pl.axialLB.d1(1));
    h.bc.bg.axialLB.d2(2) = uibuttongroup(BV{:},h.bc.pl.axialLB.d2(2));
    h.bc.bg.axialLB.d2(3) = uibuttongroup(BH{:},h.bc.pl.axialLB.d2(1));
    % --------------------------------------------------------------------
    h.bc.bg.radialLB.rc(1) = uibuttongroup(BV{:},h.bc.pl.radialLB.rc(1),'pos',[15 66-10*w 195 80]);
    h.bc.bg.radialLB.rc(2) = uibuttongroup(BV{:},h.bc.pl.radialLB.rc(2));
    h.bc.bg.radialLB.rc(3) = uibuttongroup(BH{:},h.bc.pl.radialLB.rc(1));
    h.bc.bg.radialLB.ri(1) = uibuttongroup(BV{:},h.bc.pl.radialLB.ri(1),...
        'pos',[15 16 400 125]);
    h.bc.bg.radialLB.ri(2) = uibuttongroup(BV{:},h.bc.pl.radialLB.ri(2),...
        'pos',[15 42-5*w 195 80]);
    % --------------------------------------------------------------------
    h.bc.buttonGroups.velocity = {h.bc.bg.axialLB.d1(2) ...
        h.bc.bg.axialLB.d2(2) h.bc.bg.radialLB.rc(2)};
    h.bc.buttonGroups.temperature = {h.bc.bg.axialLB.d1(3) ...
        h.bc.bg.axialLB.d2(3) h.bc.bg.radialLB.rc(3) h.bc.bg.radialLB.ri(2)};
    
% radio buttons
    uicontrol(RW{:},h.bc.pl.axialLB.d1(1));
    uicontrol(RW{:},h.bc.pl.axialLB.d2(1));
    % --------------------------------------------------------------------
    h.bc.rb.axialLB.d1.slip    = uicontrol(h.bc.bg.axialLB.d1(2),RS{:});
    h.bc.rb.axialLB.d2.slip    = uicontrol(h.bc.bg.axialLB.d2(2),RS{:});
    h.bc.rb.radialLB.rc.slip   = uicontrol(h.bc.bg.radialLB.rc(2),RS{:});
    h.bc.rb.radialLB.rc.noSlip = uicontrol(h.bc.bg.radialLB.rc(2),RN{:});
    h.bc.rb.axialLB.d1.noSlip  = uicontrol(h.bc.bg.axialLB.d1(2),RN{:});
    h.bc.rb.axialLB.d2.noSlip  = uicontrol(h.bc.bg.axialLB.d2(2),RN{:});
    % --------------------------------------------------------------------
    h.bc.rb.axialLB.d1.adiabatic  = uicontrol(h.bc.bg.axialLB.d1(3),RA{:});
    h.bc.rb.axialLB.d1.conductive = uicontrol(h.bc.bg.axialLB.d1(3),RC{:});
    h.bc.rb.axialLB.d2.adiabatic  = uicontrol(h.bc.bg.axialLB.d2(3),RA{:});
    h.bc.rb.axialLB.d2.conductive = uicontrol(h.bc.bg.axialLB.d2(3),RC{:});
    % --------------------------------------------------------------------
    h.bc.rb.radialLB.rc.wall = uicontrol(h.bc.bg.radialLB.rc(1),...
        RS{:},'string','Wall','enable','off');
    h.bc.rb.radialLB.rc.axis = uicontrol(h.bc.bg.radialLB.rc(1),...
        RN{:},'string','Symmetry','enable','inactive');
    
    h.bc.rb.radialLB.rc.adiabatic  = uicontrol(h.bc.bg.radialLB.rc(3),RA{:});
    h.bc.rb.radialLB.rc.conductive = uicontrol(h.bc.bg.radialLB.rc(3),RC{:});
    h.bc.rb.radialLB.ri.rigid  = uicontrol(h.bc.bg.radialLB.ri(1),...
        RA{:},'string','Straight indeformable surface shape','pos',[10 90-5*w 350 20]);
    h.bc.rb.radialLB.ri.static  = uicontrol(h.bc.bg.radialLB.ri(1),...
        RA{:},'string','Indeformable hydrostatic surface shape','pos',[10 50-5*w 350 20]);
    h.bc.rb.radialLB.ri.dynamic  = uicontrol(h.bc.bg.radialLB.ri(1),...
        RA{:},'string','Dynamically deformed surface shape','pos',[10 10-5*w 350 20]);
    if strcmp(h.b1.bc.z(2,2),'s')
        set(h.bc.rb.radialLB.ri.dynamic,'value',1)
    elseif strcmp(h.b1.bc.z(2,2),'i')
        set(h.bc.rb.radialLB.ri.static,'value',1)
    else
        set(h.bc.rb.radialLB.ri.rigid,'value',1)
        set(h.gl.et.Vr,'String','1','Enable','off')
        h.V_r = 1;
    end
    h.bc.rb.radialLB.ri.conductive = uicontrol(h.bc.bg.radialLB.ri(2),...
        RC{:},'value',0,'string','Newton''s law');
    h.bc.rb.radialLB.ri.adiabatic  = uicontrol(h.bc.bg.radialLB.ri(2),...
        RA{:},'value',1);
    set([h.bc.rb.radialLB.rc.noSlip h.bc.rb.radialLB.rc.slip h.bc.rb.radialLB.rc.adiabatic h.bc.rb.radialLB.rc.conductive],...
        'value',0,'enable','off')
    
% edit texts
    h.bc.et.axialLB.d1.temperature = uicontrol(E{:},...
        h.bc.pl.axialLB.d1(3),'string',num2str(h.T_d1l));
    h.bc.et.axialLB.d2.temperature = uicontrol(E{:},...
        h.bc.pl.axialLB.d2(3),'string',num2str(h.T_d2l));
    % --------------------------------------------------------------------
    h.bc.et.radialLB.rc.temperature = uicontrol(E{:},...
        h.bc.pl.radialLB.rc(3),'string',[num2str((h.T_d1l+h.T_d2l)*0.5) '+' num2str(h.T_d1l-h.T_d2l) '*z/d'],...
        'visible','off','enable','inactive');
    h.bc.et.radialLB.ri.Bi = uicontrol(E{:},...
        h.bc.pl.radialLB.ri(3),'string',num2str(h.Bi),'pos',[90 13 60 23],'visible','off');
    
% popup menus
    h.bc.pm.axialLB.d1.temperature = uicontrol(PM{:},h.bc.pl.axialLB.d1(3));
    h.bc.pm.axialLB.d2.temperature = uicontrol(PM{:},h.bc.pl.axialLB.d2(3));
    % --------------------------------------------------------------------
    h.bc.pm.radialLB.rc.temperature = uicontrol(PM{:},h.bc.pl.radialLB.rc(3),...
        'value',2,'visible','off');
    % --------------------------------------------------------------------
    h.bc.popupMenus.temperature = {h.bc.pm.axialLB.d1.temperature ...
        h.bc.pm.axialLB.d2.temperature h.bc.pm.radialLB.rc.temperature};

% axes
    h.bc.as.axialLB.d1(1)  = axes(AS{:},h.bc.pl.axialLB.main,'pos',[20 640-60*w 20 17]);
    h.bc.as.axialLB.d2(1)  = axes(AS{:},h.bc.pl.axialLB.main,'pos',[20 420-40*w 20 17]);
    h.bc.as.radialLB.rc(1) = axes(AS{:},h.bc.pl.radialLB.main,'pos',[20 640-60*w 20 17]);
    h.bc.as.radialLB.ri = axes(AS{:},h.bc.pl.radialLB.main,'pos',[20 420-40*w 20 17]);
    h.bc.as.radialLB.Bi = axes(AS{:},h.bc.pl.radialLB.ri(3),'pos',[90 13 60 23]);
    % --------------------------------------------------------------------
    h.bc.as.axialLB.d1(2)  = axes(AS{:},h.bc.pl.axialLB.d1(3),'pos',[90 18 142 23]);
    h.bc.as.axialLB.d2(2)  = axes(AS{:},h.bc.pl.axialLB.d2(3),'pos',[90 18 142 23]);
    h.bc.as.radialLB.rc(2)  = axes(AS{:},h.bc.pl.radialLB.rc(3),'pos',[90 18 142 23]);    

% static texts
    h.bc.st.headline1 = uicontrol(ST{:},[60 687-75*w 310 17],'parent',...
        h.pl.bcs,'HorizontalAlignment','right','string',...
        'Liquid Phase: Boundary Conditions in');
    h.bc.st.headline2 = uicontrol(ST{:},[464 687-75*w 100 17],'parent',...
        h.pl.bcs,'string','Direction');
    % --------------------------------------------------------------------
    h.bc.st.axialLB.d1.temperature = uicontrol(TT{:},h.bc.pl.axialLB.d1(3));
    h.bc.st.axialLB.d2.temperature = uicontrol(TT{:},h.bc.pl.axialLB.d2(3));
    % --------------------------------------------------------------------
    h.bc.st.radialLB.rc.temperature = uicontrol(TT{:},h.bc.pl.radialLB.rc(3),...
        'visible','off');
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%                                                           ____          %
%  __    |    |\  /|   |    |   |        /\    _____   |   |    |   |\  | %
% |__    |    | \/ |   |    |   |       /__\     |     |   |    |   | \ | %
%  __|   |    |    |   |____|   |___   /    \    |     |   |____|   |  \| %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save writings
    % panels
    P = {'units','pix','parent',h.pl.simulation,'pos'};
    % edit texts
    E = {'style','edit','units','pix','pos',[185 21 75 23],'parent'};
    % static texts
    S = {'units','pix','style','text','HorizontalAlignment'};
    
% subpanels, subpanels
    h.sn.pl.convergence = uipanel(P{:},[20 640-80*w 680 75],'title','Convergence Criteria');
    h.sn.pl.initialization = uipanel(P{:},[20 525-80*w 680 105],'title','Initialization');
    h.sn.pl.simulation(1) = uipanel(P{:},[20 440-80*w 680 75],'title','Simulation');
    h.sn.pl.simulation(2) = uipanel('units','pix','parent', ...
        h.sn.pl.simulation(1),'pos',[340 -5 355 85]);
    
% button groups
    h.sn.bg.initialization = uibuttongroup('units','pix','BorderType',...
        'none','pos',[15 14 400 67],'parent',h.sn.pl.initialization);

% radio buttons
    h.sn.rb.standard = uicontrol('style','rad','unit','pix','pos',...
        [10 40 450 20],'string','Fluid at rest. Constant mean temperature.',...
        'value',1,'parent',h.sn.bg.initialization);
    h.sn.rb.previous = uicontrol('style','rad','unit','pix','pos',...
        [10 7 190 20],'string','Use previous simulation.','parent',h.sn.bg.initialization);
    
% edit texts
    h.sn.et.residuals(1) = uicontrol(E{:},h.sn.pl.convergence,'string','1.0e-06');
    h.sn.et.residuals(2) = uicontrol(E{:},h.sn.pl.convergence,...
        'string','1.0e-06','pos',[530 21 75 23]);
    h.sn.et.steps = uicontrol(E{:},h.sn.pl.simulation(1),...
        'string','100','pos',[110 21 75 23]);
    
% pushbuttons
    B = {'parent',h.sn.pl.simulation(2),'units','pix','pos',[25 24+1*w 76 25],'visible','on','string'};
    h.sn.pb.run = uicontrol(B{:},'Run','parent',h.sn.pl.simulation(1),...
        'pos',[224 20 76 25],'ForegroundColor',h.selectedTabColor);
    h.sn.pb.pause = uicontrol(B{:},'Pause','parent',h.sn.pl.simulation(1),...
        'pos',[224 20 76 25],'ForegroundColor',h.unselectedTabColor,'visible','off');
    h.sn.pb.load = uicontrol(B{:},'Load');
    h.sn.pb.save = uicontrol(B{:},'Save','pos',[130 24+1*w 76 25]);
    h.sn.pb.plot = uicontrol(B{:},'Plot','pos',[235 24+1*w 76 25]);
    
% axes
    h.sn.as.res = axes(AS{:},h.sn.pl.convergence,'pos',h.sn.et.residuals(1).Position);
    h.sn.as.newton = axes(AS{:},h.sn.pl.convergence,'pos',h.sn.et.residuals(2).Position);
    h.sn.as.residuals = axes(AS{:},h.pl.simulation,'pos',[75 83 610 324]);
    if ispc
        set(h.sn.as.residuals,'pos',[110 55 530 285])
    end
    
% static texts
    uicontrol(S{:},'right','pos',[20 23 80 17],'string','Residuals:',...
        'parent',h.sn.pl.convergence);
    uicontrol(S{:},'right','pos',[283 23 200 17],'string',...
        'Newton''s increment:','parent',h.sn.pl.convergence);
    uicontrol(S{:},'right','pos',[5 23 100 17],'string','Iteration Steps',...
        'parent',h.sn.pl.simulation(1));
    h.sn.st.run = uicontrol(S{:},'left','pos',[5 4 100 17],'string',...
        'running ...','parent',h.pl.simulation,'visible','off');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%                                                                         %
%                            optical ray                                  %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save writings
    % panels
    P = {'units','pix','parent',h.pl.ray,'pos'};
    % edit texts
    E = {'style','edit','units','pix','pos',[48 50 70 23],'parent'};
    
% subpanels, subpanels
    h.or.pl.ray = uipanel(P{:},[20 615-80*w 158 100],'title','Entering Position');
    h.or.pl.particle = uipanel(P{:},[186 615-80*w 158 100],'title','Particle''s Position');
    h.or.pl.N = uipanel(P{:},[352 647-80*w 349 68],'title','Index of Refraction');
    h.or.pl.plot(1) = uipanel(P{:},[20 265-40*w 680 350-40*w],'BorderType','none');
    h.or.pl.plot(2) = uipanel(P{:},[200 0 320 265-40*w],'BorderType','none');

% edit texts
    h.or.et.ray = uicontrol(E{:},h.or.pl.ray,'string',num2str(h.r_0));
    h.or.et.particle = uicontrol(E{:},h.or.pl.particle,'pos',[48 10 70 23],'string',num2str(h.z_p));
    h.or.et.N(1) = uicontrol(E{:},h.or.pl.N,'pos',[72 17 70 23],'string',num2str(h.N_coeff(1)));
    h.or.et.N(2) = uicontrol(E{:},h.or.pl.N,'pos',[160 17 70 23],'string',num2str(h.N_coeff(2)));
    h.or.et.N(3) = uicontrol(E{:},h.or.pl.N,'pos',[273 17 32 23],'string',num2str(h.N_coeff(3)));
    
% pushbuttons
    h.or.pb.run = uicontrol('parent',h.pl.ray,'units','pix',...
        'pos',[624 615-80*w 76 25],'ForegroundColor',h.selectedTabColor,'string','Run');
    
% axes
    h.or.as.ray(1) = axes(AS{:},h.or.pl.ray,'pos',h.or.et.ray.Position);
    h.or.as.ray(2) = axes(AS{:},h.or.pl.ray,'pos',h.or.et.particle.Position);
    h.or.as.particle(1) = axes(AS{:},h.or.pl.particle,'pos',h.or.et.ray.Position);
    h.or.as.particle(2) = axes(AS{:},h.or.pl.particle,'pos',h.or.et.particle.Position);
    h.or.as.N(1) = axes(AS{:},h.or.pl.N,'pos',h.or.et.N(1).Position);
    h.or.as.N(2) = axes(AS{:},h.or.pl.N,'pos',h.or.et.N(2).Position);
    h.or.as.N(3) = axes(AS{:},h.or.pl.N,'pos',h.or.et.N(3).Position);
    h.or.as.plot(1) = axes(AS{:},h.or.pl.plot(1),'OuterPosition',[0 0 680 350-40*w]);
    h.or.as.plot(2) = axes(AS{:},h.or.pl.plot(2),'OuterPosition',[0 2 320 263-40*w]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%             _____          ____                _____                    %
%         __    |     /\    |    |  |  |      |    |     \  /             %
%        |__    |    /__\   |  __|  |  |      |    |      \/              %
%         __|   |   /    \  |    |  |  |      |    |      |               %
%               |  /      \ |____|  |  |____  |    |      |               %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save writings
    % panels
    P = {'units','pix','parent',h.pl.lsa,'pos'};
    % edit texts
    E = {'style','edit','units','pix','pos',[125 70 60 23],'parent'};
    % static texts
    S = {'units','pix','style','text','HorizontalAlignment'};
    
% subpanels, subpanels
    h.sy.pl.n_eig = uipanel(P{:},[20 585-80*w 220 130],'title','Number of Eigenvalues');
    h.sy.pl.m = uipanel(P{:},[250 585-80*w 220 130],'title','Wavenumbers');
    h.sy.pl.parameter = uipanel(P{:},[480 585-80*w 220 130],'title','Parameter of Variation');
    h.sy.pl.run(1) = uipanel(P{:},[20 500-80*w 680 75],'title','Linear Stability Analysis');
    h.sy.pl.run(2) = uipanel('units','pix','parent', ...
        h.sy.pl.run(1),'pos',[340 -5 355 85]);
    
% button groups
    h.sy.bg.parameter = uibuttongroup('units','pix','BorderType',...
        'none','pos',[20 10 30 100],'parent',h.sy.pl.parameter);
    
% radio buttons
    h.sy.rb.delta_T = uicontrol('style','rad','unit','pix','pos',...
        [10 69 15 20],'string',' ','value',1,'parent',h.sy.bg.parameter);
    h.sy.rb.T_d1 = uicontrol('style','rad','unit','pix','pos',...
        [10 41 15 20],'string','  ','parent',h.sy.bg.parameter);
    h.sy.rb.T_d2 = uicontrol('style','rad','unit','pix','pos',...
        [10 13 15 20],'string','   ','parent',h.sy.bg.parameter);
    
% edit texts
    h.sy.et.n_eig = uicontrol(E{:},h.sy.pl.n_eig,'string','12');
    h.sy.et.n_eig_c = uicontrol(E{:},h.sy.pl.n_eig,...
        'string','5','pos',[125 25 60 23]);
    h.sy.et.m_start = uicontrol(E{:},h.sy.pl.m,...
        'string',num2str(h.m_start),'pos',[125 80 60 23]);
    h.sy.et.m_delta = uicontrol(E{:},h.sy.pl.m,...
        'string','1','pos',[125 45 60 23],'enable','off');
    h.sy.et.m_end = uicontrol(E{:},h.sy.pl.m,...
        'string',h.m_end,'pos',[125 10 60 23]);
    
% pushbuttons
    B = {'parent',h.sy.pl.run(2),'units','pix','pos',[25 24+1*w 76 25],'visible','on','string'};
    h.sy.pb.run(1) = uicontrol(B{:},'Most dangerous mode','parent',h.sy.pl.run(1),...
        'pos',[27+11*w 20 160-30*w 25],'ForegroundColor',h.unselectedTabColor);
    h.sy.pb.run(2) = uicontrol(B{:},'Critical mode','parent',h.sy.pl.run(1),...
        'pos',[213-3*w 20 100-16*w 25],'ForegroundColor',h.selectedTabColor);
    h.sy.pb.stop = uicontrol(B{:},'Stop','parent',h.sy.pl.run(1),...
        'pos',h.sy.pb.run(2).Position,'ForegroundColor',h.selectedTabColor,...
        'Visible','off');
    h.sy.pb.load = uicontrol(B{:},'Load');
    h.sy.pb.save = uicontrol(B{:},'Save','pos',[130 24+1*w 76 25]);
    h.sy.pb.plot = uicontrol(B{:},'Plot','pos',[235 24+1*w 76 25]);
    
% axes
    h.sy.as.n = axes(AS{:},h.sy.pl.n_eig,'pos',h.sy.et.n_eig.Position);
    h.sy.as.n_cay = axes(AS{:},h.sy.pl.n_eig,'pos',h.sy.et.n_eig_c.Position);
    h.sy.as.m_start = axes(AS{:},h.sy.pl.m,'pos',h.sy.et.m_start.Position);
    h.sy.as.m_delta = axes(AS{:},h.sy.pl.m,'pos',h.sy.et.m_delta.Position);
    h.sy.as.m_end = axes(AS{:},h.sy.pl.m,'pos',h.sy.et.m_end.Position);
    h.sy.as.delta_T = axes(AS{:},h.sy.pl.parameter,'pos',[50 86 85 20]);
    h.sy.as.T_d1 = axes(AS{:},h.sy.pl.parameter,'pos',[50 58 85 20]);
    h.sy.as.T_d2 = axes(AS{:},h.sy.pl.parameter,'pos',[50 30 85 20]);
    h.sy.as.eigs = axes(AS{:},h.pl.lsa,'pos',[75 83-20*w 610 354-60*w]);
    h.sy.as.dependent = axes(AS{:},h.pl.lsa,'pos',[75 460-80*w 60 23]);
    
% static texts
    h.sy.st.run = uicontrol(S{:},'left','pos',[5 4 100 17],'string',...
        'running ...','parent',h.pl.lsa,'visible','off');
    h.sy.st.m = uicontrol(S{:},'left','pos',[20 475-80*w 100 17],'string',...
        'm = 1 ...','parent',h.pl.lsa,'visible','off');
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%                          ___     ___            _____                   %
%                  /\     |   |   |   |   |   |     |                     %
%                 /__\    | __|   |   |   |   |     |                     %
%                /    \   |___|   |___|   |___|     |                     %
%                                                                         %
%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%_%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axes
    axes(AS{:},h.pl.about,'pos',[120 285-1*w 120 80]);
    imshow(imread('license.png'));
    
% static texts
    uicontrol(S{:},'left','pos',[15 680-80*w 680 24-2*w],'parent',h.pl.about',...
        'FontSize',14,'FontName',font_name,'string','History','FontWeight','Bold','ForegroundColor',h.unselectedTabColor);
    uicontrol(S{:},'left','pos',[30 627-68*w 660 48-12*w],'parent',h.pl.about',...
        'FontSize',12,'FontName',font_name,'string',...
        'MaranStable was developed at TU Wien in the Institute of Fluid Mechanics and Heat Transfer under the guidance of Professor Hendrik C. Kuhlmann.');
    uicontrol(S{:},'left','pos',[30 574-56*w 660 48-12*w],'parent',h.pl.about',...
        'FontSize',12,'FontName',font_name,'string',...
        '07/2009 - 12/2010: Michael Lukasser developed MaranStable 1.0 in the framework of the project Engineering Marangoni Flows (EMA), which was supported by FFG.');
    uicontrol(S{:},'left','pos',[30 497-38*w 660 72-18*w],'parent',h.pl.about',...
        'FontSize',12,'FontName',font_name,'string',...
        '03/2011 - 01/2012: Michael Lukasser continued working on liquid bridges in the framework of the project Technical assistance for the definition and preparation of the JEREMI project on the ISS, which was supported by ESA.');
    uicontrol(S{:},'left','pos',[30 372-8*w 660 120-30*w],'parent',h.pl.about',...
        'FontSize',12,'FontName',font_name,'string',...
        '10/2018 - 04/2022: Mario Stojanovic took over in the framework of the project SAJE, which is the acronym of the FFG Project Stability Analysis for the JEREMI Experiment. He developed MaranStable 2.x and MaranStable 3.x with the help of Francesco Romano, where they removed several bugs, extended the source code by many features, and created a graphical user interface (GUI) for MaranStable.');
    
    uicontrol(S{:},'left','pos',[15 313-1*w 90 24-2*w],'parent',h.pl.about',...
        'FontSize',14,'FontName',font_name,'string','Licensing','FontWeight','Bold','ForegroundColor',h.unselectedTabColor);
    uicontrol(S{:},'left','pos',[30 212+23*w 660 96-24*w],'parent',h.pl.about',...
        'FontSize',12,'FontName',font_name,'string',...
        'This program is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. It comes as is and WITHOUT ANY WARRANTY. By using the program you agree to the above license conditions. For the full text of the license, please visit: https://creativecommons.org/licenses/by-nc-sa/4.0.');
    
    uicontrol(S{:},'left','pos',[15 153+30*w 680 24-2*w],'parent',h.pl.about',...
        'FontSize',14,'FontName',font_name,'string','Contact','FontWeight','Bold','ForegroundColor',h.unselectedTabColor);
    uicontrol(S{:},'left','pos',[30 28+60*w 660 120-30*w],'parent',h.pl.about',...
        'FontSize',12,'FontName',font_name,'string',...
        {'Address: Getreidemarkt 9-BA, 1060 Vienna, Austria.','Homepage: http://www.fluid.tuwien.ac.at','Download pages: https://github.com/fromano88, http://www.fluid.tuwien.ac.at/SAJE','E-mail: maranstable@tuwien.ac.at','Help us to improve MaranStable by reporting bugs and/or by giving us feedback.'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all edit texts summarized    
h.editTextBoxes = {h.gl.et.g h.gl.et.Vr h.gy.et.llb h.gy.et.rc ...
    h.gy.et.ri h.mh.et.radialLB.spacing(1) h.mh.et.radialLB.spacing(2) ...
    h.mh.et.radialLB.spacing(3) h.mh.et.axialLB.spacing(1) ...
    h.mh.et.axialLB.spacing(2) h.mh.et.axialLB.spacing(3) ...
    h.mh.et.radialLB.stretchingC(1) h.mh.et.radialLB.stretchingC(3) ...
    h.mh.et.axialLB.stretchingC(1) h.mh.et.axialLB.stretchingC(3) ...
    h.bc.et.axialLB.d1.temperature h.bc.et.axialLB.d2.temperature ...
    h.bc.et.radialLB.rc.temperature h.bc.et.radialLB.ri.Bi h.sn.et.residuals(1) ...
    h.sn.et.residuals(2) h.sn.et.steps h.or.et.ray h.or.et.particle ...
    h.or.et.N(1) h.or.et.N(2) h.or.et.N(3)...
    h.sy.et.n_eig h.sy.et.n_eig_c h.sy.et.m_start h.sy.et.m_delta h.sy.et.m_end};