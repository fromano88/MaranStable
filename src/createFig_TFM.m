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
    h.pl.sketch = uipanel('unit','pix','pos',[726 4 566 727],...
        'BackgroundColor',[1 1 1]);
    h.as.sketch = axes(AS{:},h.pl.sketch,'pos',[19 8 529 709]);
    h.as.logo   = axes('unit','pix','visible','off','pos',[1250 733 40 25]);
    if ispc
        set(h.pl.sketch,'pos',[726 4 536 647])
        set(h.as.sketch,'pos',[19 8 499 629])
        set(h.as.logo,'pos',[1220 653 40 25])
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
    h.pl.info(1) = uipanel(CP{:},'pos',[890 350 401 410],...
        'BackgroundColor',[1 1 1]);
    h.pl.info(2) = uipanel('parent',h.pl.info(1),'units','pix',...
        'BackgroundColor',h.selectedTabColor,'pos',[0 379 401 48]);
    if ispc
        set(h.pl.info(1),'pos',[890 230 370 420])
        set(h.pl.info(2),'pos',[0 385 368 34])
    end
    uicontrol('parent',h.pl.info(2),'unit','pix','style','text','pos',...
        [1 6+2*w 398-15*w 22],'BackgroundColor',h.selectedTabColor,'ForegroundColor',[1 1 1],...
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
    uicontrol('parent',h.pl.info(1),SA{:},[15 15 100-15*w 17],'string','CC-BY-NC-SA');
    axes(AS{:},h.pl.info(1),'pos',[80 328+2*w 30 30]);
        imshow(imread('logo_TU.png'));
    axes(AS{:},h.pl.info(1),'pos',[292-20*w 262+1*w 45 45]);
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
        'unit','pix','pos',[16 90 260 20],'string','Oberbeck-Boussinesq Approximation');
    h.gl.rb.NS2 = uicontrol(h.gl.bg.model,'style','rad','unit','pix',...
        'pos',[16 55 260 20],'string','Linearly Temperature Dependent');
    h.gl.rb.NS3 = uicontrol(h.gl.bg.model,'style','rad','unit','pix',...
        'pos',[16 20 260 20],'string','Fully Temperature Dependent');
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
        'pix','pos',[62 27 60 23],'string',num2str(h.flowopt.g/9.81));
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
    % --------------------------------------------------------------------
    str = {'Air','Ar','Ne','He','02','H2','Others ...'};
    h.gl.pm.gas = uicontrol('parent',h.gl.pl.fluidProperties,'style',...
        'pop','units','pix','pos',[35 56 130 24],'string',str);
    if max(contains(str,h.selectedGas))
        idx = find(ismember(str,h.selectedGas));
        set(h.gl.pm.gas,'Value',idx)
    else
        str=[str(1:end-1), h.selectedGas, str(end)];
        set(h.gl.pm.gas,'String',str,'value',length(str)-1);
    end
    
% axes (for latex texts)
    h.gl.as.g = axes(AS{:},h.gl.pl.g,'pos',h.gl.et.g.Position);
    h.gl.as.Vr = axes(AS{:},h.gl.pl.Vr,'pos',h.gl.et.Vr.Position);
    
% static texts
    h.gl.st.liquid = uicontrol(ST{:},[26 155 100 17],'parent',...
        h.gl.pl.fluidProperties,'string','Select Liquid');
    h.gl.st.gas = uicontrol(ST{:},[26 85 100 17],'parent', ...
        h.gl.pl.fluidProperties,'string','Select Gas');
    
    
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
        'pos',[20 530-80*w 680 170],'title','Heights');
    h.gy.pl.radii = uipanel('parent',h.pl.geometry,'units','pix','pos',...
        [20 330-80*w 680 170],'title','Radii');
    
% edit texts
    h.gy.et.ld1 = uicontrol('parent',h.gy.pl.lengths,'style','edit',...
        'units','pix','pos',[90 110 60 23],'string',num2str(h.l_d1*1000));
    h.gy.et.llb = uicontrol('parent',h.gy.pl.lengths,'style','edit',...
        'units','pix','pos',[90 67 60 23],'string',num2str(h.l_lb*1000));
    h.gy.et.ld2 = uicontrol('parent',h.gy.pl.lengths,'style','edit',...
        'units','pix','pos',[90 24 60 23],'string',num2str(h.l_d2*1000));
    h.gy.et.rc = uicontrol('parent',h.gy.pl.radii,'style','edit',...
        'units','pix','pos',[90 110 60 23],'string',num2str(h.r_c*1000));
    h.gy.et.ri = uicontrol('parent',h.gy.pl.radii,'style','edit',...
        'units','pix','pos',[90 67 60 23],'string',num2str(h.r_i*1000));
    h.gy.et.ro = uicontrol('parent',h.gy.pl.radii,'style','edit',...
        'units','pix','pos',[90 24 60 23],'string',num2str(h.r_o*1000));
    
% axes (for latex texts)
    h.gy.as.ld1 = axes(AS{:},h.gy.pl.lengths,'pos',h.gy.et.ld1.Position);
    h.gy.as.llb = axes(AS{:},h.gy.pl.lengths,'pos',h.gy.et.llb.Position);
    h.gy.as.ld2 = axes(AS{:},h.gy.pl.lengths,'pos',h.gy.et.ld2.Position);
    h.gy.as.rc = axes(AS{:},h.gy.pl.radii,'pos',h.gy.et.rc.Position);
    h.gy.as.ri = axes(AS{:},h.gy.pl.radii,'pos',h.gy.et.ri.Position);
    h.gy.as.ro = axes(AS{:},h.gy.pl.radii,'pos',h.gy.et.ro.Position);
    
    
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
    % subpanels: coefficient, axialSG
    SC = {'units','pix','pos',[234 -5 233 150],'parent'};
    SA = {'units','pix','pos',[350 -5 380 170],'parent'};
    % button (groups): start, fit, end, axialSG
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
    % radio (button) axialSG: tangent hyp., geometric progression, linear
    RAT = {'style','rad','unit','pix','pos',[5 2 150 20],'string',...
        'Hyperbolic Tangent'};
    RAG = {'style','rad','unit','pix','pos',[180 2 170 20],'string',...
        'Geometric Progression'};
    RAL = {'style','rad','unit','pix','pos',[380 2 60 20],'string','Linear'};
    % start - fit - end
    SFT = {'start', 'fit', 'end'};
    
% subpanels
    h.mh.pl.radialLB.main = uipanel(P{:},h.pl.mesh,'visible','on');
    h.mh.pl.axialLB.main  = uipanel(P{:},h.pl.mesh);
    h.mh.pl.radialSG.main = uipanel(P{:},h.pl.mesh);
    h.mh.pl.axialSG.main  = uipanel(P{:},h.pl.mesh,'pos',[0 50-65*w 720 625]);
    
% subsubpanels
    h.mh.pl.radialLB.spacing   = uipanel(PS{:},h.mh.pl.radialLB.main);
    h.mh.pl.axialLB.spacing    = uipanel(PS{:},h.mh.pl.axialLB.main);
    h.mh.pl.radialSG.spacing   = uipanel(PS{:},h.mh.pl.radialSG.main);
    h.mh.pl.axialSG.spacing(1) = uipanel(PS{:},h.mh.pl.axialSG.main,'pos',[10 464 700 158]);
    % --------------------------------------------------------------------
    h.mh.pl.radialLB.stretchingC    = uipanel(PC{:},h.mh.pl.radialLB.main);
    h.mh.pl.axialLB.stretchingC     = uipanel(PC{:},h.mh.pl.axialLB.main);
    h.mh.pl.radialSG.stretchingC    = uipanel(PC{:},h.mh.pl.radialSG.main);
    h.mh.pl.axialSG.stretchingC(1)  = uipanel(PC{:},h.mh.pl.axialSG.main,'pos',[10 301 700 158]);
    % --------------------------------------------------------------------
    h.mh.pl.radialLB.stretchingF(1) = uipanel(PF{:},h.mh.pl.radialLB.main);
    h.mh.pl.axialLB.stretchingF(1)  = uipanel(PF{:},h.mh.pl.axialLB.main);
    h.mh.pl.radialSG.stretchingF(1) = uipanel(PF{:},h.mh.pl.radialSG.main);
    h.mh.pl.axialSG.stretchingF     = uipanel(PF{:},h.mh.pl.axialSG.main,'pos',[10 28 700 268]);
    
% subsubsubpanels
    h.mh.pl.radialLB.stretchingF(2) = uipanel(SC{:},h.mh.pl.radialLB.stretchingF(1));
    h.mh.pl.axialLB.stretchingF(2)  = uipanel(SC{:},h.mh.pl.axialLB.stretchingF(1));
    h.mh.pl.radialSG.stretchingF(2) = uipanel(SC{:},h.mh.pl.radialSG.stretchingF(1));
    h.mh.pl.axialSG.spacing(2)      = uipanel(SA{:},h.mh.pl.axialSG.spacing(1));
    h.mh.pl.axialSG.stretchingC(2)  = uipanel(SA{:},h.mh.pl.axialSG.stretchingC(1));
    
% button groups
    for i = 1:3
        try
            h.mh.bg.radialLB(i) = uibuttongroup(BG{i,:},h.mh.pl.radialLB.stretchingF(i));
            h.mh.bg.axialLB(i) = uibuttongroup(BG{i,:},h.mh.pl.axialLB.stretchingF(i));
            h.mh.bg.radialSG(i) = uibuttongroup(BG{i,:},h.mh.pl.radialSG.stretchingF(i));
        catch
            h.mh.bg.radialLB(i) = uibuttongroup(BG{i,:},h.mh.pl.radialLB.stretchingF(i-2));
            h.mh.bg.axialLB(i) = uibuttongroup(BG{i,:},h.mh.pl.axialLB.stretchingF(i-2));
            h.mh.bg.radialSG(i) = uibuttongroup(BG{i,:},h.mh.pl.radialSG.stretchingF(i-2));
        end
    end
    % --------------------------------------------------------------------
    BA = {h.mh.pl.axialSG.stretchingF,'units','pix','BorderType','none','pos'};
    for i = 1:7
        h.mh.bg.axialSG(i) = uibuttongroup(BA{:},[220 (260-35*i) 460 25]);
    end
    h.mh.buttonGroups.radial = {h.mh.bg.radialLB(1) h.mh.bg.radialLB(3) ...
        h.mh.bg.radialSG(1) h.mh.bg.radialSG(3)};
    
    h.mh.buttonGroups.axial = {h.mh.bg.axialLB(1) h.mh.bg.axialLB(3)...
        h.mh.bg.axialSG(1) h.mh.bg.axialSG(3) h.mh.bg.axialSG(5) ...
        h.mh.bg.axialSG(7)};
    
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
        h.mh.rb.radialSG(i).tanh = uicontrol(h.mh.bg.radialSG(i),RT{:});
        h.mh.rb.radialSG(i).gp = uicontrol(h.mh.bg.radialSG(i),RG{:});
        h.mh.rb.radialSG(i).lin = uicontrol(h.mh.bg.radialSG(i),RL{:});
        % ----------------------------------------------------------------
        if strcmp(h.mesh.r.sf,'tanh')
            set([h.mh.rb.radialLB(i).tanh h.mh.rb.radialSG(i).tanh],'value',1)
        else
            set([h.mh.rb.radialLB(i).gp h.mh.rb.radialSG(i).gp],'value',1)
        end
        if strcmp(h.mesh.z.sf,'tanh')
            set(h.mh.rb.axialLB(i).tanh,'value',1)
        else
            set(h.mh.rb.axialLB(i).gp,'value',1)
        end
    end
    % --------------------------------------------------------------------
    set([h.mh.rb.radialLB(2).tanh h.mh.rb.axialLB(2).tanh h.mh.rb.radialSG(2).tanh],'enable','off','value',0)
    set([h.mh.rb.radialLB(2).gp h.mh.rb.axialLB(2).gp h.mh.rb.radialSG(2).gp],'enable','off')
    set([h.mh.rb.radialLB(2).lin h.mh.rb.axialLB(2).lin h.mh.rb.radialSG(2).lin],'value',1)
    % --------------------------------------------------------------------
    for i = 1:7
        h.mh.rb.axialSG(i).tanh = uicontrol(h.mh.bg.axialSG(i),RAT{:});
        h.mh.rb.axialSG(i).gp = uicontrol(h.mh.bg.axialSG(i),RAG{:});
        h.mh.rb.axialSG(i).lin = uicontrol(h.mh.bg.axialSG(i),RAL{:});
        if strcmp(h.mesh.z.sf,'tanh')
            set(h.mh.rb.axialSG(i).tanh,'value',1)
        else
            set(h.mh.rb.axialSG(i).gp,'value',1)
        end
        if mod(i,2) == 0
            set(h.mh.rb.axialSG(i).tanh,'enable','off','value',0)
            set(h.mh.rb.axialSG(i).gp,'enable','off')
            set(h.mh.rb.axialSG(i).lin,'value',1)
        end
    end
    
% edit texts
    for i = 1:3
        h.mh.et.radialLB.spacing(i) = uicontrol('pos',...
            [325 (120-35*i) 60 23],ET{:},h.mh.pl.radialLB.spacing,...
            'string',eval(['num2str(h.mesh.r.delta.' SFT{i} '(1))']));
        h.mh.et.axialLB.spacing(i) = uicontrol('pos',...
            [325 (120-35*i) 60 23],ET{:},h.mh.pl.axialLB.spacing,...
            'string',eval(['num2str(h.mesh.z.delta.' SFT{i} '(2))']));
        h.mh.et.radialSG.spacing(i) = uicontrol('pos',...
            [325 (120-35*i) 60 23],ET{:},h.mh.pl.radialSG.spacing,...
            'string',eval(['num2str(h.mesh.r.delta.' SFT{i} '(2))']));
        % ----------------------------------------------------------------
        h.mh.et.axialSG.spacing(i) = uicontrol('pos',...
            [220 (150-35*i) 60 23],ET{:},h.mh.pl.axialSG.spacing(1),...
            'string',eval(['num2str(h.mesh.z.delta.' SFT{i} '(1))']));
    end
    for i = 6:7
        h.mh.et.axialSG.spacing(i) = uicontrol('pos',...
            [215 (154-35*(i-4)) 60 23],ET{:},h.mh.pl.axialSG.spacing(2),...
            'string',eval(['num2str(h.mesh.z.delta.' SFT{i-4} '(3))']));
    end
    h.mh.et.axialSG.spacing(4) = uicontrol('pos',[220 10 60 23],ET{:},...
        h.mh.pl.axialSG.spacing(1),'string',num2str(h.mesh.z.delta.fit(2)));
    h.mh.et.axialSG.spacing(5) = uicontrol('pos',[215 119 60 23],ET{:},...
        h.mh.pl.axialSG.spacing(2),'string',num2str(h.mesh.z.delta.end(2)));
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
        h.mh.pl.axialLB.stretchingC,'string',num2str(h.mesh.z.f(2)));
    h.mh.et.axialLB.stretchingC(2) = uicontrol('pos',[325 50 60 23],ET{:},...
        h.mh.pl.axialLB.stretchingC,'string',num2str(1),'enable','off');
    h.mh.et.axialLB.stretchingC(3) = uicontrol('pos',[325 15 60 23],ET{:},...
        h.mh.pl.axialLB.stretchingC,'string',num2str(h.mesh.z.f(3)));
    % --------------------------------------------------------------------
    h.mh.et.radialSG.stretchingC(1) = uicontrol('pos',[325 85 60 23],ET{:},...
        h.mh.pl.radialSG.stretchingC,'string',num2str(h.mesh.r.f(3)));
    h.mh.et.radialSG.stretchingC(2) = uicontrol('pos',[325 50 60 23],ET{:},...
        h.mh.pl.radialSG.stretchingC,'string',num2str(1),'enable','off');
    h.mh.et.radialSG.stretchingC(3) = uicontrol('pos',[325 15 60 23],ET{:},...
        h.mh.pl.radialSG.stretchingC,'string',num2str(h.mesh.r.f(4)));
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    h.mh.et.axialSG.stretchingC(1) = uicontrol('pos',[220 115 60 23],ET{:},...
        h.mh.pl.axialSG.stretchingC(1),'string',num2str(h.mesh.z.f(1)));
    h.mh.et.axialSG.stretchingC(2) = uicontrol('pos',[220 80 60 23],ET{:},...
        h.mh.pl.axialSG.stretchingC(1),'string',num2str(1),'enable','off');
    h.mh.et.axialSG.stretchingC(3) = uicontrol('pos',[220 45 60 23],ET{:},...
        h.mh.pl.axialSG.stretchingC(1),'string',num2str(h.mesh.z.f(2)));
    h.mh.et.axialSG.stretchingC(4) = uicontrol('pos',[220 10 60 23],ET{:},...
        h.mh.pl.axialSG.stretchingC(1),'string',num2str(1),'enable','off');
    h.mh.et.axialSG.stretchingC(5) = uicontrol('pos',[215 119 60 23],ET{:},...
        h.mh.pl.axialSG.stretchingC(2),'string',num2str(h.mesh.z.f(3)));
    h.mh.et.axialSG.stretchingC(6) = uicontrol('pos',[215 84 60 23],ET{:},...
        h.mh.pl.axialSG.stretchingC(2),'string',num2str(1),'enable','off');
    h.mh.et.axialSG.stretchingC(7) = uicontrol('pos',[215 49 60 23],ET{:},...
        h.mh.pl.axialSG.stretchingC(2),'string',num2str(h.mesh.z.f(4)));
    
    
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
        h.mh.as.radialSG.spacing(i) = axes(AS{:},...
            h.mh.pl.radialSG.spacing,'pos',h.mh.et.radialSG.spacing(i).Position);
        % ----------------------------------------------------------------
        h.mh.as.radialLB.stretchingC(i) = axes(AS{:},...
            h.mh.pl.radialLB.stretchingC,'pos',h.mh.et.radialLB.stretchingC(i).Position);
        h.mh.as.axialLB.stretchingC(i) = axes(AS{:},...
            h.mh.pl.axialLB.stretchingC,'pos',h.mh.et.axialLB.stretchingC(i).Position);
        h.mh.as.radialSG.stretchingC(i) = axes(AS{:},...
            h.mh.pl.radialSG.stretchingC,'pos',h.mh.et.radialSG.stretchingC(i).Position);
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
    % --------------------------------------------------------------------
    h.mh.as.radialSG.stretchingF(1) = axes(AS{:},...
        h.mh.pl.radialSG.stretchingF(1),'pos',[12 92 20 17]);
    h.mh.as.radialSG.stretchingF(2) = axes(AS{:},...
        h.mh.pl.radialSG.stretchingF(2),'pos',[12 96 20 17]);
    h.mh.as.radialSG.stretchingF(3) = axes(AS{:},...
        h.mh.pl.radialSG.stretchingF(1),'pos',[477 92 20 17]);
    for i = 1:4
        h.mh.as.axialSG.spacing(i) = axes(AS{:},...
            h.mh.pl.axialSG.spacing(1),'pos',h.mh.et.axialSG.spacing(i).Position);
        h.mh.as.axialSG.stretchingC(i) = axes(AS{:},...
            h.mh.pl.axialSG.stretchingC(1),'pos',h.mh.et.axialSG.stretchingC(i).Position);
    end
    for i = 5:7
        h.mh.as.axialSG.spacing(i) = axes(AS{:},...
            h.mh.pl.axialSG.spacing(2),'pos',h.mh.et.axialSG.spacing(i).Position);
        h.mh.as.axialSG.stretchingC(i) = axes(AS{:},...
            h.mh.pl.axialSG.stretchingC(2),'pos',h.mh.et.axialSG.stretchingC(i).Position);
    end
    for i = 1:7
        h.mh.as.axialSG.stretchingF(i) = axes(AS{:},...
            h.mh.pl.axialSG.stretchingF,'pos',[10 (261-35*i) 20 17]);
    end
    
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
    % button (groups): velocity, heat
    BV = {'units','pix','BorderType','none','pos',[15 70-10*w 195 80],'parent'};
    BH = {'units','pix','BorderType','none','pos',[467 66-10*w 195 80],'parent'};
    % radio (button): wall, noSlip, slip, adiabatic, conductive
    RW = {'style','rad','unit','pix','pos',[25 95-10*w 140 20],'string',...
        'No Penetration','value',1,'enable','inactive','parent'};
    RN = {'style','rad','unit','pix','pos',[10 50 110 20],'string',...
        'No Slip','value',1};
    RS = {'style','rad','unit','pix','pos',[10 10 140 20],'string',...
        'Slip','value',0};
    RA = {'style','rad','unit','pix','pos',[10 50 110 20],'string','Adiabatic'};
    RC = {'style','rad','unit','pix','pos',[10 10 110 20],'string',...
        'Conductive','value',1};
    % edit texts
    E = {'style','edit','units','pix','pos',[90 18 142 23],'parent'};
    % popup menu
    PM = {'style','pop','units','pix','pos',[540 17 88 24],...
        'string',{'constant','variable'},'parent'};
    % (statix) texts: at, temperature,
    TT = {'units','pix','style','text','HorizontalAlignment','left',...
            'pos',[400+25*w 20 135-25*w 17],'string','Temperature Profile:','parent'};
    
% subpanels
    h.bc.pl.radialLB.main = uipanel(P{:},h.pl.bcs,'visible','on');
    h.bc.pl.axialLB.main  = uipanel(P{:},h.pl.bcs);
    h.bc.pl.radialSG.main = uipanel(P{:},h.pl.bcs);
    h.bc.pl.axialSG.main  = uipanel(P{:},h.pl.bcs);
    
% subsubpanels
    h.bc.pl.axialLB.d1(1) = uipanel(SP1{:},h.bc.pl.axialLB.main);
    h.bc.pl.axialLB.d2(1) = uipanel(SP2{:},h.bc.pl.axialLB.main);
    h.bc.pl.radialLB.rc(1) = uipanel(SP1{:},h.bc.pl.radialLB.main);
    h.bc.pl.radialLB.ri(1) = uipanel(SP2{:},h.bc.pl.radialLB.main);
    h.bc.pl.axialSG.d1(1) = uipanel(SP1{:},h.bc.pl.axialSG.main,'pos',[20 430-60*w 680 210]);
    h.bc.pl.axialSG.d2(1) = uipanel(SP1{:},h.bc.pl.axialSG.main,'pos',[20 160-60*w 680 210]);
    h.bc.pl.radialSG.d1(1) = uipanel(SP1{:},h.bc.pl.radialSG.main);
    h.bc.pl.radialSG.d2(1) = uipanel(SP2{:},h.bc.pl.radialSG.main);
    h.bc.pl.radialSG.ro(1) = uipanel(SP1{:},h.bc.pl.radialSG.main,'pos',[20 40 680 160-20*w]);
    
% subsubsubpanels
    h.bc.pl.axialLB.d1(2) = uipanel(SSP1{:},h.bc.pl.axialLB.d1(1));
    h.bc.pl.axialLB.d1(3) = uipanel(SSP2{:},h.bc.pl.axialLB.d1(1));
    h.bc.pl.axialLB.d2(2) = uipanel(SSP1{:},h.bc.pl.axialLB.d2(1));
    h.bc.pl.axialLB.d2(3) = uipanel(SSP2{:},h.bc.pl.axialLB.d2(1));
    % --------------------------------------------------------------------
    h.bc.pl.radialLB.rc(2) = uipanel(SSP1{:},h.bc.pl.radialLB.rc(1));
    h.bc.pl.radialLB.rc(3) = uipanel(SSP2{:},h.bc.pl.radialLB.rc(1));
    h.bc.pl.radialLB.ri(2) = uipanel(SSP1{:},h.bc.pl.radialLB.ri(1),'pos',[451 -5 232 170]);
    % --------------------------------------------------------------------
    h.bc.pl.axialSG.d1(2) = uipanel(SSP1{:},h.bc.pl.axialSG.d1(1),'pos',[227 50 226 170]);
    h.bc.pl.axialSG.d1(3) = uipanel(SSP1{:},h.bc.pl.axialSG.d1(1),'pos',[-5 48 690 52]);
    % --------------------------------------------------------------------
    h.bc.pl.axialSG.d2(2) = uipanel(SSP1{:},h.bc.pl.axialSG.d2(1),'pos',[227 50 226 170]);
    h.bc.pl.axialSG.d2(3) = uipanel(SSP1{:},h.bc.pl.axialSG.d2(1),'pos',[-5 48 690 52]);
    % --------------------------------------------------------------------
    h.bc.pl.radialSG.d1(2) = uipanel(SSP1{:},h.bc.pl.radialSG.d1(1));
    h.bc.pl.radialSG.d1(3) = uipanel(SSP2{:},h.bc.pl.radialSG.d1(1));
    h.bc.pl.radialSG.d2(2) = uipanel(SSP1{:},h.bc.pl.radialSG.d2(1));
    h.bc.pl.radialSG.d2(3) = uipanel(SSP2{:},h.bc.pl.radialSG.d2(1));
    h.bc.pl.radialSG.ro(2) = uipanel(SSP1{:},h.bc.pl.radialSG.ro(1));
    h.bc.pl.radialSG.ro(3) = uipanel(SSP2{:},h.bc.pl.radialSG.ro(1));
    
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
    h.bc.bg.axialSG.d1(1) = uibuttongroup(BV{:},h.bc.pl.axialSG.d1(1),...
        'pos',[15 111 195 90]);
    h.bc.bg.axialSG.d1(2) = uibuttongroup(BV{:},h.bc.pl.axialSG.d1(2),...
        'pos',[15 65 195 80]);
    h.bc.bg.axialSG.d1(3) = uibuttongroup(BH{:},h.bc.pl.axialSG.d1(1),...
        'pos',[467 116 195 80]);
    h.bc.bg.axialSG.d2(1) = uibuttongroup(BV{:},h.bc.pl.axialSG.d2(1),...
        'pos',[15 111 195 90]);
    h.bc.bg.axialSG.d2(2) = uibuttongroup(BV{:},h.bc.pl.axialSG.d2(2),...
        'pos',[15 65 195 80]);
    h.bc.bg.axialSG.d2(3) = uibuttongroup(BH{:},h.bc.pl.axialSG.d2(1),...
        'pos',[467 116 195 80]);
    % --------------------------------------------------------------------
    h.bc.bg.radialSG.d1(2) = uibuttongroup(BV{:},h.bc.pl.radialSG.d1(2));
    h.bc.bg.radialSG.d1(3) = uibuttongroup(BH{:},h.bc.pl.radialSG.d1(1));
    h.bc.bg.radialSG.d2(2) = uibuttongroup(BV{:},h.bc.pl.radialSG.d2(2));
    h.bc.bg.radialSG.d2(3) = uibuttongroup(BH{:},h.bc.pl.radialSG.d2(1));
    h.bc.bg.radialSG.ro(2) = uibuttongroup(BV{:},h.bc.pl.radialSG.ro(2));
    h.bc.bg.radialSG.ro(3) = uibuttongroup(BH{:},h.bc.pl.radialSG.ro(1));
    % --------------------------------------------------------------------
    h.bc.buttonGroups.velocity = {h.bc.bg.axialLB.d1(2) ...
        h.bc.bg.axialLB.d2(2) h.bc.bg.radialLB.rc(2) h.bc.bg.axialSG.d1(1) ...
        h.bc.bg.axialSG.d2(1) h.bc.bg.axialSG.d1(2) h.bc.bg.axialSG.d2(2) ...
        h.bc.bg.radialSG.d1(2) h.bc.bg.radialSG.d2(2) h.bc.bg.radialSG.ro(2)};
    h.bc.buttonGroups.temperature = {h.bc.bg.axialLB.d1(3) ...
        h.bc.bg.axialLB.d2(3) h.bc.bg.radialLB.rc(3) h.bc.bg.radialLB.ri(2) ...
        h.bc.bg.axialSG.d1(3) h.bc.bg.axialSG.d2(3) h.bc.bg.radialSG.d1(3) ...
        h.bc.bg.radialSG.d2(3) h.bc.bg.radialSG.ro(3)};
    
% radio buttons
    uicontrol(RW{:},h.bc.pl.axialLB.d1(1));
    uicontrol(RW{:},h.bc.pl.axialLB.d2(1));
    uicontrol(RW{:},h.bc.pl.radialSG.d1(1));
    uicontrol(RW{:},h.bc.pl.radialSG.d2(1));
    uicontrol(RW{:},h.bc.pl.radialSG.ro(1));
    % --------------------------------------------------------------------
    h.bc.rb.axialLB.d1.slip  = uicontrol(h.bc.bg.axialLB.d1(2),RS{:});
    h.bc.rb.axialLB.d2.slip  = uicontrol(h.bc.bg.axialLB.d2(2),RS{:});
    h.bc.rb.radialLB.rc.slip = uicontrol(h.bc.bg.radialLB.rc(2),RS{:});
    h.bc.rb.axialSG.d1.slip  = uicontrol(h.bc.bg.axialSG.d1(2),RS{:});
    h.bc.rb.axialSG.d2.slip  = uicontrol(h.bc.bg.axialSG.d2(2),RS{:});
    h.bc.rb.radialSG.d1.slip = uicontrol(h.bc.bg.radialSG.d1(2),RS{:});
    h.bc.rb.radialSG.d2.slip = uicontrol(h.bc.bg.radialSG.d2(2),RS{:});
    h.bc.rb.radialSG.ro.slip = uicontrol(h.bc.bg.radialSG.ro(2),RS{:});
    % --------------------------------------------------------------------
    h.bc.rb.axialLB.d1.noSlip  = uicontrol(h.bc.bg.axialLB.d1(2),RN{:});
    h.bc.rb.axialLB.d2.noSlip  = uicontrol(h.bc.bg.axialLB.d2(2),RN{:});
    h.bc.rb.radialLB.rc.noSlip = uicontrol(h.bc.bg.radialLB.rc(2),RN{:});
    h.bc.rb.axialSG.d1.noSlip  = uicontrol(h.bc.bg.axialSG.d1(2),RN{:});
    h.bc.rb.axialSG.d2.noSlip  = uicontrol(h.bc.bg.axialSG.d2(2),RN{:});
    h.bc.rb.radialSG.d1.noSlip = uicontrol(h.bc.bg.radialSG.d1(2),RN{:});
    h.bc.rb.radialSG.d2.noSlip = uicontrol(h.bc.bg.radialSG.d2(2),RN{:});
    h.bc.rb.radialSG.ro.noSlip = uicontrol(h.bc.bg.radialSG.ro(2),RN{:});
    set([h.bc.rb.axialSG.d1.noSlip h.bc.rb.axialSG.d1.slip h.bc.rb.axialSG.d2.noSlip h.bc.rb.axialSG.d2.slip],'value',0,'enable','off')
    % --------------------------------------------------------------------
    h.bc.rb.axialLB.d1.adiabatic  = uicontrol(h.bc.bg.axialLB.d1(3),RA{:});
    h.bc.rb.axialLB.d1.conductive = uicontrol(h.bc.bg.axialLB.d1(3),RC{:});
    h.bc.rb.axialLB.d2.adiabatic  = uicontrol(h.bc.bg.axialLB.d2(3),RA{:});
    h.bc.rb.axialLB.d2.conductive = uicontrol(h.bc.bg.axialLB.d2(3),RC{:});
    % --------------------------------------------------------------------
    h.bc.rb.radialLB.rc.wall       = uicontrol(h.bc.bg.radialLB.rc(1),...
        RS{:},'string','No Penetration','enable','on');
    h.bc.rb.radialLB.rc.axis       = uicontrol(h.bc.bg.radialLB.rc(1),...
        RN{:},'string','Symmetry','enable','on');
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
    h.bc.rb.radialLB.ri.adiabatic  = uicontrol(h.bc.bg.radialLB.ri(2),RA{:});
    h.bc.rb.radialLB.ri.conductive = uicontrol(h.bc.bg.radialLB.ri(2),RC{:});
    set([h.bc.rb.radialLB.rc.noSlip h.bc.rb.radialLB.rc.slip h.bc.rb.radialLB.rc.adiabatic h.bc.rb.radialLB.rc.conductive h.bc.rb.radialLB.ri.adiabatic],'value',0,'enable','off')
    % --------------------------------------------------------------------
    h.bc.rb.axialSG.d1.wall       = uicontrol(h.bc.bg.axialSG.d1(1),...
        RA{:},'string','No Penetration','pos',[10 70 140 20]);
    h.bc.rb.axialSG.d1.inflow     = uicontrol(h.bc.bg.axialSG.d1(1),...
        RA{:},'string','Inflow','pos',[10 35 70 20]);
    h.bc.rb.axialSG.d1.outflow    = uicontrol(h.bc.bg.axialSG.d1(1),...
        RA{:},'string','Outflow','pos',[10 0 70 20],'value',1);
    h.bc.rb.axialSG.d2.wall       = uicontrol(h.bc.bg.axialSG.d2(1),...
        RA{:},'string','No Penetration','pos',[10 70 140 20]);
    h.bc.rb.axialSG.d2.inflow     = uicontrol(h.bc.bg.axialSG.d2(1),...
        RA{:},'string','Inflow','pos',[10 35 70 20],'value',1);
    h.bc.rb.axialSG.d2.outflow    = uicontrol(h.bc.bg.axialSG.d2(1),...
        RA{:},'string','Outflow','pos',[10 0 70 20]);
    h.bc.rb.axialSG.d1.conductive = uicontrol(h.bc.bg.axialSG.d1(3),RC{:});
    h.bc.rb.axialSG.d1.adiabatic  = uicontrol(h.bc.bg.axialSG.d1(3),RA{:});
    h.bc.rb.axialSG.d2.adiabatic  = uicontrol(h.bc.bg.axialSG.d2(3),RA{:});
    h.bc.rb.axialSG.d2.conductive = uicontrol(h.bc.bg.axialSG.d2(3),RC{:});
    % --------------------------------------------------------------------
    h.bc.rb.radialSG.d1.adiabatic  = uicontrol(h.bc.bg.radialSG.d1(3),RA{:});
    h.bc.rb.radialSG.d1.conductive = uicontrol(h.bc.bg.radialSG.d1(3),RC{:});
    h.bc.rb.radialSG.d2.adiabatic  = uicontrol(h.bc.bg.radialSG.d2(3),RA{:});
    h.bc.rb.radialSG.d2.conductive = uicontrol(h.bc.bg.radialSG.d2(3),RC{:});
    h.bc.rb.radialSG.ro.adiabatic  = uicontrol(h.bc.bg.radialSG.ro(3),RA{:});
    h.bc.rb.radialSG.ro.conductive = uicontrol(h.bc.bg.radialSG.ro(3),RC{:});
    
% edit texts
    h.bc.et.axialLB.d1.temperature = uicontrol(E{:},...
        h.bc.pl.axialLB.d1(3),'string',num2str(h.T_d1l));
    h.bc.et.axialLB.d2.temperature = uicontrol(E{:},...
        h.bc.pl.axialLB.d2(3),'string',num2str(h.T_d2l));
    % --------------------------------------------------------------------
    h.bc.et.radialLB.rc.temperature = uicontrol(E{:},...
        h.bc.pl.radialLB.rc(3),'string',[num2str((h.T_d1l+h.T_d2l)*0.5) '+' num2str(h.T_d1l-h.T_d2l) '*z/d'],...
        'visible','off','enable','inactive');
    % --------------------------------------------------------------------
    h.bc.et.axialSG.d1.temperature = uicontrol(E{:},...
        h.bc.pl.axialSG.d1(3),'string',num2str(h.T_p1),'pos',[90 15 142 23]);
    h.bc.et.axialSG.d2.temperature = uicontrol(E{:},...
        h.bc.pl.axialSG.d2(3),'string',num2str(h.T_p2),'pos',[90 15 142 23]);
    h.bc.et.axialSG.d1.velocity = uicontrol(E{:},...
        h.bc.pl.axialSG.d1(1),'string',num2str(h.w_in),...
        'pos',[86 13 142 23],'visible','off');
    h.bc.et.axialSG.d2.velocity = uicontrol(E{:},...
        h.bc.pl.axialSG.d2(1),'string',num2str(-h.w_in),'pos',[86 13 142 23]);
    % --------------------------------------------------------------------
    h.bc.et.radialSG.d1.temperature = uicontrol(E{:},...
        h.bc.pl.radialSG.d1(3),'string',num2str(h.T_d1g));
    h.bc.et.radialSG.d2.temperature = uicontrol(E{:},...
        h.bc.pl.radialSG.d2(3),'string',num2str(h.T_d2g));
    h.bc.et.radialSG.ro.temperature = uicontrol(E{:},...
        h.bc.pl.radialSG.ro(3),'string',num2str((h.T_d1g+h.T_d2g)/2));
    
% popup menus
    h.bc.pm.headline = uicontrol('parent',h.pl.bcs,'style','pop',...
        'pos',[373 683-75*w 92 24],'string',{'inner/outer','top/bottom'});
    % --------------------------------------------------------------------
    h.bc.pm.axialLB.d1.temperature = uicontrol(PM{:},h.bc.pl.axialLB.d1(3));
    h.bc.pm.axialLB.d2.temperature = uicontrol(PM{:},h.bc.pl.axialLB.d2(3));
    % --------------------------------------------------------------------
    h.bc.pm.radialLB.rc.temperature = uicontrol(PM{:},h.bc.pl.radialLB.rc(3),...
        'value',2,'visible','off');
    % --------------------------------------------------------------------
    h.bc.pm.axialSG.d1.temperature = uicontrol(PM{:},...
        h.bc.pl.axialSG.d1(3),'pos',[540 14 88 24]);
    h.bc.pm.axialSG.d2.temperature = uicontrol(PM{:},...
        h.bc.pl.axialSG.d2(3),'pos',[540 14 88 24]);
    h.bc.pm.axialSG.d1.velocity = uicontrol(PM{:},...
        h.bc.pl.axialSG.d1(1),'pos',[536 11 88 24],'visible','off');
    h.bc.pm.axialSG.d2.velocity = uicontrol(PM{:},...
        h.bc.pl.axialSG.d2(1),'pos',[536 11 88 24]);
    % --------------------------------------------------------------------
    h.bc.pm.radialSG.d1.temperature = uicontrol(PM{:},h.bc.pl.radialSG.d1(3));
    h.bc.pm.radialSG.d2.temperature = uicontrol(PM{:},h.bc.pl.radialSG.d2(3));
    h.bc.pm.radialSG.ro.temperature = uicontrol(PM{:},h.bc.pl.radialSG.ro(3));
    % --------------------------------------------------------------------
    h.bc.popupMenus.temperature = {h.bc.pm.axialLB.d1.temperature ...
        h.bc.pm.axialLB.d2.temperature h.bc.pm.radialLB.rc.temperature ...
        h.bc.pm.axialSG.d1.temperature h.bc.pm.axialSG.d2.temperature ...
        h.bc.pm.radialSG.d1.temperature h.bc.pm.radialSG.d2.temperature ...
        h.bc.pm.radialSG.ro.temperature};
    h.bc.popupMenus.velocity = {h.bc.pm.axialSG.d1.velocity h.bc.pm.axialSG.d2.velocity};

% axes
    h.bc.as.axialLB.d1(1)  = axes(AS{:},h.bc.pl.axialLB.main,'pos',[20 640-60*w 20 17]);
    h.bc.as.axialLB.d2(1)  = axes(AS{:},h.bc.pl.axialLB.main,'pos',[20 420-40*w 20 17]);
    h.bc.as.radialLB.rc(1) = axes(AS{:},h.bc.pl.radialLB.main,'pos',[20 640-60*w 20 17]);
    h.bc.as.radialLB.ri = axes(AS{:},h.bc.pl.radialLB.main,'pos',[20 420-40*w 20 17]);
    h.bc.as.axialSG.d1(1)  = axes(AS{:},h.bc.pl.axialSG.main,'pos',[20 640-60*w 20 17]);
    h.bc.as.axialSG.d2(1)  = axes(AS{:},h.bc.pl.axialSG.main,'pos',[20 370-60*w 20 17]);
    h.bc.as.radialSG.d1(1) = axes(AS{:},h.bc.pl.radialSG.main,'pos',[20 640-60*w 20 17]);
    h.bc.as.radialSG.d2(1) = axes(AS{:},h.bc.pl.radialSG.main,'pos',[20 420-40*w 20 17]);
    h.bc.as.radialSG.ro(1) = axes(AS{:},h.bc.pl.radialSG.main,'pos',[20 200-20*w 20 17]);
    % --------------------------------------------------------------------
    h.bc.as.axialLB.d1(2)  = axes(AS{:},h.bc.pl.axialLB.d1(3),'pos',[90 18 142 23]);
    h.bc.as.axialLB.d2(2)  = axes(AS{:},h.bc.pl.axialLB.d2(3),'pos',[90 18 142 23]);
    h.bc.as.radialLB.rc(2)  = axes(AS{:},h.bc.pl.radialLB.rc(3),'pos',[90 18 142 23]);
    h.bc.as.axialSG.d1(2)  = axes(AS{:},h.bc.pl.axialSG.d1(3),'pos',[90 15 142 23]);
    h.bc.as.axialSG.d2(2)  = axes(AS{:},h.bc.pl.axialSG.d2(3),'pos',[90 15 142 23]);
    h.bc.as.radialSG.d1(2)  = axes(AS{:},h.bc.pl.radialSG.d1(3),'pos',[90 18 142 23]);
    h.bc.as.radialSG.d2(2)  = axes(AS{:},h.bc.pl.radialSG.d2(3),'pos',[90 18 142 23]);
    h.bc.as.radialSG.ro(2)  = axes(AS{:},h.bc.pl.radialSG.ro(3),'pos',[90 18 142 23]);
    % --------------------------------------------------------------------
    h.bc.as.axialSG.d1(3)  = axes(AS{:},h.bc.pl.axialSG.d1(1),'pos',[90 15 142 23]);
    h.bc.as.axialSG.d2(3)  = axes(AS{:},h.bc.pl.axialSG.d2(1),'pos',[90 15 142 23]);
    

% static texts
    h.bc.st.headline1 = uicontrol(ST{:},[60 687-75*w 310 17],'parent',...
        h.pl.bcs,'HorizontalAlignment','right','string',...
        'Liquid Phase: Boundary Conditions on');
    h.bc.st.headline2 = uicontrol(ST{:},[469 687-75*w 100 17],'parent',...
        h.pl.bcs,'string','Boundary');
    % --------------------------------------------------------------------
    h.bc.st.axialLB.d1.temperature = uicontrol(TT{:},h.bc.pl.axialLB.d1(3));
    h.bc.st.axialLB.d2.temperature = uicontrol(TT{:},h.bc.pl.axialLB.d2(3));
    % --------------------------------------------------------------------
    h.bc.st.radialLB.rc.temperature = uicontrol(TT{:},h.bc.pl.radialLB.rc(3),...
        'visible','off');
    % --------------------------------------------------------------------
    h.bc.st.axialSG.d1.temperature = uicontrol(TT{:},...
        h.bc.pl.axialSG.d1(3),'pos',[400+25*w 17 135-25*w 17]);
    h.bc.st.axialSG.d2.temperature = uicontrol(TT{:},...
        h.bc.pl.axialSG.d2(3),'pos',[400+25*w 17 135-25*w 17]);
    h.bc.st.axialSG.d1.velocity = uicontrol(TT{:},h.bc.pl.axialSG.d1(1),...
        'pos',[396+25*w 14 135-25*w 17],'string','Velocity Profile: ','visible','off');
    h.bc.st.axialSG.d2.velocity = uicontrol(TT{:},h.bc.pl.axialSG.d2(1),...
        'pos',[396+25*w 14 135-25*w 17],'string','Velocity Profile: ');
    % --------------------------------------------------------------------
    h.bc.st.radialSG.d1.temperature = uicontrol(TT{:},h.bc.pl.radialSG.d1(3));
    h.bc.st.radialSG.d2.temperature = uicontrol(TT{:},h.bc.pl.radialSG.d2(3));
    h.bc.st.radialSG.ro.temperature = uicontrol(TT{:},h.bc.pl.radialSG.ro(3));
    
    
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
        [10 40 450 20],'string','Fluids at rest. Constant mean temperature.',...
        'value',1,'parent',h.sn.bg.initialization);
    h.sn.rb.previous = uicontrol('style','rad','unit','pix','pos',...
        [10 7 190 20],'string','Use previous simulation.','parent',h.sn.bg.initialization);
    
% edit texts
    h.sn.et.residuals(1) = uicontrol(E{:},h.sn.pl.convergence,'string',num2str(h.flowopt.tolerance.residuals,12));
    h.sn.et.residuals(2) = uicontrol(E{:},h.sn.pl.convergence,...
        'string',num2str(h.flowopt.tolerance.newton,12),'pos',[530 21 75 23]);
    h.sn.et.steps = uicontrol(E{:},h.sn.pl.simulation(1),...
        'string',num2str(h.steps),'pos',[110 21 75 23]);
    
% pushbuttons
    B = {'parent',h.sn.pl.simulation(2),'units','pix','pos',[25 24+1*w 76 25],'visible','on','string'};
    h.sn.pb.run = uicontrol(B{:},'Run','parent',h.sn.pl.simulation(1),...
        'pos',[224 20 76 25],'ForegroundColor',h.selectedTabColor);
    h.sn.pb.stop = uicontrol(B{:},'Stop','parent',h.sn.pl.simulation(1),...
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
    E = {'style','edit','units','pix','pos',[60 80 60 23],'parent'};
    % static texts
    S = {'units','pix','style','text','HorizontalAlignment'};
    
% subpanels, subpanels
    if ispc
        txt_wavenumber = 'Wave Numbers along φ';
    else
        txt_wavenumber = 'Wave Numbers in φ';
    end
    h.sy.pl.eigs = uipanel(P{:},[20 585-80*w 325 130],'title','Eigenvalue Solver');
    h.sy.pl.m = uipanel(P{:},[355 585-80*w 145 130],'title',txt_wavenumber);
    h.sy.pl.conv = uipanel(P{:},[510 585-80*w 190 65],'title','Zero Growth Rate');
    h.sy.pl.run(1) = uipanel(P{:},[20 500-80*w 680 75],'title','Linear Stability Analysis');
    h.sy.pl.run(2) = uipanel('units','pix','parent', ...
        h.sy.pl.run(1),'pos',[340 -5 355 85]);
    
% button groups
    h.sy.bg.parameter = uibuttongroup(P{:},[510 653-80*w 190 62],'title','Parameter of Variation');
    
% radio buttons
    h.sy.rb.delta_T = uicontrol('style','rad','unit','pix','pos',...
        [12 15 15 20],'string',' ','value',1,'parent',h.sy.bg.parameter);
    h.sy.rb.T_d1 = uicontrol('style','rad','unit','pix','pos',...
        [76 15 15 20],'string','  ','parent',h.sy.bg.parameter);
    h.sy.rb.T_d2 = uicontrol('style','rad','unit','pix','pos',...
        [138 15 15 20],'string','   ','parent',h.sy.bg.parameter);
    
% edit texts
    h.sy.et.n_eig = uicontrol(E{:},h.sy.pl.eigs,'string',num2str(h.flowopt.eigs.n));
    h.sy.et.n_eig_c = uicontrol(E{:},h.sy.pl.eigs,...
        'string',num2str(h.flowopt.eigs.n_cayley),'pos',[60 45 60 23]);
    h.sy.et.kryl = uicontrol(E{:},h.sy.pl.eigs,...
        'string',num2str(h.flowopt.eigs.krylov),'pos',[250 10 60 23]);
    h.sy.et.tol_eigs = uicontrol(E{:},h.sy.pl.eigs,...
        'string',num2str(h.flowopt.eigs.tol),'pos',[250 80 60 23]);
    h.sy.et.maxit = uicontrol(E{:},h.sy.pl.eigs,...
        'string',num2str(h.flowopt.eigs.maxit),'pos',[250 45 60 23]);
    h.sy.et.m_start = uicontrol(E{:},h.sy.pl.m,...
        'string',num2str(h.m_start),'pos',[70 80 60 23]);
    h.sy.et.m_delta = uicontrol(E{:},h.sy.pl.m,...
        'string','1','pos',[70 45 60 23],'enable','off');
    h.sy.et.m_end = uicontrol(E{:},h.sy.pl.m,...
        'string',num2str(h.m_end),'pos',[70 10 60 23]);
    h.sy.et.conv = uicontrol(E{:},h.sy.pl.conv,...
        'string',num2str(h.flowopt.tolerance.growth),'pos',[89 15 60 23]);
    
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
    h.sy.as.n = axes(AS{:},h.sy.pl.eigs,'pos',h.sy.et.n_eig.Position);
    h.sy.as.n_cay = axes(AS{:},h.sy.pl.eigs,'pos',h.sy.et.n_eig_c.Position);
    h.sy.as.kryl = axes(AS{:},h.sy.pl.eigs,'pos',h.sy.et.kryl.Position);
    h.sy.as.tol_eigs = axes(AS{:},h.sy.pl.eigs,'pos',h.sy.et.tol_eigs.Position);
    h.sy.as.maxit = axes(AS{:},h.sy.pl.eigs,'pos',h.sy.et.maxit.Position);
    h.sy.as.m_start = axes(AS{:},h.sy.pl.m,'pos',h.sy.et.m_start.Position);
    h.sy.as.m_delta = axes(AS{:},h.sy.pl.m,'pos',h.sy.et.m_delta.Position);
    h.sy.as.m_end = axes(AS{:},h.sy.pl.m,'pos',h.sy.et.m_end.Position);
    h.sy.as.delta_T = axes(AS{:},h.sy.bg.parameter,'pos',[29 23 85 20]);
    h.sy.as.T_d1 = axes(AS{:},h.sy.bg.parameter,'pos',[93 23 85 20]);
    h.sy.as.T_d2 = axes(AS{:},h.sy.bg.parameter,'pos',[155 23 85 20]);
    h.sy.as.conv = axes(AS{:},h.sy.pl.conv,'pos',h.sy.et.conv.Position);
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
        'FontSize',14,'FontName',font_name,'string','History','FontWeight','Bold','ForegroundColor',h.selectedTabColor);
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
        'FontSize',14,'FontName',font_name,'string','Licensing','FontWeight','Bold','ForegroundColor',h.selectedTabColor);
    uicontrol(S{:},'left','pos',[30 212+23*w 660 96-24*w],'parent',h.pl.about',...
        'FontSize',12,'FontName',font_name,'string',...
        'This program is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. It comes as is and WITHOUT ANY WARRANTY. By using the program you agree to the above license conditions. For the full text of the license, please visit: https://creativecommons.org/licenses/by-nc-sa/4.0.');
    
    uicontrol(S{:},'left','pos',[15 153+30*w 680 24-2*w],'parent',h.pl.about',...
        'FontSize',14,'FontName',font_name,'string','Contact','FontWeight','Bold','ForegroundColor',h.selectedTabColor);
    uicontrol(S{:},'left','pos',[30 28+60*w 660 120-30*w],'parent',h.pl.about',...
        'FontSize',12,'FontName',font_name,'string',...
        {'Address: Getreidemarkt 9-BA, 1060 Vienna, Austria.','Homepage: http://www.fluid.tuwien.ac.at','Download pages: https://github.com/fromano88/MaranStable, http://www.fluid.tuwien.ac.at/SAJE','E-mail: maranstable@tuwien.ac.at','Help us to improve MaranStable by reporting bugs and/or by giving us feedback.'});

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all edit texts summarized    
h.editTextBoxes = {h.gl.et.g h.gl.et.Vr h.gy.et.ld1 h.gy.et.llb ...
    h.gy.et.ld2 h.gy.et.rc h.gy.et.ri h.gy.et.ro ...
    h.mh.et.radialLB.spacing(1) h.mh.et.radialLB.spacing(2) ...
    h.mh.et.radialLB.spacing(3) h.mh.et.axialLB.spacing(1) ...
    h.mh.et.axialLB.spacing(2) h.mh.et.axialLB.spacing(3) ...
    h.mh.et.radialSG.spacing(1) h.mh.et.radialSG.spacing(2) ...
    h.mh.et.radialSG.spacing(3) h.mh.et.axialSG.spacing(1) ...
    h.mh.et.axialSG.spacing(2) h.mh.et.axialSG.spacing(3) ...
    h.mh.et.axialSG.spacing(4) h.mh.et.axialSG.spacing(5) ...
    h.mh.et.axialSG.spacing(6) h.mh.et.axialSG.spacing(7) ...
    h.mh.et.radialLB.stretchingC(1) h.mh.et.radialLB.stretchingC(3) ...
    h.mh.et.axialLB.stretchingC(1) h.mh.et.axialLB.stretchingC(3) ...
    h.mh.et.radialSG.stretchingC(1) h.mh.et.radialSG.stretchingC(3) ...
    h.mh.et.axialSG.stretchingC(1) h.mh.et.axialSG.stretchingC(3) ...
    h.mh.et.axialSG.stretchingC(5) h.mh.et.axialSG.stretchingC(7) ...
    h.bc.et.axialLB.d1.temperature h.bc.et.axialLB.d2.temperature ...
    h.bc.et.radialLB.rc.temperature h.bc.et.axialSG.d1.temperature ...
    h.bc.et.axialSG.d2.temperature h.bc.et.radialSG.d1.temperature ...
    h.bc.et.radialSG.d2.temperature h.bc.et.radialSG.ro.temperature ...
    h.bc.et.axialSG.d1.velocity h.bc.et.axialSG.d2.velocity ...
    h.sn.et.residuals(1) h.sn.et.residuals(2) h.sn.et.steps ...
    h.or.et.ray h.or.et.particle h.or.et.N(1) h.or.et.N(2) h.or.et.N(3) ...
    h.sy.et.n_eig h.sy.et.n_eig_c h.sy.et.kryl h.sy.et.tol_eigs h.sy.et.maxit h.sy.et.m_start h.sy.et.m_delta h.sy.et.m_end h.sy.et.conv};
