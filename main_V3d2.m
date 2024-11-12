function main_V3d2

evalin('base','clearvars'); close all; clc

% comment next line for standalone applications and put all files in the one folder
% evalin('base','addpath(''./src'');')

MS_version = 'MaranStable 3.2';

f_open = figure('units','pixels','menubar','none','name',...
    MS_version,'numbertitle','off','resize','off','Visible','Off');

w = ispc; % w = 1 --> windows
if ~ispc && ~isunix % neither windows nor linux
    waitfor(errordlg('MaranStable is only supported by Windows and Linux!','Error 062'))
    delete(f_open)
    return
end
if ispc
    f_size = [399 399];
else
    f_size = [388 399];
end
f_open.Position(3:4) = f_size;
movegui(f_open,'center')

% background
    c_table  = [[128 184 222]; [255 255 255]; [236 169 140]]./255;
    [x0, y0] = meshgrid(1:3,linspace(0,1,3));
    [xn, yn] = meshgrid(1:3,linspace(0,1));
    c_map    = interp2(x0,y0,c_table,xn,yn);
    axes('unit','pix','visible','off','pos',[0 0 400 400]);
    fill([0 200 400 400 200 0],[0 0 0 400 400 400],[1 0.5 0 0 0.5 1],'LineStyle','none')
    colormap(c_map)
    axis off

% MS logo and welcome text
    axes('unit','pix','visible','off','pos',[8 20 384 360]);
    logo_MS
    text(-3.35,5.4,['Welcome to ' MS_version],'interpreter','Latex','FontSize',17);

% TU Wien logo
    axes('unit','pix','visible','off','pos',[357 5 40 40]);
    imshow(imread('logo_TU_splash.png'));
    
% isw Wien logo
    axes('unit','pix','visible','off','pos',[5 357 45 45]);
    imshow(imread('logo_isw_splash.png'));

h.model = 'TFM';

% button group: model
    h.bg.model = uibuttongroup('units','pix','pos',[102 168 180 92],...
        'BackgroundColor',[1 1 1],'SelectionChangedFcn',{@bg_model,h},'visible','off');
    h.rb.SFM = uicontrol(h.bg.model,'style','rad','unit','pix','pos',...
        [15+12*w 53 150 20],'string','Single-Fluid Model','BackgroundColor',[1 1 1]);
    h.rb.TFM = uicontrol(h.bg.model,'style','rad','unit','pix','pos',...
        [15+12*w 17 150 20],'string','Two-Fluid Model','value',1,'BackgroundColor',[1 1 1]);

    
% pushbutton
    h.pb.ok = uicontrol('units','pix','pos',[194 143 88 25],'string','Go','call',{@pb_ok,h},'visible','off');

% correct figure size
    pause(0.1)
    f_open.Position(3:4) = f_size;

% little pause before the button group and the push button are shown
set(f_open,'Visible','On')
pause(4)
set([h.bg.model h.pb.ok],'Visible','on')

uiwait(f_open)

    function bg_model(varargin)
        selection = get(get(h.bg.model,'SelectedObject'),'String');

        switch selection
            case 'Single-Fluid Model'
                h.model = 'SFM';

            case 'Two-Fluid Model'
                h.model = 'TFM';
        end
    end

    function pb_ok(varargin)
        delete(f_open)

        if strcmp(h.model,'TFM')
            evalin('base','TwoFluidModel_V3d2')
        else
            evalin('base','SingleFluidModel_V3d2')
        end
    end
end
