function CMap = myColormap(mapName,nColors)
% possible names: BlueWhiteRed, Cold2Warm, Cool2Warm, Rainbow, Thermal, Thermal2, Thermal3,

if nargin == 1
    nColors = 256;
end

switch mapName
    case 'Cold2Warm'
        CTable = [0.142 0.000 0.850;
                  0.097 0.112 0.970;
                  0.160 0.342 1.000;
                  0.240 0.531 1.000;
                  0.340 0.692 1.000;
                  0.460 0.829 1.000;
                  0.600 0.920 1.000;
                  0.740 0.978 1.000;
                  0.920 1.000 1.000;
                  1.000 1.000 0.920;
                  1.000 0.948 0.740;
                  1.000 0.840 0.600;
                  1.000 0.676 0.460;
                  1.000 0.472 0.340;
                  1.000 0.240 0.240;
                  0.970 0.155 0.210;
                  0.850 0.085 0.187;
                  0.650 0.000 0.130];
              
    case 'Cool2Warm' 
          CTable = diverging_map(linspace(0,1,257),[0.230 0.299 0.754],[0.706 0.016 0.150]);
              
    case 'Thermal'
        CTable = [0.0  0.0  0.0;
                  0.5  0.0  0.0;
                  1.0  0.2  0.0;
                  1.0  1.0  0.0;
                  1.0  1.0  1.0];
              
    case 'Thermal2'
        CTable = [0.0  0.0  0.0;
                  0.3  0.0  0.7;
                  1.0  0.2  0.0;
                  1.0  1.0  0.0;
                  1.0  1.0  1.0];
              
    case 'Thermal3'
        CTable = [1.0  1.0  0.0;
                  1.0  0.2  0.0;
                  0.3  0.0  0.7;
                  0.0  0.0  0.0];

    case 'BlueWhiteRed'
        CTable = [0.2300 0.2990 0.7540;
                  0.5532 0.6890 0.9954;
                  1.0000 1.0000 1.0000;
                  0.9583 0.6030 0.4819;
                  0.7060 0.0161 0.1500];

    case 'Rainbow'
        CTable = [  5  97 255;
                    5 108 247;
                    5 119 240;
                    5 130 232;
                    5 139 223;
                    5 149 213;
                    5 158 203;
                    5 166 191;
                    5 175 179;
                    5 184 168;
                    5 193 154;
                    5 202 140;
                    5 212 127;
                    5 220 109;
                    6 229  92;
                    4 237  75;
                   70 243  39;
                  126 245  28;
                  164 249  12;
                  194 251   9;
                  225 253   6;
                  255 255   3;
                  255 244  20;
                  255 232  38;
                  255 221  55;
                  255 209  55;
                  255 196  55;
                  255 184  55;
                  255 172  55;
                  255 159  55;
                  255 147  55;
                  255 133  55;
                  255 119  55;
                  255 104  55;
                  254  85  54;
                  252  66  48;
                  253  38  54;
                  242  30  64;
                  230  20  74;
                  218  10  84;
                  204  11  91;
                  189  12  98;
                  174  13 106];
              
        CTable = CTable./255;
              
    otherwise
        error('Chosen colormap not implemented. Choose among:''Cold2Warm'',''BlueWhiteRed'',''Thermal'',''Thermal2''')
end

if strcmp(mapName,'Blue_White_Red')
    lims = get(gca, 'CLim');
    % Find ratio of negative to positive
    if lims(1)*lims(2) < 0
        % Find ratio of negative to positive
        ratio = abs(lims(1))/(abs(lims(1)) + lims(2));
        nCol(1) = round(nColors*ratio);
        nCol(2) = nColors-nCol(1);
        for i = 1:2
            newTable{i} = [CTable(2*i-1,:); CTable(2*i,:); CTable(2*i+1,:)]; %#ok<AGROW>
            [l,~] = size(newTable{i});
            [x0,y0] = meshgrid(1:3,linspace(0,1,l));
            [xn,yn] = meshgrid(1:3,linspace(0,1,nCol(i)));
            newMap{i} = interp2(x0,y0,newTable{i},xn,yn); %#ok<AGROW>
        end
        
        CMap = [newMap{1}; newMap{2}];
        
    else
        error('Colormap ''BlueWhiteRed'' can only be used if the values are partially positive and negative. Use colormap ''Cold2Warm'' instead.')
        
    end
else
    [l,~] = size(CTable);
    [x0,y0] = meshgrid(1:3,linspace(0,1,l));
    [xn,yn] = meshgrid(1:3,linspace(0,1,nColors));

    CMap = interp2(x0,y0,CTable,xn,yn);
end