function [nCells, aspectRatio, skewness] = gridQuality(blocks)
% gridQuality returns grid properties such as:
    %  - number of cells
    %  - maximal aspect ratio
    %  - maximal skewness

nCells = 0;
aspectRatio = 0;
skewness = 0;

b1 = evalin('base','b1'); %#ok<*NASGU> 
if length(blocks) == 2
    b2 = evalin('base','b2');
end
    
for n=1:length(blocks)
    b = eval(blocks{n});
    % numer of cells
        nCells = nCells + numel(b.R.p);
    
    % axial distances between two adjacent cells
        axialDistance = diff(b.Z.v');
        
    % lengths of the cells: north and south
        lengths = sqrt(diff(b.R.v').^2+diff(b.Z.v').^2);
    
    % heights of the cells: east and west
        heights = diff(b.R.v);
        
    % size of the matrix
        [l, k] = size(lengths);
    
    for j = 1:k-1
        for i = 1:l
            % aspect ratio
                aR = max([lengths(i,j) lengths(i,j+1) heights(j,i) heights(j,i+1)])/min([lengths(i,j) lengths(i,j+1) heights(j,i) heights(j,i+1)]);
                aspectRatio = max(aR,aspectRatio);
            
            % skewness = max((teta_max-90°)/90°,(90°-teta_min)/90°)
            % cells are in general trapezoids because radial lines are always parallel
            % alpha - beta - gamma - delta of each cell not as important as deviation of these angles from 90° -> denoted by '_D'
                alphas_D = acos(axialDistance(i,j)/lengths(i,j));
                betas_D  = pi/2 - asin(axialDistance(i,j)/lengths(i,j));
                gammas_D = pi/2 - asin(axialDistance(i,j)/lengths(i,j+1));
                deltas_D = acos(axialDistance(i,j)/lengths(i,j+1));
                sk       = max([alphas_D betas_D gammas_D deltas_D])/(pi/2);
                skewness = max(sk, skewness);
        end
    end
end

end