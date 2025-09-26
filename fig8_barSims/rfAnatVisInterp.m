function [varargout] = rfAnatVisInterp(interParams,varargin)

if length(varargin) == 2
    interpHSVA = 1;
    hsvImageMatrix = varargin{1};
    alphaImageMatrix = varargin{2};
elseif length(varargin) == 1
    interpHSVA = 0;
    valuesImageMatrix = varargin{1};
else
    error('(rfAnatVisInterp) Wrong number of varargins (1 for standard value interp, 2 for HSVA interp).')
end

x_anat_anchor = interParams.xAnchorsAnat;
y_anat_anchor = interParams.yAnchorsAnat;
x_vis_anchor = interParams.xAnchorsVis;
y_vis_anchor = interParams.yAnchorsVis;
visSpacing = interParams.stepSizeVis;
visMin = interParams.rfRange(1);
visMax = interParams.rfRange(2);
anatCoordsX = interParams.anatCoordsX;
anatCoordsY = interParams.anatCoordsY;

%% scatteredInterpolant
% estimate functions that predict visual coordinates given anatomical coordinates
xvisEst = scatteredInterpolant(x_anat_anchor', y_anat_anchor', x_vis_anchor', 'natural'); % est x-coord
yvisEst = scatteredInterpolant(x_anat_anchor', y_anat_anchor', y_vis_anchor', 'natural'); % est y-coord

% estimate x and y coords in visual space given a set of anatomical coordinates
for i = 1:size(anatCoordsX,1)
    xVisEstCoords(i) = xvisEst(anatCoordsX(i),anatCoordsY(i));
    yVisEstCoords(i) = yvisEst(anatCoordsX(i),anatCoordsY(i));
end

xVisEstCoords = xvisEst(anatCoordsX(:),anatCoordsY(:));
yVisEstCoords = yvisEst(anatCoordsX(:),anatCoordsY(:));

if interpHSVA == 1
    
    % not used
    %{
    % find interpolated alpha and color values from visual coordinates' RF properties
    xVisCoords = linspace(visMax,visMin,size(hsvImageMatrix,1));
    yVisCoords = linspace(visMin,visMax,size(hsvImageMatrix,1));
    [yVisGrid, xVisGrid] = ndgrid(xVisCoords, yVisCoords);

    % Temporary storage to hold values computed in parfor
    nRows = size(anatCoordsX,1);
    nCols = size(anatCoordsX,2);
    tempValues = cell(1, nRows*nCols);

    parfor index = 1:nRows*nCols
        iX = ceil(index / nCols); % Calculate row index
        iY = mod(index - 1, nCols) + 1; % Calculate column index

        % Perform the interpolation and store in temporary array
        tempValues{index} = InterpolateColorAlpha(xVisEstCoords(iX,iY), yVisEstCoords(iX,iY), xVisGrid, yVisGrid, visSpacing, hsvImageMatrix, alphaImageMatrix);
    end

    % Reshape temporary array back into the desired 2D cell array
    for index = 1:nRows*nCols
        iX = ceil(index / nCols);
        iY = mod(index - 1, nCols) + 1;
        anatColorAlpha{iX, iY} = tempValues{index};
    end

    % Assign colors and alpha values to the grid
    k = 1; % Index for anatColorAlpha
    for j = 1:length(anatCoordsY) % Iterate over y-dimension
        for i = 1:length(anatCoordsX) % Iterate over x-dimension
            row = i;
            col = j;

            % Assign color and alpha value
            colorImage(row, col, :) = anatColorAlpha{k}(1:3);
            alphaImage(row, col) = anatColorAlpha{k}(4);

            k = k + 1; % Increment the index for anatColorAlpha
        end
    end

    varargout{1} = colorImage;
    varargout{2} = alphaImage;
    %}

else

    % find interpolated values from visual coordinates' RF properties
    xVisCoords = linspace(visMin,visMax,size(valuesImageMatrix,1));
    yVisCoords = linspace(visMax,visMin,size(valuesImageMatrix,1));
    
    [xVisGrid, yVisGrid] = meshgrid(xVisCoords, yVisCoords);
    
    % % find interpolated values from visual coordinates' RF properties
    % xVisCoords = linspace(visMax,visMin,size(valuesImageMatrix,1));
    % yVisCoords = linspace(visMin,visMax,size(valuesImageMatrix,1));
    % 
    % [yVisGrid, xVisGrid] = ndgrid(xVisCoords, yVisCoords);
    % yVisGrid = flipud(yVisGrid);
    % xVisGrid = fliplr(xVisGrid);

    valuesImage = interp2(xVisGrid,yVisGrid,valuesImageMatrix,xVisEstCoords,yVisEstCoords);
    
    varargout{1} = valuesImage;

    if any(isnan(valuesImage))
        % keyboard
    end

end

end
