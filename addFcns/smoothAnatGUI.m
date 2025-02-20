% smoothAnatGUI - (use with rf simulation GUI) interpolates RF responses from visual coordinates to given anatomical coordinates
%                 also smooths responses if they are provided in grid format
% 
% Usage - [allRandRespsAnat,allRandRespsAnatSmoothed,mmScaleIdxs] = smoothAnatGUI(respMat,rfParams,anatSampling,smoothSigma,coordSpace,customFile);
%
% Input - respMat: nOris x nXLocs x nYLocs matrix of responses from spatially distributed RFs to different orientations
%         rfParams: rfParams struct created by generating rf responses from rf simulator
%         anatSampling: scalar (>0) value specifying desired sampling density along both x and y dims of anatomical grid locations
%                       (required if coordSpace = "grid")
%         smoothSigma: standard deviation of the gaussian smoothing kernel to locally average responses, in pixels
%                      (irrelevant if coordSpace = "nonuniform")
%         coordSpace: either "grid" or "nonuniform"
%                     "grid" creates a square of sampled locations, with sampling density specified by anatSampling
%                     "nonuniform" assumes user-inputted anatomical rf locations (see customFile field)
%         customFile: either a table with header "x" for desired x locations and header "y" for desired y locations
%                     or a string specifying path to a custom .csv file that can be loaded as a table as specified above
%         customStarts: a 2-vector containing x and y offset values to add to custom X and Y coordinates
%
% Output - allRespsAnat: nStimOris x nXLocs x nYLocs matrix containing responses to each stimulus orientation at each location
%          allRespsAnatSmoothed: see above but spatially smoothed
%          scaleIdxs: gives left and right indices of scale bar to plot in the anatomical grid
%
% Austin Kuo - last update: 11/15/24

function [allRespsAnat,allRespsAnatSmoothed,scaleCoords] = smoothAnatGUI(respMat,rfParams,anatSampling,smoothSigma,coordSpace,customFile,customStarts)

p = gcp('nocreate');
if isempty(p)
    parpool('threads');
end

if ~exist("anatSpacing","var")
    interParams.anatSampling = 81; % default
else
    interParams.anatSampling = anatSampling;
end
interParams.stepSizeVis = min(diff(rfParams.xCenter));
if min(diff(rfParams.xCenter)) ~= min(diff(rfParams.yCenter))
    error('(smoothAnatGUI) Differences in x and y sampling sizes has not been implemented yet.')
end
interParams.stepSizeVis = min(diff(rfParams.xCenter));
interParams.rfRange = [min(rfParams.xCenter) max(rfParams.xCenter)];
if ~isequal([min(rfParams.xCenter) max(rfParams.xCenter)],[min(rfParams.xCenter) max(rfParams.xCenter)])
    error('(smoothAnatGUI) Having differences in x and y RF ranges has not been implemented yet.')
end
if ~exist("coordSpace","var")
    coordSpace = "grid"; % default
elseif strcmpi(coordSpace,"nonuniform")

end

% corresponding coordinates in visual and anatomical spaces
% x_anat_anchor = [-214.7799, 85.8491, 126.1006, -53.7735, -16.0377, 188.9937, 98.4277, -51.2579, 217.9245]; % anat coords x - anchor points
% y_anat_anchor = [-97.5083, 99.0161, -34.983, -18.2177, -137.0511, 82.0783, 198.9129, -226.7704, -78.237]; % anat coords y - anchor points
interParams.xAnchorsAnat = [0.122, 1.126, 1.242, 0.655, 0.789, 1.387, 1.178, 0.655, 1.515]; % anat coords x - anchor points in mm units, bounds are 0-1.658, cartesian coords, origin at bottom left
interParams.yAnchorsAnat = [0.520, 1.156, 0.710, 0.754, 0.377, 1.070, 1.467, 0.091, 0.599]; % anat coords y - anchor points in mm units, bounds are 0-1.653, cartesian coords, origin at bottom left
interParams.xAnchorsVis = [-3, -1, 1, -1,  1,  0, -2,  2, 2]*40/3; % visual space coords x - anchor points; FOV is 80x80 deg, center is (0,0)
interParams.yAnchorsVis = [-3,  1, 1, -1, -1,  2,  2, -2, 2]*40/3; % visual space coords y - anchor points

% set anatomical sample locations
interParams.xAnatGridLims = [0.1,1.75]; % this range has been chosen to center in anatomical space the 80x80deg visual FOV that Liang et al. use
interParams.yAnatGridLims = [0 1.65]; % this range has been chosen to center in anatomical space the 80x80deg visual FOV that Liang et al. use

switch lower(coordSpace)
    case "grid"
        xAnats = linspace(min(interParams.xAnatGridLims),max(interParams.xAnatGridLims),anatSampling);
        yAnats = linspace(max(interParams.yAnatGridLims),min(interParams.yAnatGridLims),anatSampling);
        [xAnatGrid, yAnatGrid] = meshgrid(xAnats, yAnats);
        interParams.anatCoordsX = xAnatGrid(:);
        interParams.anatCoordsY = yAnatGrid(:);
        
        % scale bar
        endVal = 0.95*(max(interParams.xAnatGridLims) - min(interParams.xAnatGridLims)) + min(interParams.xAnatGridLims); % end scale bar 5 percent from the right side
        heightVal = 0.025*(max(interParams.yAnatGridLims) - min(interParams.yAnatGridLims)) + min(interParams.yAnatGridLims); % set height of scale bar to be 5 percent from bottom
        scaleCoords = [endVal-0.25,endVal,heightVal];

    case "nonuniform"
        
        if ~exist("customStarts","var")
            nonuniformStartX = 0;
            nonuniformStartY = 0;
        else
            nonuniformStartX = customStarts(1);
            nonuniformStartY = customStarts(2);
        end

        if ischar(customFile) || isstring(customFile)
            temp = readtable(customFile);
            temp.x = temp.x/3;
            temp.y = temp.y/3;
            interParams.anatCoordsX = temp.x + (nonuniformStartX - min(temp.x));
            interParams.anatCoordsY = temp.y + (nonuniformStartY - min(temp.y));
        elseif istable(customFile)
            interParams.anatCoordsX = customFile.x + (nonuniformStartX - min(customFile.x));
            interParams.anatCoordsY = customFile.y + (nonuniformStartY - min(customFile.y));
        end

        % scale bar
        endVal = 0.95*(max(interParams.anatCoordsX) - min(interParams.anatCoordsX)) + min(interParams.anatCoordsX); % end scale bar 5 percent from the right side
        minY = max(interParams.anatCoordsY) - range(interParams.anatCoordsX);
        totalHeightY = max(interParams.anatCoordsY) - minY;
        heightVal = 0.05*totalHeightY + minY; % set height of scale bar to be 5 percent from bottom
        scaleCoords = [endVal-0.05,endVal,heightVal]; % make scale bar 50 microns
        
end

for anatOriIdx = 1:size(respMat,1) % to-do: figure out a way to unify the fact that you either put in a 2d matrix of coordinates or a vector (which can't be smoothed)
    % convert responses from visual to anatomical coords
    if strcmpi(coordSpace,"grid")
        t = rfAnatVisInterp(interParams,squeeze(respMat(anatOriIdx,:,:)));
        allRespsAnat(anatOriIdx,:,:) = reshape(t,sqrt(length(t)),sqrt(length(t)));
        allRespsAnatSmoothed(anatOriIdx,:,:) = imgaussfilt(squeeze(allRespsAnat(anatOriIdx,:,:)),smoothSigma);
    elseif strcmpi(coordSpace,"nonuniform") % a vector of coordinates is fed in, so we need to save out x and y values
        allRespsAnat(anatOriIdx,:) = rfAnatVisInterp(interParams,squeeze(respMat(anatOriIdx,:,:)));
        allRespsAnatSmoothed(anatOriIdx,:) = nan;
    end

end