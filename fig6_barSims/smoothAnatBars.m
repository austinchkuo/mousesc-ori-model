% smoothAnatBars - (use with bar simulations) interpolates RF responses from visual coordinates to given anatomical coordinates
%                 also smooths responses if they are provided in grid format
% 
% Usage - [allRandRespsAnat,allRandRespsAnatSmoothed] = smoothAnatBars(respMat,rfParams,anatSampling,smoothSigma,coordSpace,customFile);
%
% Input - respMat: nBarOris x nXLocs x nYLocs matrix of responses from spatially distributed RFs to different bar stimuli
%         rfParams: rfParams struct created by generating rf responses from rf simulator
%         anatSampling: scalar (>0) value specifying desired sampling density along both x and y dims of anatomical grid locations
%                       (required if coordSpace = "grid")
%         smoothSigma: standard deviation of the gaussian smoothing kernel to locally average responses, in pixels
%                      (irrelevant if coordSpace = "nonuniform")
%
% Output - allRespsAnat: nStimOris x nXLocs x nYLocs matrix containing responses to each stimulus orientation at each location
%          allRespsAnatSmoothed: see above but spatially smoothed
%
% Austin Kuo - last update: 9/24/2025

function [allRespsAnat,allRespsAnatSmoothed] = smoothAnatBars(respMat,rfParams,anatSampling,smoothSigma)

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
    error('(smoothAnatBars) Differences in x and y sampling sizes has not been implemented yet.')
end
interParams.stepSizeVis = min(diff(rfParams.xCenter));
interParams.rfRange = [min(rfParams.xCenter) max(rfParams.xCenter)];
if ~isequal([min(rfParams.xCenter) max(rfParams.xCenter)],[min(rfParams.xCenter) max(rfParams.xCenter)])
    error('(smoothAnatBars) Having differences in x and y RF ranges has not been implemented yet.')
end

% corresponding coordinates in visual and anatomical spaces
interParams.xAnchorsAnat = [0.122, 1.126, 1.242, 0.655, 0.789, 1.387, 1.178, 0.655, 1.515]; % anat coords x - anchor points in mm units, bounds are 0-1.658, cartesian coords, origin at bottom left
interParams.yAnchorsAnat = [0.520, 1.156, 0.710, 0.754, 0.377, 1.070, 1.467, 0.091, 0.599]; % anat coords y - anchor points in mm units, bounds are 0-1.653, cartesian coords, origin at bottom left
interParams.xAnchorsVis = [-3, -1, 1, -1,  1,  0, -2,  2, 2]*40/3; % visual space coords x - anchor points; FOV is 80x80 deg, center is (0,0)
interParams.yAnchorsVis = [-3,  1, 1, -1, -1,  2,  2, -2, 2]*40/3; % visual space coords y - anchor points

% set anatomical sample locations
interParams.xAnatGridLims = [0.1,1.75]; % this range has been chosen to center in anatomical space the 80x80deg visual FOV that Liang et al. use
interParams.yAnatGridLims = [0 1.65]; % this range has been chosen to center in anatomical space the 80x80deg visual FOV that Liang et al. use
xAnats = linspace(min(interParams.xAnatGridLims),max(interParams.xAnatGridLims),anatSampling);
yAnats = linspace(max(interParams.yAnatGridLims),min(interParams.yAnatGridLims),anatSampling);
[xAnatGrid, yAnatGrid] = meshgrid(xAnats, yAnats);
interParams.anatCoordsX = xAnatGrid(:);
interParams.anatCoordsY = yAnatGrid(:);

for anatOriIdx = 1:size(respMat,1) 
    % convert responses from visual to anatomical coords
    t = rfAnatVisInterp(interParams,squeeze(respMat(anatOriIdx,:,:)));
    allRespsAnat(anatOriIdx,:,:) = reshape(t,sqrt(length(t)),sqrt(length(t)));
    allRespsAnatSmoothed(anatOriIdx,:,:) = imgaussfilt(squeeze(allRespsAnat(anatOriIdx,:,:)),smoothSigma);
end