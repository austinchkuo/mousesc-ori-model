% createQuadratureLogGaussianFiltersSpace.m

function [bpP1,bpfP1,bpP2,bpfP2,xCenterPix,yCenterPix] = createQuadratureLogGaussianFiltersSpace(fovPix,fovDeg,g1_sigma,g2_sigma,xCenter,yCenter,g2_amp)

if mod(fovPix,2) == 0
    fovPix = fovPix + 1;
end

if ~exist('g2_amp','var')
    g2_amp = 1;
end

if g1_sigma > g2_sigma % always assign g2_sigma to be the larger of the two inputs
    temp = g1_sigma;
    g1_sigma = g2_sigma;
    g2_sigma = temp;
end

deltaS = 1/(fovPix-1);
x = (0:deltaS:1)-0.5;
y = fliplr((0:deltaS:1)-0.5);
[X,Y] = meshgrid(x,y);

%%% gaussian envelope %%%
g1_sigma = g1_sigma/fovDeg;
g2_sigma = g2_sigma/fovDeg;
% generate envelopes in the center of the FOV
g1 = (1/(2*pi*g1_sigma^2)) * exp(-(X.^2 / (2 * (g1_sigma^2)) + Y.^2 / (2 * (g1_sigma^2))));
g2 = g2_amp * (1/(2*pi*g2_sigma^2)) * exp(-(X.^2 / (2 * (g2_sigma^2)) + Y.^2 / (2 * (g2_sigma^2))));

% calculate number of indices to shift based on X and Y center
fovDegSpacing = min(diff(linspace(-fovDeg/2,fovDeg/2,fovPix)));
indsPerDeg = 1/fovDegSpacing; % indices/deg
xShift = xCenter*indsPerDeg;
yShift = yCenter*indsPerDeg;

dogPreShift = g1-g2;
% dogPreMask = circshift(dogPreShift, round([-yShift, xShift]));
dogPreMask = shiftImage2D(dogPreShift,round(xShift),round(yShift));
mask = abs(dogPreMask) > 0.001;
dog = dogPreMask.*mask;

xCenterPix = round(xShift) + round(fovPix/2);
yCenterPix = -round(yShift) + round(fovPix/2);

% deg2pix = (fovPix-1)/fovDeg; % conversion from deg to pix (scale factor of pix/deg)
% xCenterPix = xCenter*deg2pix;
% yCenterPix = (-1)*yCenter*deg2pix; % -1 is needed because of matlab's top to bottom indexing for low to high y-values
% sig1pix = g1_sigma * deg2pix;
% sig2pix = g2_sigma * deg2pix;
% 
% g1 = fspecial('Gaussian',2*ceil(3*sig2pix)+1,sig1pix);
% g2 = g2_amp * fspecial('Gaussian',2*ceil(3*sig2pix)+1,sig2pix);
% 
% dog = g1 - g2;
% 
% shiftValueX = ceil(fovPix/2);
% 
% xShift = round(xCenterPix + shiftValueX);
% yShift = round(yCenterPix + shiftValueX);
% halfSize = floor(size(g2,1)/2);
% xStart = xShift - halfSize;
% xEnd = xShift + halfSize;
% yStart = yShift - halfSize;
% yEnd = yShift + halfSize;
% fov = zeros(fovPix,fovPix);
% fov(yStart:yEnd, xStart:xEnd) = dog;

d1 = getHalfFourier(dog);
[bpP1,bpfP1] = reconstructFromHalfFourier(d1);

d2 = getHalfFourier(-dog);
% d2.phase = wrapValue(d2.phase + pi,-pi,pi);
[bpP2,bpfP2] = reconstructFromHalfFourier(d2);
