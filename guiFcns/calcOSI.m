% calcOSI - calculate orientation preferences, osi, gOSI, and maximum
%           response from a set of responses to oriented grating stimuli
% 
% Usage - [oriPref,osi,gOSI,maxResp] = calcOSI(stimThetas,respMat,normExtent,customNorm);
%
% Input - stimThetas: vector of stimulus orientations that were used, corresponding to the 1st dim of respMat (vec length and dim1 of respMat should match)
%         respMat: nOris x nXLocs x nYLocs matrix of responses from spatially distributed RFs to different orientations
%         normExtent: "single" (default): normalizes gOSI from 0-1 based on each RF at each location's responses
%                     "all": normalizes gOSI from 0-1 across responses from RFs across all spatial locations and orientations
%         customNorm: set to any value user desires to normalize gOSIs by (not required)
%
% Output - oriPref: nXLocs x nYLocs matrix containing orientation preference indices at each location
%          osi: nXLocs x nYLocs matrix containing OSIs at each location
%          gOSI: nXLocs x nYLocs matrix containing gOSIs at each location
%          maxResp: nXLocs x nYLocs matrix containing the value of the maximum response at each location
%                   (i.e., the response to the orientation specified by oriPref)
%
% Austin Kuo - last update: 11/15/24

function [oriPref,osi,gOSI,maxResp] = calcOSI(stimThetas,respMat,normExtent,customNorm)

if ~exist("normExtent","var"); normExtent = "single"; end

[maxResp,maxIdx] = max(respMat);
maxResp = squeeze(maxResp);
maxIdx = squeeze(maxIdx);
[minResp,~] = min(respMat);
minResp = squeeze(minResp);
oriPref = stimThetas(maxIdx);
osi = (maxResp-minResp)./(maxResp+minResp);

tempMat = reshape(respMat,[size(respMat,1),size(respMat,2)*size(respMat,3)]);
numer_gOSI_flat = abs(tempMat'*exp(1i*2*deg2rad(stimThetas(:))));
numer_gOSI = reshape(numer_gOSI_flat,[size(respMat,2),size(respMat,3)]);
if numer_gOSI < 1e-10; numer_gOSI = 0; end % eliminate numerical error because abs(exp(i*2*theta) + exp(i*2*(theta+pi/2))) should be EXACTLY 0, not some infinitesimally number that gets inflated the more terms there are
if strcmp(normExtent,"all") % normalize each RF's gOSI across all RF spatial locations (typically results in very low gOSIs overall)
    denom_gOSI = sum(abs(respMat),"all");
elseif strcmp(normExtent,"single") % normalize each RF's gOSI within its own responses across stimulus orientations
    denom_gOSI = squeeze(sum(abs(respMat)));
end
if exist("customNorm","var") % user sets their own gOSI normalization factor
    denom_gOSI = customNorm;
end
if iscolumn(numer_gOSI)
    numer_gOSI = numer_gOSI';
end
gOSI = numer_gOSI./denom_gOSI;
gOSI(isnan(gOSI)) = 0;

maxResp(isnan(maxResp)) = 0;
