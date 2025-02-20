% inject noise pulled from a normal distribution before calculating OSI
% 
% there are miniscule artifacts biasing towards vertical and horizontal oris simply due to the nature of of square pixels
% if you rescale without setting the minimum to 0, you will magnify these effects thousands of times in magnitude
% 0 is the minimum possible response for squared and summed responses

function [noisyResps] = injectNoiseSimGUI(rawResps,deltaFMax,noiseStd)

tempResps = [rawResps(:);0];
tempResps = rescale(tempResps,0,deltaFMax);
tempResps(end) = [];
rescaledResps = reshape(tempResps,size(rawResps));
noise = normrnd(1,noiseStd,size(rescaledResps));
noisyResps = rescaledResps+noise;