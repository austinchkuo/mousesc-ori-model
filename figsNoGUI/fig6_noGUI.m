% fig6_noGUI.m

% (center-surround): use fig4_noGUI
% RF parameters:
% C/S size ratio: 0.8 to 4.7 deg (resulting center frequency: 0.08 cpd)
% surround amplitude: 0.4

% (V1-like):
% RF parameters:
% LW ratio: 4 to 2
% # subregions: 2
% center frequency: 0.08 cpd

clear
close all

% simulated response parameters
deltaFMax = 1;
noiseStd = 0.01;
smoothingSigmaVis = 0.44; % std of visual smoothing kernel, in deg
smoothingSigmaAnat = 0.0083; % std of anatomical smoothing kernel, in mm
anatSamples = 241; % density of the anatomical grid (anatSamples x anatSamples); NOTE: main variable determining function runtime
sigmaAnatPix = anatSamples*smoothingSigmaAnat / 1.66; % 1.66mm is the FOV of the anatomy (see Liang et al., 2023; Fig. 1c)

% visualization parameters
gOSIThresh = 0.1;
alphaExp = 0.00001;
disableLabels = 1;

% specify data file
dirstring = "~/Documents/MATLAB/data/rfSimData/";
dataFile = "rfGaborSimData_25-01-22_1546.mat"; % 0.04 cpd, 30 deg radius circle
loadFilename = strcat(dirstring,dataFile);
fprintf("Loading data from: %s\n",loadFilename)
load(loadFilename,"d","p")
visSamples = length(p.rfParams.xCenter);
visFOV = abs(p.rfParams.xCenter(end) - p.rfParams.xCenter(1));
rfGridScaleFactor = visSamples/visFOV; % number of steps/deg
sigmaVisPix = rfGridScaleFactor*smoothingSigmaVis;

% get RF responses averaged across stimulus phases
fprintf('Averaging responses across %0.00f stimulus phases... ',length(p.stimParams.phase))
rfRespPhaseAvg = avgResponsesAcrossStimPhase(p.rfParams,p.stimParams,d.responses);
fprintf('...done.\n\n')

% normalize and inject noise
fprintf('Rescaling stimulus-driven responses between 0-1, adding baseline response, and adding gaussian noise...\n')
rfRespNormIndVis = injectNoiseSimGUI(rfRespPhaseAvg,deltaFMax,noiseStd);
fprintf('done.\n\n')

% fix rotated response matrix (see printed statements)
fprintf('The code that generates responses rotates the coordinates by exactly 90 degrees...\n')
fprintf('Post-hoc undo for now since initial rotation doesn''t affect overall responses, but low priority fix later...\n')
fprintf('Rotating responses spatially by -90 degrees...\n')
for i = 1:size(rfRespNormIndVis,4)
    rfRespNormIndVisTemp(:,:,:,i) = fixRotation(rfRespNormIndVis(:,:,:,i));
end
rfRespNormIndVis = rfRespNormIndVisTemp;
clear rfRespNormIndVisTemp
fprintf('...done.\n\n')

% choose random individual units (vis)
fprintf('Choosing random RF orientations (visual)...\n')
rng(10) % but choose the same ones to make comparisons between conditions
randIdx = randi(size(rfRespNormIndVis,4),size(rfRespNormIndVis,2),size(rfRespNormIndVis,3));
for stimOriIdx = 1:size(rfRespNormIndVis,1)
    for rfPosIdx_x = 1:size(rfRespNormIndVis,2)
        for rfPosIdx_y = 1:size(rfRespNormIndVis,3)
            rfRespNormIndVisRand(stimOriIdx,rfPosIdx_x,rfPosIdx_y) = rfRespNormIndVis(stimOriIdx,rfPosIdx_x,rfPosIdx_y,randIdx(rfPosIdx_x,rfPosIdx_y));
        end
    end
end
fprintf('...done.\n\n')

fprintf('Smoothing single unit RF responses in visual coordinates...\n')
rfRespNormPopVisRand = spatialSmoothResponses(p,rfRespNormIndVisRand,sigmaVisPix);
fprintf('...done.\n\n')

% calculate orientation preferences, gOSIs, max responses in visual coordinates
fprintf('Calculating orientation preferences, gOSIs, and max responses for RF responses in visual coordinates...\n')
[oriPrefIndVis,~,gOSIIndVis,maxRespIndVis] = calcOSI(p.stimParams.theta,rfRespNormIndVisRand);
[oriPrefPopVis,~,gOSIPopVis,maxRespPopVis] = calcOSI(p.stimParams.theta,rfRespNormPopVisRand);
gOSINormIndVis = rescaleOSIs(gOSIIndVis);
gOSINormPopVis = rescaleOSIs(gOSIPopVis);
fprintf('...done.\n\n')

% generate anatomical coordinates and estimate responses in anat coords
tic
fprintf('Generating anatomical coordinates and estimating responses...\n')
for rfOri = 1:size(rfRespNormIndVis,4)
    [rfRespNormIndAnat(:,:,:,rfOri),~,mmScale] = smoothAnatGUI(rfRespNormIndVis(:,:,:,rfOri),p.rfParams,anatSamples,sigmaAnatPix);
end
toc
fprintf('...done.\n\n')

% choose random individual units (anat)
fprintf('Choosing random RF orientations (anat)...\n')
randIdx = randi(size(rfRespNormIndAnat,4),size(rfRespNormIndAnat,2),size(rfRespNormIndAnat,3));
for stimOriIdx = 1:size(rfRespNormIndAnat,1)
    for rfPosIdx_x = 1:size(rfRespNormIndAnat,2)
        for rfPosIdx_y = 1:size(rfRespNormIndAnat,3)
            rfRespNormIndAnatRand(stimOriIdx,rfPosIdx_x,rfPosIdx_y) = rfRespNormIndAnat(stimOriIdx,rfPosIdx_x,rfPosIdx_y,randIdx(rfPosIdx_x,rfPosIdx_y));
        end
    end
end
fprintf('...done.\n\n')

fprintf('Smoothing single unit RF responses in anatomical coordinates...\n')
rfRespNormPopAnatRand = spatialSmoothResponses(p,rfRespNormIndAnatRand,sigmaAnatPix);
% for stimOri = 1:size(rfRespNormIndAnatRand,1)
%     rfRespNormPopAnatRand(stimOri,:,:) = imgaussfilt(squeeze(rfRespNormIndAnatRand(stimOri,:,:)),sigmaAnatPix);
% end
fprintf('...done.\n\n')

% calculate orientation preferences, gOSIs, max responses in anatomical coordinates
fprintf('Calculating orientation preferences, gOSIs, and max responses for RF responses in anatomical coordinates...\n')
[oriPrefIndAnat,~,gOSIIndAnat,maxRespIndAnat] = calcOSI(p.stimParams.theta,rfRespNormIndAnatRand);
[oriPrefPopAnat,~,gOSIPopAnat,maxRespPopAnat] = calcOSI(p.stimParams.theta,rfRespNormPopAnatRand);
gOSINormIndAnat = rescaleOSIs(gOSIIndAnat);
gOSINormPopAnat = rescaleOSIs(gOSIPopAnat);
fprintf('...done.\n\n')

plotOriPrefAnat(oriPrefPopAnat,gOSINormPopAnat,anatSamples,gOSIThresh,alphaExp,mmScale,disableLabels)

%% %%%%%%%%%%%%%%%%%% %%
%%% helper functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%
function cr_phase_avg = avgResponsesAcrossStimPhase(rfParams,stimParams,responses)
    % complex responses averaged across stimulus phases
    reshapeSize = nonzeros([length(stimParams.theta),length(rfParams.xCenter),length(rfParams.yCenter),length(rfParams.theta)])';
    cr_phase_avg = reshape(mean(responses.complex_response,2),reshapeSize);

end

function smoothedResps = spatialSmoothResponses(p,unsmoothedResps,smoothingSigma)
    % unsmoothedResps should be size (nStimOris, nXLocations, nYlocations)
    
    % ensure non-singleton dimension is found
    xyDims = [length(p.rfParams.xCenter),length(p.rfParams.yCenter)];
    if isequal(xyDims,[1,1])
        error("You can't smooth over a single RF")
    else
        [nSteps,XorY] = max(xyDims);
        if XorY == 1
            rfSpan = abs(p.rfParams.xCenter(end) - p.rfParams.xCenter(1));
        elseif XorY == 2
            rfSpan = abs(p.rfParams.yCenter(end) - p.rfParams.yCenter(1));
        end
        rfGridScaleFactor = nSteps/rfSpan; % number of steps/deg
    end
    
    for stimOriIdx = 1:size(unsmoothedResps,1)
        % smooth RF responses for every stimulus orientation
        smoothedResps(stimOriIdx,:,:,:) = imgaussfilt(squeeze(unsmoothedResps(stimOriIdx,:,:,:)),smoothingSigma,'Padding',0);
    end

end

function rotResp = fixRotation(preResp)
    tempResp = permute(preResp,[2,3,1]);
    tempResp = rot90(tempResp,-1);
    rotResp = permute(tempResp,[3,1,2]);

end

function scaledOSIs = rescaleOSIs(preOSIs)
    temp = rescale([preOSIs(:);0]);
    temp(end) = [];
    scaledOSIs = reshape(temp,size(preOSIs,1),size(preOSIs,2));

end

function [] = plotOriPrefAnat(oriPrefs,OSIs,anatSamples,osiThreshold,alphaExponent,mmScale,disableLabels)

    xAnats = linspace(0.1,1.75,anatSamples); % these endpoints are chosen so that the image is centered
    yAnats = linspace(1.65,0,anatSamples);   % 1.66mm is the FOV of the anatomy (see Liang et al., 2023; Fig. 1c)
    [X, Y] = meshgrid(xAnats, yAnats);
    xLocs = X(:);
    yLocs = Y(:);
    sColor = oriPrefs(:);
    tempTable = table(xLocs,yLocs,sColor);
    s_oriPref = scatter(tempTable,'xLocs','yLocs','Filled','ColorVariable','sColor');
    hold("on")
    % quick hack for scalebar since grid and custom points are handled differently
    % line([mmScale(1) mmScale(2)],[mmScale(3) mmScale(3)],'LineWidth',2,'Color','black');
    colormap('hsv')
    clim([0 180])
    thisAlpha = OSIs(:);
    thisAlpha(thisAlpha<osiThreshold) = 0;
    thisAlpha = rescale(thisAlpha);
    s_oriPref.AlphaData = thisAlpha.^alphaExponent;
    s_oriPref.MarkerFaceAlpha = 'flat';
    s_oriPref.AlphaDataMapping = 'none';
    s_oriPref.SizeData = 10;
    set(gca,'xlim',[min(xLocs) max(xLocs)])
    set(gca,'ylim',[max(yLocs)-range(xLocs) max(yLocs)])
    axis('square')
    
    % use this to disable axis lines and labels
    if disableLabels
        set(gca,'xColor','w')
        set(gca,'yColor','w')
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        set(get(gca,'XAxis'), 'Visible', 'off')
        set(get(gca,'YAxis'), 'Visible', 'off')
        xlabel('')
        ylabel('')
    end
    
    hold('off')

end