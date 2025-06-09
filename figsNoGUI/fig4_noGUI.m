% fig4_noGUI.m

% RF parameters:
% C/S size ratio: 0.8 : 4.7 deg (resulting center frequency: 0.08 cpd)
% surround amplitude: 0.4

clear
close all

% simulated response parameters
deltaFMax = 1;
noiseStd = 0.01;
smoothingSigmaVis = 0.44; % std of visual smoothing kernel, in deg
smoothingSigmaAnat = 0.0083; % std of anatomical smoothing kernel, in mm
anatSamples = 161; % density of the anatomical grid (anatSamples x anatSamples); NOTE: main variable determining function runtime
sigmaAnatPix = anatSamples*smoothingSigmaAnat / 1.66; % 1.66mm is the FOV of the anatomy (see Liang et al., 2023; Fig. 1c)

% visualization parameters
gOSIThresh = 0.5;
alphaExp = 1.5;
disableLabels = 1;

% specify data file (different files: change stimulus parameters, same RF parameters)
dirstring = "~/Documents/MATLAB/rfSimData/fig4/";
dataFile = "rfGaborSimData_24-12-11_1327.mat";
% data for figure 4
% dataFile = "rfGaborSimData_24-12-11_1410.mat"; % 0.01 cpd, 30 deg radius circle
% dataFile = "rfGaborSimData_24-12-11_1323.mat"; % 0.02 cpd, 30 deg radius circle
% dataFile = "rfGaborSimData_24-12-11_1327.mat"; % 0.04 cpd, 30 deg radius circle
% dataFile = "rfGaborSimData_24-12-11_1331.mat"; % 0.08 cpd, 30 deg radius circle
% dataFile = "rfGaborSimData_24-12-11_1414.mat"; % 0.16 cpd, 30 deg radius circle
% dataFile = "rfGaborSimData_24-12-11_1418.mat"; % 0.32 cpd, 30 deg radius circle

loadFilename = strcat(dirstring,dataFile);
fprintf("Loading data from: %s\n",loadFilename)
load(loadFilename,"d","p")

% get RF responses averaged across stimulus phases
fprintf('Averaging responses across %0.00f stimulus phases... ',length(p.stimParams.phase))
rfRespPhaseAvg = avgResponsesAcrossStimPhase(p.rfParams,p.stimParams,d.responses);
fprintf('...done.\n\n')

% normalize and inject noise
fprintf('Rescaling stimulus-driven responses between 0-1, adding baseline response, and adding gaussian noise...\n')
rfRespNormIndVis = injectNoiseSimGUI(rfRespPhaseAvg,deltaFMax,noiseStd);
fprintf('done.\n\n')

% spatially smooth RF responses to simulate population measurements (in visual coords)
fprintf('Smoothing single unit RF responses in visual coordinates to simulate population responses...\n')
rfRespNormPopVis = spatialSmoothResponses(p,rfRespNormIndVis,smoothingSigmaVis);
fprintf('...done.\n\n')

% fix rotated response matrix (see printed statements)
fprintf('The code that generates responses rotates the coordinates by exactly 90 degrees...\n')
fprintf('Post-hoc undo for now since initial rotation doesn''t affect overall responses, but low priority fix later...\n')
fprintf('Rotating responses spatially by -90 degrees...\n')
rfRespNormIndVis = fixRotation(rfRespNormIndVis);
rfRespNormPopVis = fixRotation(rfRespNormPopVis);
fprintf('...done.\n\n')

% calculate orientation preferences, gOSIs, max responses in visual coordinates
fprintf('Calculating orientation preferences, gOSIs, and max responses for RF responses in visual coordinates...\n')
[oriPrefIndVis,~,gOSIIndVis,maxRespIndVis] = calcOSI(p.stimParams.theta,rfRespNormIndVis);
[oriPrefPopVis,~,gOSIPopVis,maxRespPopVis] = calcOSI(p.stimParams.theta,rfRespNormPopVis);
gOSINormIndVis = rescaleOSIs(gOSIIndVis);
gOSINormPopVis = rescaleOSIs(gOSIPopVis);
fprintf('...done.\n\n')

% generate anatomical coordinates and estimate responses in anat coords
tic
fprintf('Generating anatomical coordinates, estimating responses, and smoothing responses... ')
[rfRespIndAnat,rfRespPopAnat,mmScale] = smoothAnatGUI(rfRespNormIndVis,p.rfParams,anatSamples,sigmaAnatPix);
toc
fprintf('...done.\n\n')

% calculate orientation preferences, gOSIs, max responses in anatomical coordinates
fprintf('Calculating orientation preferences, gOSIs, and max responses for RF responses in anatomical coordinates...\n')
[oriPrefIndAnat,~,gOSIIndAnat,maxRespIndAnat] = calcOSI(p.stimParams.theta,rfRespIndAnat);
[oriPrefPopAnat,~,gOSIPopAnat,maxRespPopAnat] = calcOSI(p.stimParams.theta,rfRespPopAnat);
gOSINormIndAnat = rescaleOSIs(gOSIIndAnat);
gOSINormPopAnat = rescaleOSIs(gOSIPopAnat);
fprintf('...done.\n\n')


%% plotting
% plots the orientation preference maps in figure 4 (depending on loaded on stimulus file)
plotOriPrefAnat(oriPrefPopAnat,gOSINormPopAnat,anatSamples,gOSIThresh,alphaExp,mmScale,disableLabels)

savePlot = 0;
if savePlot
    figuresDir = "~/Desktop/";
    savefile = "popOriMap.png";
    fullsave = strcat(figuresDir,savefile);
    fprintf('Saving figure to %s...\n',fullsave)
    exportgraphics(gca,fullsave,'Resolution',450,"Padding","figure");
    fprintf('...done.\n\n')
end

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
        smoothedResps(stimOriIdx,:,:) = imgaussfilt(squeeze(unsmoothedResps(stimOriIdx,:,:)),smoothingSigma * rfGridScaleFactor,'Padding',0);
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

