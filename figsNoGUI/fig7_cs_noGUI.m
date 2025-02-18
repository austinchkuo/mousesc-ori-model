% fig7_cs_noGUI.m

% RF parameters:
% C/S size ratio: 0.8 to 4.7 deg (resulting center frequency: 0.08 cpd)
% surround amplitude: 0.4

clear
close all

% simulated response parameters
deltaFMax = 1;
noiseStd = 0.01;

% visualization parameters
gOSIThresh = 0;
alphaExp = 0.00001;
disableLabels = 1;
customFile = "~/Documents/MATLAB/population-model/rfSimulations/anatSamples/plot-data.csv"; % specifies spatial locations of units matching Fig. 4 from Liang et al., 2023
customStarts = [0.81 0.68]; % starting x and y values of the custom RF locations

% specify data file
dirstring = "~/Documents/MATLAB/data/rfSimData/";
dataFile = "rfGaborSimData_25-02-07_1040.mat";
% dataFile = "rfGaborSimData_25-02-07_1040.mat"; % 0.04 cpd, vertical edge
% dataFile = "rfGaborSimData_25-02-07_1045.mat"; % 0.04 cpd, horizontal edge
% dataFile = "rfGaborSimData_25-02-07_1141.mat"; % 0.04 cpd, full field
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

% fix rotated response matrix (see printed statements)
fprintf('The code that generates responses rotates the coordinates by exactly 90 degrees...\n')
fprintf('Post-hoc undo for now since initial rotation doesn''t affect overall responses, but low priority fix later...\n')
fprintf('Rotating responses spatially by -90 degrees...\n')
rfRespNormIndVis = fixRotation(rfRespNormIndVis);
fprintf('...done.\n\n')

% generate anatomical coordinates and estimate responses in anat coords
tic
fprintf('Generating anatomical coordinates and estimating responses...\n')
[rfRespNormIndAnat,~,mmScale] = smoothAnatGUI(rfRespNormIndVis,p.rfParams,[],[],"nonuniform",customFile,customStarts);
toc
fprintf('...done.\n\n')

% calculate orientation preferences, gOSIs, max responses in anatomical coordinates
fprintf('Calculating orientation preferences, gOSIs, and max responses for RF responses in anatomical coordinates...\n')
[oriPrefIndAnat,~,gOSIIndAnat,maxRespIndAnat] = calcOSI(p.stimParams.theta,rfRespNormIndAnat);
gOSINormIndAnat = rescaleOSIs(gOSIIndAnat);
fprintf('...done.\n\n')

plotOriPrefAnat(oriPrefIndAnat,gOSINormIndAnat,[],gOSIThresh,alphaExp,mmScale,disableLabels,customFile,customStarts)

%% %%%%%%%%%%%%%%%%%% %%
%%% helper functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%
function cr_phase_avg = avgResponsesAcrossStimPhase(rfParams,stimParams,responses)
    % complex responses averaged across stimulus phases
    reshapeSize = nonzeros([length(stimParams.theta),length(rfParams.xCenter),length(rfParams.yCenter),length(rfParams.theta)])';
    cr_phase_avg = reshape(mean(responses.complex_response,2),reshapeSize);

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

function [] = plotOriPrefAnat(oriPrefs,OSIs,anatSamples,osiThreshold,alphaExponent,mmScale,disableLabels,customFile,customStarts)

    % x and y coordinates are given by the file specified
    temp = readtable(customFile);
    temp.x = temp.x/2;
    temp.y = temp.y/2;
    xLocs = temp.x + (customStarts(1) - min(temp.x));
    yLocs = temp.y + (customStarts(2) - min(temp.y));
    liangNonOS = zeros(1,53); liangNonOS(4) = 1; liangNonOS(22) = 1; liangNonOS(29) = 1; liangNonOS(32) = 1; liangNonOS(41) = 1;
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
    s_oriPref.SizeData = 100;
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