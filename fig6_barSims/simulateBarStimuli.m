function [response,stimParams,rfParams] = simulateBarStimuli(barOrientation,barType,boundType,stimParams,rfParams)

pixSize = 1200;
% set stimulus parameters
if ieNotDefined('stimParams')
    % defaults if not user-specified
    stimParams.size = pixSize;
    stimParams.barWidth = 2; % in degrees
    stimParams.stepSize = 1;
    stimParams.barLocExtent = 30;
    stimParams.barLocs = -stimParams.barLocExtent:stimParams.stepSize:stimParams.barLocExtent; % repeat for x and y
    stimParams.fovDeg = 60;
end
stimParams.barOrientations = barOrientation;
if ieNotDefined('barType'); barType = "boundary"; end
stimParams.barType = barType;
if ieNotDefined('boundType'); boundType = "edge"; end
stimParams.boundType = boundType;
if strcmp(stimParams.boundType,"box")
    stimParams.edgeOri = "box";
else
    stimParams.edgeOri = "vertical"; %"vertical", "box"
end

% set RF parameters
if ieNotDefined('rfParams')
    rfParams.size = pixSize;
    rfParams.fovDeg = 60;
    rfParams.SigmaCenter = 6; % 0.8;
    rfParams.SigmaSurround = 15;% 4.7;
    rfParams.surroundAmp = 0.4; % 0.4;
    rfParams.xMin = -20;
    rfParams.xMax = 20;
    rfParams.yMin = -20;
    rfParams.yMax = 20;
    rfParams.stepSize = 1;
    % define RF tiling
    nStepsX = abs(rfParams.xMin - rfParams.xMax)/rfParams.stepSize + 1;
    nStepsY = abs(rfParams.yMin - rfParams.yMax)/rfParams.stepSize + 1;
    rfParams.xCenter = linspace(rfParams.xMin,rfParams.xMax,nStepsX);
    rfParams.yCenter = linspace(rfParams.yMin,rfParams.yMax,nStepsY);
end

parallelize = 1;
verbose = 1;
response = guiRFStimRespWrapper(stimParams,rfParams,parallelize,verbose);

%% helper functions
function [responseP1] = guiRFStimRespWrapper(stimParams,rfParams,parallelize,verbose)

responseP1 = nan(length(stimParams.barOrientations),length(stimParams.barLocs),length(rfParams.xCenter),length(rfParams.yCenter));
responseP2 = nan(length(stimParams.barOrientations),length(stimParams.barLocs),length(rfParams.xCenter),length(rfParams.yCenter));
complex_response = nan(length(stimParams.barOrientations),length(stimParams.barLocs),length(rfParams.xCenter),length(rfParams.yCenter));

totalTimeStart = tic;
for xIdx = 1:length(rfParams.xCenter)
    if verbose
        fprintf('\nGenerating RFs and responses - x location %d of %d (x = %0.2f deg)...',xIdx,length(rfParams.xCenter),rfParams.xCenter(xIdx))
    end
    for yIdx = 1:length(rfParams.yCenter)
            % space
            tempRF = generateRFSpace(rfParams.size,rfParams.fovDeg,[rfParams.SigmaCenter rfParams.SigmaSurround],...
                                     rfParams.xCenter(xIdx),rfParams.yCenter(yIdx),rfParams.surroundAmp); % if you error here, check that your FOV is large enough for your RF tiling
            rfMat_P1{yIdx} = tempRF.combP1;
            rfMat_P2{yIdx} = tempRF.combP2;

    end

    % generate stimuli and calculate responses
    [responseP1(:,:,xIdx,:),responseP2(:,:,xIdx,:),complex_response(:,:,xIdx,:)] = calcRFRespLoop(stimParams,rfParams,rfMat_P1,rfMat_P2,parallelize);
    
    % calculate progress
    percentDone = xIdx/length(rfParams.xCenter)*100;
    totalTimeEnd = toc(totalTimeStart);
    tempTimeStr = seconds(totalTimeEnd);
    tempTimeStr.Format = 'hh:mm:ss';
    estimatedTotalTime = seconds(100*totalTimeEnd/percentDone);
    estimatedTotalTime.Format = 'hh:mm:ss';
    timeLeft = estimatedTotalTime - tempTimeStr;
    timeLeft.Format = 'hh:mm:ss';
    if verbose
        fprintf('\n(guiRFStimRespWrapper) Total time elapsed (%0.2f%%): %s\n(guiRFStimRespWrapper) Estimated total time: %s. Estimated time to completion: %s\n',percentDone,char(tempTimeStr),char(estimatedTotalTime),char(timeLeft))
    else
        fprintf('\n(guiRFStimRespWrapper) Time elapsed (%0.2f%%): %s',percentDone,char(tempTimeStr))
    end
    
end
end

function [responsesP1,responsesP2,complex_responses] = calcRFRespLoop(stimParams,rfParams,rfsP1,rfsP2,parallelize)

if ~exist('splitByColumns','var'); splitByColumns = 1; end

verbose = 2;
stimOris = stimParams.barOrientations;
stimLocs = stimParams.barLocs;
rfXCenters = rfParams.xCenter;
rfYCenters = rfParams.yCenter;

nStimOris = length(stimOris);
nStimLocs = length(stimLocs);
nRFXLocs = length(rfXCenters);
nRFYLocs = length(rfYCenters);

if splitByColumns
    nRFXLocs = 1;
end

% must preallocate arrays for parallelization
if verbose > 0
    fprintf('\n(calcRFRespLoop) Generating stimuli and calculating responses to %0.0f RFs...',nRFXLocs*nRFYLocs)
end
% Preallocate arrays for center-surround RFs
responsesP1 = nan(nStimOris, nStimLocs, nRFXLocs, nRFYLocs);
responsesP2 = nan(nStimOris, nStimLocs, nRFXLocs, nRFYLocs);
complex_responses = nan(nStimOris, nStimLocs, nRFXLocs, nRFYLocs);

p = gcp('nocreate');
if isempty(p)
    parpool('threads');
end

% Total number of iterations needed for the linear index
totalIterations = nStimOris * nStimLocs;

if verbose > 1
    tic
    fprintf('\n(calcRFRespLoop) Calculating responses for %0.0f (%0.0f oris x %0.0f locations) stimuli... ',totalIterations,nStimOris,nStimLocs)
end
parfor linearIdx = 1:totalIterations % parfor
    % Convert linear index back to sub indices
    [stimOriIdx, stimLocIdx] = ind2sub([nStimOris,nStimLocs],linearIdx);

    % Generate stimuli and calculate response, store in temporary cell arrays
    [tempRespP1{linearIdx},tempRespP2{linearIdx},tempComplexResp{linearIdx}] = genStimCalcRespBar(stimParams,rfParams,rfsP1,rfsP2,stimOriIdx,stimLocIdx);
end

% Reassign to the main output arrays
for linearIdx = 1:totalIterations
    [stimOriIdx, stimLocIdx] = ind2sub([nStimOris, nStimLocs], linearIdx);
    responsesP1(stimOriIdx,stimLocIdx,:,:) = tempRespP1{linearIdx};
    responsesP2(stimOriIdx,stimLocIdx,:,:) = tempRespP2{linearIdx};
    complex_responses(stimOriIdx,stimLocIdx,:,:) = tempComplexResp{linearIdx};
end
if verbose > 1
    toc
end
end

end % end simulateBarStimuli()










