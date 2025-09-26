% genStimCalcResp.m - wrapper function for generating a set of grating stimuli and calculating a set of RF responses to generated stimuli

function [responseP1, responseP2, complex_response] = genStimCalcRespBar(stimParams,rfParams,rfsP1,rfsP2,stimOriIdx,stimLocIdx,splitByColumns)

if ~exist('splitByColumns','var'); splitByColumns = 1; end

stimOris = stimParams.barOrientations;
stimLocs = stimParams.barLocs;
stimWidth = stimParams.barWidth;
fovPix = stimParams.size;
fovDeg = stimParams.fovDeg;
rfXCenters = rfParams.xCenter;
rfYCenters = rfParams.yCenter;

nRFXLocs = length(rfXCenters);
nRFYLocs = length(rfYCenters);

if splitByColumns
    nRFXLocs = 1;
end

responseP1 = zeros(nRFXLocs, nRFYLocs);
responseP2 = zeros(nRFXLocs, nRFYLocs);
complex_response = zeros(nRFXLocs, nRFYLocs);

% generate stimuli (space)
stimSpace = zeros(fovPix + 1);
% generate bar stimuli across locations
pixPerDeg = fovPix/fovDeg;
barWidthPix = stimWidth * pixPerDeg;
zeroIdx = 1+round((fovDeg/2)*pixPerDeg);
if strcmp(stimParams.boundType,"edge")
    if strcmp(stimParams.edgeOri,"vertical")
        stimSpace(:,round(zeroIdx):end) = 1;
        edgeStartIdx = [1,zeroIdx]; % [row,col], not [x,y]
    elseif strcmp(stimParams.edgeOri,"horizontal")
        stimSpace(round(zeroIdx):end,:) = 1;
        edgeStartIdx = [zeroIdx,1]; % [row,col], not [x,y]
    end
elseif strcmp(stimParams.boundType,"corner")
    stimSpace(round(zeroIdx):end,round(zeroIdx):end) = 1;
    edgeStartIdx = [zeroIdx,zeroIdx];
elseif strcmp(stimParams.boundType,"box")
    stimSpace(round(zeroIdx/2):round(3*zeroIdx/2),round(zeroIdx/2):round(3*zeroIdx/2)) = 1;
    edgeStartIdx = [zeroIdx/2,zeroIdx/2]; % [row,col], not [x,y]
    edgeEndIdx = [3*zeroIdx/2,3*zeroIdx/2]; % [row,col], not [x,y]
end
% thisOri = stimOris(stimOriIdx);
thisLoc = stimLocs(stimLocIdx);
barSpace = zeros(fovPix + 1);

if strcmp(stimOris(stimOriIdx),"horizontal")
    barCenterPix = zeroIdx + -(thisLoc*pixPerDeg); % because matlab uses the top left corner as (0,0) and positive y is down
    barLimsSpace = [barCenterPix - barWidthPix/2, barCenterPix + barWidthPix/2];
    barFlag = 1;
    % crop for bar going out of bounds
    if barLimsSpace(1) > size(stimSpace,1)
        barLimsSpace(1) = size(stimSpace,1);
    end
    if barLimsSpace(2) > size(stimSpace,2)
        barLimsSpace(2) = size(stimSpace,2);
    end
    if barLimsSpace(1) < 1
        barLimsSpace(1) = 1;
    end
    if barLimsSpace(2) < 1
        barLimsSpace(2) = 1;
    end
    if strcmp(stimParams.barType,"boundary")
        % do not draw a bar unless bar is fully in bounds of the "display"
        if strcmp(stimParams.edgeOri,"box")
            if barLimsSpace(1) < edgeStartIdx(1) || barLimsSpace(2) > edgeEndIdx(1)
                barFlag = 0;
            end
        elseif barLimsSpace(1) < edgeStartIdx(1)
            barFlag = 0;
        end
    end
    if barFlag
        barSpace(round(min(barLimsSpace)):round(max(barLimsSpace)),:) = 1; % make a horizontal bar
    end

elseif strcmp(stimOris(stimOriIdx),"vertical")
    barCenterPix = zeroIdx + thisLoc*pixPerDeg; % because matlab uses the top left corner as (0,0) and positive x is right
    barLimsSpace = [barCenterPix - barWidthPix/2, barCenterPix + barWidthPix/2];
    barFlag = 1;
    % crop for bar going out of bounds
    if barLimsSpace(1) > size(stimSpace,1)
        barLimsSpace(1) = size(stimSpace,1);
    end
    if barLimsSpace(2) > size(stimSpace,2)
        barLimsSpace(2) = size(stimSpace,2);
    end
    if barLimsSpace(1) < 1
        barLimsSpace(1) = 1;
    end
    if barLimsSpace(2) < 1
        barLimsSpace(2) = 1;
    end
    if strcmp(stimParams.barType,"boundary")
        % do not draw a bar unless bar is fully in bounds of the "display"
        if strcmp(stimParams.edgeOri,"box")
            if barLimsSpace(1) < edgeStartIdx(2) || barLimsSpace(2) > edgeEndIdx(2)
                barFlag = 0;
            end
        elseif barLimsSpace(1) < edgeStartIdx(2)
            barFlag = 0;    
        end
    end
    if barFlag
        barSpace(:,round(barLimsSpace(1)):round(barLimsSpace(2))) = 1; % make a vertical bar
    end

end
finalStimSpace = stimSpace.*barSpace;
if strcmp(stimParams.edgeOri,"vertical")
    finalStimSpace(:,1:round(zeroIdx)-1) = -1;
    % calculate baseline
    noStimSpace = finalStimSpace;
    noStimSpace(:,zeroIdx:end) = 0;
elseif strcmp(stimParams.edgeOri,"horizontal")
    finalStimSpace(1:round(zeroIdx)-1,:) = -1;
    % calculate baseline
    noStimSpace = finalStimSpace;
    noStimSpace(zeroIdx:end,:) = 0;
elseif strcmp(stimParams.edgeOri,"box")
    % vertical blank spaces
    finalStimSpace(:,1:round(zeroIdx/2)-1) = -1;
    finalStimSpace(:,round(3*zeroIdx/2)+1:end) = -1;
    % horizontal blank spaces
    finalStimSpace(1:round(zeroIdx/2)-1,:) = -1;
    finalStimSpace(round(3*zeroIdx/2)+1:end,:) = -1;
    noStimSpace = stimSpace-1;
end

% loop through all RFs and calculate responses
for rfPosIdx_x = 1:nRFXLocs
    for rfPosIdx_y = 1:nRFYLocs
        %{
        figure(2)
        subplot(1,2,1)
        imagesc(stimSpace)
        axis square
        subplot(1,2,2)
        imagesc(rfsP1{1,1})
        axis square
        pause(0.05)
        %}
        % fprintf('\nBar lims space: %0.2f:%0.2f',round(barLimsSpace(1)),round(barLimsSpace(2)))
        baselineP1 = dot(rfsP1{rfPosIdx_x, rfPosIdx_y}(:), noStimSpace(:));
        baselineP2 = dot(rfsP2{rfPosIdx_x, rfPosIdx_y}(:), noStimSpace(:));
        responseP1(rfPosIdx_x, rfPosIdx_y) = dot(rfsP1{rfPosIdx_x, rfPosIdx_y}(:), finalStimSpace(:)) - baselineP1;
        responseP2(rfPosIdx_x, rfPosIdx_y) = dot(rfsP2{rfPosIdx_x, rfPosIdx_y}(:), finalStimSpace(:)) - baselineP2;
        complex_response(rfPosIdx_x, rfPosIdx_y) = combineOnOffResponses(responseP1(rfPosIdx_x, rfPosIdx_y),responseP2(rfPosIdx_x, rfPosIdx_y),1);

    end
end
