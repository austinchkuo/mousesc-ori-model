% simBarScript.m - script to simulate bar stimuli and calculate RF responses
%                  can also apply anatomical transformation to RF responses
%                  used to generate plots from Fig. 8 of Kuo et al., 2025
% 
% Usage - run the script with the appropriate parameters to generate plots from Fig. 8
% 
% adjustable parameters: boundType - "edge", "corner", "box"
%                            edge: draws a straight edge separating simulated "on-screen" region and "off-screen" region
%                                  the edge is centered along either the x- or y-axis, dependent on the bar stimulus being 
%                                  horizontally or vertically oriented 
%                            corner: draws a corner separating simulated "on-screen" region and "off-screen" region
%                                    the corner point is centered at (0,0) in cartesian coordinates (not matrix indices), 
%                                    making it so that 3/4 of the screen is considered "off-screen" and the remaining 1/4 
%                                    is considered "on-screen"
%                            box: draws a square "on-screen" region surrounded by a larger square "off-screen" region
%                                 the "on-screen" square has half the diameter of the "off-screen" region
%                                 both square regions are centered at (0,0) in cartesian coordinates of the FOV
%                        barType - "mask", "boundary"
%                            mask: treats the simulated display edges as a mask (see manuscript methods for details)
%                            boundary: treats the simulated display edges as a boundary (see manuscript methods for details)
%                        stimParams - struct that contains the following fields
%                            size: FOV in pixels
%                            barWidth: width of bar in degrees
%                            stepSize: how far the bar moves each timestep along the FOV in degrees
%                            barLocExtent: extent that the bar travels across the FOV (in degrees, relative to center = (0,0))
%                            fovDeg: size of the FOV, in degrees
%                        rfParams - struct that contains the following fields
%                            size: FOV in pixels (should be kept the same as stimParams.size)
%                            fovDeg: size of the FOV, in degrees (should be kept the same as stimParams.fovDeg)
%                            SigmaCenter: std dev of the RF center, in degrees
%                            SigmaSurround: std dev of the RF surround, in degrees
%                            surroundAmp: relative ratio of volumes of the 2d surround gaussian to the 2d center gaussian
%                            xMin/xMax: extent of RF center locations (x-coord), in degrees, relative to the FOV center (0,0)
%                            yMin/yMax: extent of RF center locations (y-coord), in degrees, relative to the FOV center (0,0)
%                            stepSize: how far apart RF centers are spaced along extent, in degrees

clear
close all

boundType = "box"; % "edge", "corner", "box"
if strcmp(boundType,"box")
    barType = "boundary";
else
    barType = ["mask","boundary"];
end

pixSize = 1200;
% set stimulus parameters
stimParams.size = pixSize;
stimParams.barWidth = 2; % in degrees
stimParams.stepSize = 1;
stimParams.barLocExtent = 30;
stimParams.barLocs = -stimParams.barLocExtent:stimParams.stepSize:stimParams.barLocExtent; % repeat for x and y
stimParams.fovDeg = 60;

% set RF parameters
rfParams.size = pixSize;
rfParams.fovDeg = 60;
rfParams.SigmaCenter = 0.8; 
rfParams.SigmaSurround = 4.7;
rfParams.surroundAmp = 0.4;
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

figLabels = 0;
useRelu = 1;

if strcmp(boundType,"box")
    anatXform = 1;
else
    anatXform = 0;
end

for i = 1:length(barType)
    [respH,stimParamsH,rfParamsH] = simulateBarStimuli("horizontal",barType(i),boundType,stimParams,rfParams);
    [respV,stimParamsV,rfParamsV] = simulateBarStimuli("vertical",barType(i),boundType,stimParams,rfParams);
    midRFHx = round(length(rfParamsH.xCenter)/2);
    midRFHy = round(length(rfParamsH.yCenter)/2);
    midRFVx = round(length(rfParamsV.xCenter)/2);
    midRFVy = round(length(rfParamsV.yCenter)/2);
    if midRFHx ~= midRFVx || midRFHy ~= midRFVy
        fprinf('\n(simBarScript) The middle RF for vertical and horizontal stimuli are not the same!! Is this intended??\n')
    end
    respMidH = respH(1,:,midRFHx,midRFHy);
    respMidV = respV(1,:,midRFVx,midRFVy);
    if useRelu
        respMidH_noRelu = respMidH;
        respMidV_noRelu = respMidV;
        respMidH(respMidH<0) = 0;
        respMidV(respMidV<0) = 0;
    end

    % rescale so that the numbers are in normalized and not arbitrary units
    meanHorzMid = mean(respMidH);
    meanVertMid = mean(respMidV);
    nrespMidH = respMidH / max(abs([respMidH respMidV]));
    nrespMidV = respMidV / max(abs([respMidH respMidV]));
    nmeanHorzMid = mean(nrespMidH);
    nmeanVertMid = mean(nrespMidV);
    barHorzOut{i} = meanHorzMid;
    barVertOut{i} = meanVertMid;
    
    if useRelu
        % rescale so that the numbers are in normalized and not arbitrary units
        nrespMidH_noRelu = respMidH_noRelu / max(abs([respMidH_noRelu respMidV_noRelu]));
        nrespMidV_noRelu = respMidV_noRelu / max(abs([respMidH_noRelu respMidV_noRelu]));
        nmeanHorzMid_noRelu = mean(nrespMidH_noRelu);
        nmeanVertMid_noRelu = mean(nrespMidV_noRelu);
    end
    
    if strcmp(barType(i),"mask")
        barCentersH = stimParamsH.barLocs;
        barCentersV = stimParamsV.barLocs;
    elseif strcmp(barType(i),"boundary")
        barCentersH = stimParamsH.barLocs;
        barCentersV = barCentersH;
    end
    figure(i)
    plot(barCentersH,nrespMidH,"linewidth",2,"color",[0 0 1])
    hold on
    plot(barCentersV,nrespMidV,"linewidth",2,"color",[1 0 0])
    if useRelu % plot non-rectified responses
        plot(barCentersH,nrespMidH_noRelu,"linewidth",2,"color",[0 0 1 0.2],"linestyle","--")
        plot(barCentersV,nrespMidV_noRelu,"linewidth",2,"color",[1 0 0 0.2],"linestyle","--")
    end
    if figLabels
        xlabel('Bar center along screen (deg)')
        ylabel('RF response (normalized)')
        legend({"Horizontal bar","Vertical bar"})
    end
    axis square
end

if length(barType) > 1
    i = i+1;
    figure(i)
    % rescale means to max 1
    barResps = [barHorzOut{1} barHorzOut{2}; barVertOut{1} barVertOut{2}];
    nbarResps = barResps ./ max(max(barResps));
    if figLabels
        b = bar(["Mask","Boundary"],[nbarResps(1,1),nbarResps(2,1); nbarResps(1,2),nbarResps(2,2)]);
        title(sprintf("Horizontal to vertical bar response ratio: %0.2f%%, %0.2f%%",abs(barHorzOut{1}/barVertOut{1})*100,abs(barHorzOut{2}/barVertOut{2})*100))
        xlabel('Edge behavior')
        ylabel('Mean response across bar positions (normalized)')
        legend({"Horizontal","Vertical"})
        axis square
    else
        bar([nbarResps(1,1),nbarResps(2,1); nbarResps(1,2),nbarResps(2,2)]);
        xticklabels([]);
        yticklabels([]);
        ylim([-0.2 1.1])
        yticks([-1:0.2:1])
        axis square
        i = i+1;
        figure(i)
        b = bar(["Mask","Boundary"],[nbarResps(1,1),nbarResps(2,1); nbarResps(1,2),nbarResps(2,2)]);
        ylim([-0.2 1.1])
        title(sprintf("Horizontal to vertical bar response ratio: %0.2f%%, %0.2f%%",abs(barHorzOut{1}/barVertOut{1})*100,abs(barHorzOut{2}/barVertOut{2})*100))
        xlabel('Edge behavior')
        ylabel('Mean response across bar positions (normalized)')
        legend({"Horizontal","Vertical"})
        axis square
    end
end

% make maps of responses of grid responses
if size(respH,3) > 1
    meanh = mean(squeeze(respH));
    meanv = mean(squeeze(respV));
    
    if length(size(meanh)) > 2
        hv(1,:,:) = squeeze(meanh);
        hv(2,:,:) = squeeze(meanv); 
    else % do not squeeze if the matrix is already 2-dim (i.e., single row of RFs), but transpose to make dims: RFs x location
        hv(1,:,:) = meanh';
        hv(2,:,:) = meanv';
    end

    % do the same rotation fix as for the grating stimuli
    tempResps = permute(hv(:,:,:),[2,3,1]);
    tempResps = rot90(tempResps,-1);
    hv = permute(tempResps,[3,1,2]);

    hvResp = reshape(hv,[2,size(hv,2)*size(hv,3)]);
    if useRelu
        hvResp(hvResp<0) = 0;
    end
    hvPref = zeros(size(hvResp,2)*size(hvResp,3),1);

    for j = 1:size(hvResp,2)
        if hvResp(1,j) > hvResp(2,j)
            hvPref(j) = 1;
        elseif hvResp(1,j) < hvResp(2,j)
            hvPref(j) = -1;
        end
    end
    hvSel = abs(hvResp(1,:)-hvResp(2,:))./(abs(hvResp(1,:))+abs(hvResp(2,:))); 
    hvSel(isnan(hvSel)) = 0;
    
    i = i+1;
    figure(i)
    [X,Y] = meshgrid(rfParamsH.xCenter,rfParamsH.yCenter);
    s = scatter(X(:),Y(:),50*ones(length(X(:)),1),hvPref(:),'filled','MarkerFaceAlpha','flat');
    s.AlphaData = hvSel'.^1;
    vCol = [71 204 204];
    hCol = [242 85 98];
    prefCmap = [vCol; hCol]/255;
    colormap(prefCmap)
    axis square
end

if anatXform
    % set up inputs for smoothAnatBars
    respMat = reshape(hvResp,size(hv));
    rfParams.xCenter = rfParamsH.xCenter;
    rfParams.yCenter = rfParamsH.yCenter;
    anatSampling = 81;
    smoothSigma = 1.5; % 0.0083;

    % apply anatomical transformation to the RFs
    [hvaResp,hvaRespSmoothed] = smoothAnatBars(respMat,rfParams,anatSampling,smoothSigma);
    hvaPref = zeros(size(hvaResp,2)*size(hvaResp,3),1);
    hvaResp = reshape(hvaResp,[2 size(hvaResp,2)*size(hvaResp,3)]);
    for j = 1:size(hvaResp,2)
        if hvaResp(1,j) > hvaResp(2,j)
            hvaPref(j) = 1;
        elseif hvaResp(1,j) < hvaResp(2,j)
            hvaPref(j) = -1;
        end
    end
    hvaSel = abs(hvaResp(1,:)-hvaResp(2,:))./(abs(hvaResp(1,:))+abs(hvaResp(2,:))); 
    hvaSel(isnan(hvaSel)) = 0;

    i = i+1;
    figure(i)
    % translate anatomical coordinates to mm scale
    xAnats = linspace(0.4125,1.2375,anatSampling);
    yAnats = linspace(1.3375,0.5125,anatSampling);
    [X,Y] = meshgrid(xAnats,yAnats);
    s = scatter(X(:),Y(:),10*ones(length(X(:)),1),hvaPref(:),'filled','MarkerFaceAlpha','flat');
    s.AlphaData = hvaSel'.^1;
    colormap(prefCmap)
    axis square

end