% shiftImage2D - shift image in 2 dimensions, given an x and y shift (units in indices)
% 
% behaves like circshift, but without any wrapping (any new pixels are given a value of 0)
% 
% helpfully, this deals with the indexing garbage that is soooo annoying to deal with under the hood
% 
% you do not have to worry about the fact that y indices and y coordinates are flipped (i.e., low to high index goes from top to bottom, while low to high coordinates goes from bottom to top)
% you do not have to worry about the fact that a shift in x coordinates affects the COLUMNS, while a shift in y coordinates affects the ROWS
% simply input what you intuitively think in regular cartesian coordinate space what a shift in x/y should be (negative for left/down, positive for up/right)
% 
% naturally, this means that if you are thinking in matrix coordinate space, do not use this function... try using a MATLAB-provided function directly instead
% 
% AK - 5/14/24

function xyShiftedMat = shiftImage2D(origMat,xShift,yShift)

if xShift < 0 % shift left

    if yShift < 0 % shift down
        yShiftedMat = paddata(origMat(1:end+yShift,1:end),size(origMat),Side="leading");
    else % shift up
        yShiftedMat = paddata(origMat(1+yShift:end,1:end),size(origMat));
    end

    xyShiftedMat = paddata(yShiftedMat(1:end,1-xShift:end),size(yShiftedMat));

else % shift right

    if yShift < 0 % shift down
        yShiftedMat = paddata(origMat(1:end+yShift,1:end),size(origMat),Side="leading");
    else % shift up
        yShiftedMat = paddata(origMat(1+yShift:end,1:end),size(origMat));
    end

    xyShiftedMat = paddata(yShiftedMat(1:end,1:end-xShift),size(yShiftedMat),Side="leading");

end