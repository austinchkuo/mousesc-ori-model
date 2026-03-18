%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    reconstructFromHalfFourier    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [im,imf] = reconstructFromHalfFourier(d)

d.halfFourier = d.mag.*(cos(d.phase)+1i*sin(d.phase));

% first make the last column of the half fourier space which includes
% the dc and should have the frequency components replicated corectly
halfFourier = [d.halfFourier d.dc];
halfFourier(end+1:end+floor(d.originalDims(1)/2)) = conj(d.halfFourier(end:-1:end-floor(d.originalDims(1)/2)+1));
halfFourier = reshape(halfFourier,d.originalDims(1),ceil(d.originalDims(2)/2));

% replicate the frequency components to make the negative frequencies which
% are the complex conjugate of the positive frequncies
halfFourier2 = fliplr(flipud(conj(halfFourier(:,1:floor(d.originalDims(2)/2)))));
imf = [halfFourier halfFourier2];
im = ifft2(ifftshift([halfFourier halfFourier2]));