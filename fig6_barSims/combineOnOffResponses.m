% combineOnOffResponses.m - threshold and sum on-off RF phase responses

function complexResp = combineOnOffResponses(respP1,respP2,threshold)

% default: threshold opposing responses instead of squaring and summing
if ~exist('threshold','var'); threshold = 1; end

if threshold
    respP1(respP1<0) = 0;
    respP2(respP2<0) = 0;
    complexResp = respP1 + respP2;
else
    complexResp = respP1^2 + respP2^2;
end

