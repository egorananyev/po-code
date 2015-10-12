function [ px ] = cm2px ( cm, sdims )

% screenXcm
% screenYcm
% screenXpx
% screenYpx
% screenDist

% px = (sdims.screenXpx * cm) / sdims.screenXcm;
px = (c2n(sdims,'screenXpx',1) * cm) / c2n(sdims,'screenXcm',1);

end