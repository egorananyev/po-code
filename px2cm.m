function [ cm ] = px2cm( px, sdims )

% sdims.screenXcm
% sdims.screenYcm
% sdims.screenXpx
% sdims.screenYpx
% sdims.screenDist

% cm = (sdims.screenXcm * px) / sdims.screenXpx;
cm = (c2n(sdims,'screenXcm',1) * px) / c2n(sdims,'screenXpx',1);

end