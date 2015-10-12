function [] = drawFixationBox(wPtr,xDispCenter,YCenter,boxSize,boxColour,boxFill)

if nargin<=4,
    boxColour=[255 255 255]; % the colour of the edge
end

lineWidth = 2;
fixLength = 5;
fixColour = [255 0 0];

% The fill:
if nargin>5
    Screen('FillRect', wPtr, boxFill, [xDispCenter-boxSize, ...
        YCenter-boxSize, xDispCenter+boxSize, YCenter+boxSize], lineWidth);
end
% The frame:
Screen('FrameRect', wPtr, boxColour, [xDispCenter-boxSize, ...
    YCenter-boxSize, xDispCenter+boxSize, YCenter+boxSize], lineWidth);
% The fixation cross:
Screen('DrawLine', wPtr, fixColour, xDispCenter-fixLength, ...
    YCenter ,xDispCenter+fixLength, YCenter, lineWidth);
Screen('DrawLine', wPtr, fixColour, xDispCenter, ...
    YCenter-fixLength, xDispCenter, YCenter+fixLength,lineWidth);

end