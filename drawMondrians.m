function [] = drawMondrians(img,wPtr,xDispCenter,YCenter,boxSize,maskAlpha)

mondtex = Screen('MakeTexture', wPtr, img);

Screen('DrawTexture', wPtr, mondtex, [], ...
    [xDispCenter-boxSize, YCenter-boxSize, xDispCenter, YCenter], 0, [], maskAlpha)
Screen('DrawTexture', wPtr, mondtex, [], ...
    [xDispCenter, YCenter-boxSize, xDispCenter+boxSize, YCenter], 90, [], maskAlpha)
Screen('DrawTexture', wPtr, mondtex, [], ...
    [xDispCenter-boxSize, YCenter, xDispCenter, YCenter+boxSize], 180, [], maskAlpha)
Screen('DrawTexture', wPtr, mondtex, [], ...
    [xDispCenter, YCenter, xDispCenter+boxSize, YCenter+boxSize], 270, [], maskAlpha)

end