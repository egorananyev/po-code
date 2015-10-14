function [] = drawMondrians(img,wPtr,xDispCenter,YCenter,boxSize)

mondtex = Screen('MakeTexture', wPtr, img);

Screen('DrawTexture', wPtr, mondtex, [], ...
    [xDispCenter-boxSize, YCenter-boxSize, xDispCenter, YCenter])
Screen('DrawTexture', wPtr, mondtex, [], ...
    [xDispCenter, YCenter-boxSize, xDispCenter+boxSize, YCenter])
Screen('DrawTexture', wPtr, mondtex, [], ...
    [xDispCenter-boxSize, YCenter, xDispCenter, YCenter+boxSize])
Screen('DrawTexture', wPtr, mondtex, [], ...
    [xDispCenter, YCenter, xDispCenter+boxSize, YCenter+boxSize])

end