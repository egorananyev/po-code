function equil3
% Simple version of flicker photometry

% BASIC SETUPS --------------------------------------------------
%clear all;
%close all;
%sca;
%commandwindow;
%HideCursor;
%PsychDefaultSetup(1);


% OPEN WINDOW ---------------------------------------------------
screenNo = max(Screen('Screens'));
white = WhiteIndex(screenNo);
black = BlackIndex(screenNo);
grey = white/2;
fixPoint = 4;

% Keyboard
KbName('UnifyKeyNames');
targetUsageName = 'Keyboard'; % change accordingly
targetProduct = 'Dell USB Keyboard'; % change accordingly
% targetProduct = 'Apple Keyboard'; % temp
dev = PsychHID('Devices');
deviceIndex = find(strcmpi(targetUsageName, {dev.usageName}) & ...
    strcmpi(targetProduct, {dev.product}));
% deviceIndex = -3;
KbQueueCreate(deviceIndex);
KbQueueStart(deviceIndex);

% numMultiSamples = 4;
% [window, windowRect] = PsychImaging('OpenWindow', ...
%     screenNo, black, [0 0 1920 1080], 32, 2,[], numMultiSamples, []);
[window, windowRect] = Screen('OpenWindow', screenNo, black);
Screen('Flip', window);
flipInterval = Screen('GetFlipInterval',window);
[xCenter, yCenter] = RectCenter(windowRect);
screenXpixels = windowRect(3);
screenYpixels = windowRect(4);
Screen('TextSize', window, 40);
% Blend function, anti-aliasing
Screen('BlendFunction', window, ...
    'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');


escapeKey = KbName('ESCAPE');
fastKey = KbName('RightArrow');
slowKey = KbName('LeftArrow');
equalKey = KbName('DownArrow');

%TIME -----------------------------------------------------------
durFix = 0.5;
delay = flipInterval/2;


% TEST BUFFER ---------------------------------------------------
disp('Test begin')
DrawFormattedText(window,'Trials\n\nPress Any Key to Start','center','center',grey);
Screen('Flip', window);
KbStrokeWait(deviceIndex);

%Initial rects
red = [round(255*0.05) 0 0];
green = [round(255*0.05*.89) round(255*0.05) 0];
% green = [round(165*.89) round(165*.11) 0];
% green = [round(165*.78) round(165*.22) 0];
% green = [round(165*.53) round(165*.47) 0];
baseRect = [0 0 100 100];
allRects= CenterRectOnPointd(baseRect,xCenter,yCenter);
loopFlag = 1;
vbl = Screen('Flip',window);
%----------------------------------------------------------------
%                        EXPERIMENT LOOP
%----------------------------------------------------------------
while loopFlag
    [keyIsDown, secs, keyCode] = KbCheck(deviceIndex);
    for m1 = 1:4
        Screen('FillRect', window, green, allRects);
        vbl = Screen('Flip',window,vbl+delay);
    end
    
    %present reference object
    for m1 = 1:4
        Screen('FillRect', window, red, allRects);
        vbl = Screen('Flip',window,vbl+delay);
    end
    

    if keyCode(escapeKey)||keyCode(slowKey)||keyCode(fastKey)||keyCode(equalKey)
        if keyCode(escapeKey)
            ShowCursor;
            disp('Experiment aborted');
            sca;
            return;
        elseif keyCode(slowKey)
            green = green - [0 1 0];
        elseif keyCode(fastKey)
            green = green + [0 1 0] ;
        elseif keyCode (equalKey)
            loopFlag = 0;
            disp('E');
            disp(red);
            disp(green);
            break;
        end
    end
    %Keyboard check ends
end


% CLOSE --------------------------------------------------------
% save(matDataName);

DrawFormattedText(window,...
    'Press Any Key To Exit',...
    'center', 'center', grey);
Screen('Flip', window);
KbStrokeWait(deviceIndex);
sca;
disp('Test ends.')

%===============================END=============================
end
