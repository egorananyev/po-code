function flickerV4(grbVal)

% Flicker photometry, sine wave flicker

% BASIC SETUPS --------------------------------------------------
% clear all;
% close all;
% sca;
% HideCursor;
% commandwindow;
% PsychDefaultSetup(1);

% Screen('Preference', 'SkipSyncTests', 1); % a necessary evil - only for Yosemite

% OPEN WINDOW ---------------------------------------------------
Screen('Screens')
screenNo = min(Screen('Screens'));
white = WhiteIndex(screenNo);
black = BlackIndex(screenNo);
grey = white/2;
fixPoint = 4;

% Keyboard
KbName('UnifyKeyNames');
quitkey = 'c';
deviceIndex = -3;
KbQueueCreate(deviceIndex);
KbQueueStart(deviceIndex);

% numMultiSamples = 4;
% [window, windowRect] = PsychImaging('OpenWindow', ...
%     screenNo, black, [0 0 1920 1080], 32, 2,[], numMultiSamples, []);
[window, windowRect] = Screen('OpenWindow', screenNo, black, ...
    [0 0 1920 1080]); %, 32, 2,[], numMultiSamples, []);
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
% fastKey = KbName('RightArrow');
% slowKey = KbName('LeftArrow');
slowKey = KbName('RightArrow');
fastKey = KbName('LeftArrow');
equalKey = KbName('DownArrow');

%TIME -----------------------------------------------------------
durFix = 0.5;
delay = flipInterval/2;


% TEST BUFFER ---------------------------------------------------
disp('Test begin')
DrawFormattedText(window,'Trials\n\nPress Any Key to Start','center','center',grey);
Screen('Flip', window);
KbStrokeWait(deviceIndex);

loopFlag = 1;

%Setup parameters
redLum = grbVal;   %RGB value
greenLum = grbVal;
m = 0.5;
fs =0.01;  %spatial frenquency 
ft = 15; %temporal frequency
red = zeros(1,600);
green = zeros(1,600);
R = zeros(100,600);
G = zeros(100,600);
RGB = zeros(100,600,3); 
vbl = Screen('Flip',window);
dstRects = [xCenter-300 yCenter-50 xCenter+300 yCenter+50];
RGBTexture = Screen('MakeTexture', window, zeros(100,600));

%----------------------------------------------------------------
%                        EXPERIMENT LOOP
%----------------------------------------------------------------
while loopFlag
    [keyIsDown, secs, keyCode] = KbCheck;
    
    for x=1:600
        red(x)=0.5*redLum*((1+m*sin(2*pi*fs*x)*sin(2*pi*ft*vbl))+(1+cos(2*pi*fs*x)*cos(2*pi*ft*vbl)));
        green(x)=0.5*greenLum*((1+m*sin(2*pi*fs*x)*sin(2*pi*ft*vbl))+(1-cos(2*pi*fs*x)*cos(2*pi*ft*vbl)));
    end
    R = repmat(red,100,1);
    G = repmat(green,100,1);
%     R = R/255;  %rescaling
%     G = G/255;
    RGB(:,:,1)=R;
    RGB(:,:,2)=G;
%     RGB(1:10,1:10,1)
%     RGB(1:10,1:10,2)
%     RGB(1:10,1:10,3)
    Screen('Close',RGBTexture);
    RGBTexture = Screen('MakeTexture', window, RGB);
    Screen('DrawTextures', window, RGBTexture, [], dstRects);
%     Screen('DrawTextures', window, RGBTexture);
    vbl = Screen('Flip',window,vbl+delay);
    
    if keyCode(escapeKey)||keyCode(slowKey)||keyCode(fastKey)||keyCode(equalKey)
        if keyCode(escapeKey)
            ShowCursor;
            disp('Experiment aborted');
            sca;
            return;
        elseif keyCode(slowKey)
            redLum = redLum - 1;
        elseif keyCode(fastKey)
            redLum = redLum + 1 ;
        else
            loopFlag = 0;
            disp('E');
            disp(redLum);
            disp(greenLum);
            break;
        end
    end
   %Keyboard check ends
end


% CLOSE --------------------------------------------------------

DrawFormattedText(window,...
    'Press Any Key To Exit',...
    'center', 'center', grey);
Screen('Flip', window);
KbStrokeWait(deviceIndex);
sca;
disp('Test ends.')

%===============================END=============================

end
