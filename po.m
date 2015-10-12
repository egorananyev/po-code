function po

%% Preparing the variables.

% Current test variables:
subj = 1001;
domEye = 0; % 0=right, 1=left
condFileName = 'po-cond01';
nRevs = 10; % number of reversals for the staircases
 
% Keyboard:
KbName('UnifyKeyNames');
quitkey = 'c';
space = 'space';
% targetUsageName = 'Keyboard'; % change accordingly
% targetProduct = 'Dell USB Keyboard'; % change accordingly
% targetProduct = 'Apple Keyboard'; % temp
% dev = PsychHID('Devices');
% deviceIndex = find(strcmpi(targetUsageName, {dev.usageName}) & ...
%     strcmpi(targetProduct, {dev.product}));
deviceIndex = -3;
KbQueueCreate(deviceIndex);
KbQueueStart(deviceIndex);

% Output file name:
dateNtime = datestr(now,'yyyy-mm-dd_HHMMSS');
sessionName = strcat('po_', condFileName', '_s', mat2str(subj), ...
    '_d', mat2str(domEye)', '_', dateNtime);
outFileName = strcat('../po_data/', sessionName, '.csv');

% Instructions text:
textInstr = 'Press any button to continue';
textNextTrial = 'Press spacebar to continue';

%% Preparing PsychToolBox and screen.

% Prepping PsychToolBox:
Screen('Preference', 'SkipSyncTests', 1); % a necessary evil - only for Yosemite

AssertOpenGL; % for 3D rendering 
screenid = max(Screen('Screens'));
% InitializeMatlabOpenGL(1); % for 3D rendering
% PsychImaging('PrepareConfiguration'); 

% Some verification for colour names (for the gratings):
white = WhiteIndex(screenid);
black = BlackIndex(screenid);
gray = round((white+black)/2);
if gray == white
    gray = white / 2;
end
inc = white - gray; % increment
backgroundCol = black;

try
[wPtr, rect] = Screen('OpenWindow', screenid, backgroundCol);
% wPtr = 10 % the number that designates the created window;
% rect = [0 0 1920 1080] % RectLeft=1, RectTop=2, RectRight=3, RectBottom=4
% Screen('BlendFunction', wPtr, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% Screen('CloseAll'); %temp

%% Resolution and display locations.

% screen resolution, in pixels:
disp.resX = rect(3); % 1920
disp.resY = rect(4); % 1080

% display center locations, in pixels:
disp.boxColour = [255 255 255]; % white
disp.boxSize = 90;
disp.distX = 170; % 150; % display center distance from the vertical midline (left/right)
disp.distY = -100; % display center distance from the horizontal midline (up)
disp.centX(1:2) = [disp.resX/2-disp.distX disp.resX/2+disp.distX]; % 230(R) & 730(L)
disp.centY = disp.resY/2 + disp.distY; % 540-200=340

%Screen('CloseAll'); %temp

%% Experimental conditions.

% Reading the conditions file with the settings
[~,~,condTable] = xlsread(strcat(condFileName, '.xlsx'));
numofConds = size(condTable,1)-1; % number of conditions

% Also reading screen dimensions for px2cm and cm2px conversions:
[~,~,sdims] = xlsread('screenDims.xlsx');

% Setting up the staircases.
for i=1:numofConds
    staircs(i) = PAL_AMUD_setupUD('up',1,'down',1,...
        'stepsizeup',c2n(condTable,'vcUp',i),...
        'stepsizedown',c2n(condTable,'vcDown',i),...
        'startvalue',c2n(condTable,'vcSt',i),...
        'stopcriterion','reversals','stoprule',nRevs); %#ok<AGROW>
end

%% Presenting the instructins window.
Screen('TextFont', wPtr, 'Cambria');
Screen('TextSize', wPtr, 32);
% Defining the edges of the left and right box:
boxL = [disp.centX(1)-disp.boxSize disp.centY-disp.boxSize...
    disp.centX(1)+disp.boxSize disp.centY+disp.boxSize];
boxR = [disp.centX(2)-disp.boxSize disp.centY-disp.boxSize...
    disp.centX(2)+disp.boxSize disp.centY+disp.boxSize];
% Drawing the fixation box on the left:
drawFixationBox(wPtr, disp.centX(1), disp.centY, disp.boxSize, disp.boxColour);
% Drawing the text box centered in the left box:
DrawFormattedText(wPtr, textInstr, 'center', 'center', [255 255 50], 15, ...
    [], [], [], [], boxL)
% Drawing the fixation and text boxes on the right:
drawFixationBox(wPtr, disp.centX(2), disp.centY, disp.boxSize, disp.boxColour);
DrawFormattedText(wPtr, textInstr, 'center', 'center', [255 255 50], 15, ...
    [], [], [], [], boxR)
Screen(wPtr, 'Flip');

%% Monitoring for the key presses during the instruction:
while 1,
%     [keyIsDownPrev, secsPrev, keyCodePrev] = KbQueueCheck(deviceIndex);
    [keyIsDown, ~, keyCode] = KbCheck(deviceIndex); % not sure if this is the
        % optimal command in terms of waiting.
    if keyIsDown,
        % If the quit key is pressed, quit:
        if keyCode(KbName(quitkey)),
            Screen('CloseAll');
            ShowCursor;
            return;
        else
        % If any other key is pressed, proceed after .5 seconds:
            break;
        end
    end
end
WaitSecs(0.5);

%% Preparing the gratings and screen.

% Run the movie animation for a fixed period.
% curTrial = 5; %TEMP
frameRate=Screen('FrameRate',screenid);
% If MacOSX does not know the frame rate the 'FrameRate' will return 0.
% That usually means we run on a flat panel with 60 Hz fixed refresh
% rate:
if frameRate == 0
    frameRate = 60;
end

% Use realtime priority for better timing precision:
priorityLevel=MaxPriority(wPtr);
Priority(priorityLevel);

%% Going through the staircases (running the trials).

% Making sure there are still staircases to complete:
numofStaircs = numofConds;
while numofStaircs>0 % while there are still staircs to complete

% Random-shuffling the active (remaining) staircs:
activeStaircIndcs = find([staircs.reversal]<10);
trialCond = activeStaircIndcs(randperm(numofStaircs));

% Some info on this iteration of staircs:
display('=====#=====');
display(sprintf('number of remaining staircases: %i', numofStaircs));

% Screen('CloseAll'); %temp

for curTrial=1:numofStaircs
    %% Drawing the gratings for this trial.
    display('=====');
    curStairc = trialCond(curTrial);
    display(sprintf('current trial#: %3i', curTrial));
    display(sprintf('current staircase: %i', curStairc));
    display(sprintf('staircase start value: %.2f', c2n(condTable,'vcSt',curStairc)));
    
    % Singleton config:
    singlType = c2n(condTable,'singlType',curStairc);
    singlCont = c2n(condTable,'singlCont',curStairc);
    display(['singleton type: ', singlType]);
    display(sprintf('singleton contrast: %.2f', singlCont));
    % Randomly assigned vars:
    gabOri = randi(2)-1; % 0=horizontal, 1=vertical
    gabNum = c2n(condTable,'gabNum',curStairc);
    singlLoc = randi(gabNum);
        
    % Duration, mask, and grating settings:
    stimT = c2n(condTable,'stimT',curStairc);
    postStimBlankT = c2n(condTable,'postStimBlankT',curStairc);
    gabSize = c2n(condTable,'gabSize',curStairc);
    gabSf = c2n(condTable,'gabSf',curStairc);
    maskRR = c2n(condTable,'maskRR',curStairc);
    odtTilt = c2n(condTable,'odtTilt',curStairc);
    odtT = c2n(condTable,'odtT',curStairc);
    
    % Converting Gabor settings into pixels:
    gabSize = round(cm2px(gabSize, sdims));
%     gabSf = round(cm2px(gabSf, sdims))*2*pi; % 0.0142 cy/px = .5 cy/cm

    %% Rendering the gratings:
    numFrames = round(stimT*frameRate/1000);
    maxContr = staircs(curStairc).xCurrent; % current max contrast
    for i = 1:numFrames
        curContr = ( exp(-(-1+(2*i/numFrames))^2*4) )*maxContr;
        m = renderGrating(gabSize, gabOri*90, .5, 0);
        %         m = .75*m;
%         m(m<0)=0;
        m = m - m*.5;
        display(sprintf('m min=%.2f; max=%.2f',min(min(m)),max(max(m))));
        ftwindow = renderWindow(m, 1);
%         im = gray + inc*m;
%         im = white*m; % *curContr;
%         display(sprintf('im min=%.2f; max=%.2f',min(min(im)),max(max(im))));
        im = curContr*gray + curContr*inc*m;
        im(im<0)=0;
        display(sprintf('im min=%.2f; max=%.2f',min(min(im)),max(max(im))));
%         gab(i) = Screen('MakeTexture', wPtr, ftwindow.*im ); %#ok<AGROW>
        
    end

    %% Animation loop.
    % scale this patch up and down to draw individual patches of the different
    % wanted sizes:
    si = 200;

    % Size of support in pixels, derived from si:
    tw = 2*si+1;
    th = 2*si+1;

    % Initial parameters of gabors:

    % Phase of underlying sine grating in degrees:
    phase = 0;
    % Spatial constant of the exponential "hull"
    sc = 20.0;
    % Frequency of sine grating:
    freq = .1;
    % Contrast of grating:
    contrast = 20.0;
    % Aspect ratio width vs. height:
    aspectratio = 1.0;
    % Build a procedural gabor texture for a gabor with a support of tw x th
    % pixels and the 'nonsymetric' flag set to 1 == Gabor shall allow runtime
    % change of aspect-ratio:
    rotAngles = 90;
    gab = CreateProceduralGabor(wPtr, tw, th, 1);
    mypars = [phase+180, freq, sc, contrast, aspectratio, 0, 0, 0]';
    Screen('DrawTexture', wPtr, gab, [], [], [], [], [], [], [], ...
        kPsychDontDoRotation, [phase, freq, sc, contrast, aspectratio, 0, 0, 0]);
    texrect = Screen('Rect', gab);
%     dstRects = zeros(4);
    dstRects(1:4) = CenterRectOnPoint(texrect, disp.centX(1), disp.centY)';
    
    for i=1:numFrames %frameSet
        % Draw the left image:
%         boxMulti = 1; % 1.8;
%         Screen('DrawTexture', wPtr, gab(i), [], ...
%             [disp.centX(1)-disp.boxSize*boxMulti disp.centY-disp.boxSize*boxMulti ...
%             disp.centX(1)+disp.boxSize*boxMulti disp.centY+disp.boxSize*boxMulti]);
        Screen('DrawTextures', wPtr, gab, [], dstRects, rotAngles,...
            [], [], [], [], kPsychDontDoRotation, mypars);
        drawFixationBox(wPtr, disp.centX(1), disp.centY, disp.boxSize, ...
            disp.boxColour);
        Screen('Flip', wPtr);

        %% Monitoring for keypresses.
        [keyIsDown, ~, keyCode] = KbCheck(deviceIndex);
        if keyIsDown,
            if keyCode(KbName(quitkey)),
                Screen('CloseAll');
                ShowCursor;
                return;
            end
        end
    end
    
    %% Continue screen.
    % Drawing the fixation box on the left:
    drawFixationBox(wPtr, disp.centX(1), disp.centY, disp.boxSize, disp.boxColour);
    % Drawing the text box centered in the left box:
    DrawFormattedText(wPtr, textNextTrial, 'center', 'center', [255 255 50], 15, ...
        [], [], [], [], boxL)
    % Drawing the fixation and text boxes on the right:
    drawFixationBox(wPtr, disp.centX(2), disp.centY, disp.boxSize, disp.boxColour);
    DrawFormattedText(wPtr, textNextTrial, 'center', 'center', [255 255 50], 15, ...
        [], [], [], [], boxR)
    Screen(wPtr, 'Flip');
    % Monitoring for a keypress:
    while 1,
        [keyIsDown, ~, keyCode] = KbCheck(deviceIndex); % not sure if this is the
            % ... optimal command in terms of waiting.
        if keyIsDown,
            % If the 'quit' key is pressed, quit:
            if keyCode(KbName(quitkey)),
                Screen('CloseAll');
                ShowCursor;
                return;
            elseif keyCode(KbName(space)),
            % If the 'space' key is pressed, proceed to the next trial:
                break;
            end
        end
    end

end

% Updating the number of remaining staircs for the while loop:
numofStaircs = sum([staircs.reversal]<nRevs);

end

Screen('CloseAll');

catch %#ok<CTCH>
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    Priority(0);
    psychrethrow(psychlasterror);
end