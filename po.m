function po

%% Preparing the variables.

% Current test variables:
subj = 1001;
domEye = 1; % 0=right, 1=left
condFileName = 'po-cond01';

% Some default variables:
nRevs = 10; % number of reversals for the staircases
gabPhase = 0; % pase of underlying sine grating in degrees
sc = 5; % spatial constant of the exponential "hull"
radius = 70; % for Gabor arrangement

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
% white = WhiteIndex(screenid);
% black = BlackIndex(screenid);
% gray = round((white+black)/2);
% if gray == white
%     gray = white / 2;
% end
% inc = white - gray; % increment
backgroundCol = 0;

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
disp.centX(1:2) = [disp.resX/2-disp.distX disp.resX/2+disp.distX]; % 230(L) & 730(R)
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
        'stopcriterion','reversals','stoprule',nRevs);  %#ok<AGROW>
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

%% Preparing the screen.

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


%% Preparing the Mondrian files.
mondFiles    = dir('mondrians/');
mondFileIds  = find(~[mondFiles.isdir]);
mondFileIds2 = mondFileIds;
for k=1:length(mondFileIds) % some checks, I guess
    if isempty(strfind(mondFiles(mondFileIds(k)).name, '.jpg')),
        mondFileIds2(k) = [];
    end
end
% This is the cell variable that contains all of the images (no need to
% reopen the files on the fly):
mondImg = cell(1,length(mondFileIds2));
for k=1:length(mondFileIds2)
    mondImg{k} = imread(['mondrians/' mondFiles(mondFileIds2(k)).name]);
end

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
        
    % Duration and mask settings:
    stimT = c2n(condTable,'stimT',curStairc);
    postStimBlankT = c2n(condTable,'postStimBlankT',curStairc);
    odtTilt = c2n(condTable,'odtTilt',curStairc);
    odtT = c2n(condTable,'odtT',curStairc);
    maskRR = c2n(condTable,'maskRR',curStairc);
    
    % Gabor settings:
    gabNum = c2n(condTable,'gabNum',curStairc); % number of Gabors
    gabSize = c2n(condTable,'gabSize',curStairc);
    gabSize = round(cm2px(gabSize, sdims));
    gabSf = c2n(condTable,'gabSf',curStairc);
    
    % Randomly assigned vars:
    gabOri(1:gabNum) = repmat((randi(2)-1)*90, 1, gabNum); % 0=hori, 1=vert
    singlLoc = randi(gabNum);
    display(sprintf('singleton ID: %.2f', singlLoc));
    primCol = randi(2); % 1=red, 2=green
    
    maxContr = staircs(curStairc).xCurrent; % current max contrast
    display(sprintf('singleton brightness: %.2f', maxContr));


    %% Rendering the gratings:
    numFrames = round(stimT*frameRate/1000);

    % scale this patch up and down to draw individual patches of the different
    % wanted sizes:
    si = gabSize; % 200;
    % Size of support in pixels, derived from si:
    tw = 2*si+1;
    th = 2*si+1;

    % Build a procedural gabor texture for a gabor with a support of tw x th
    % pixels and the 'nonsymetric' flag set to 1 == Gabor shall allow runtime
    % change of aspect-ratio:
    gab = CreateProceduralGabor(wPtr, tw, th, 1);
    % The rectangles in which the Gabors will be drawn
    texrect = Screen('Rect', gab);
    
    % Locations:
    fiStep = 360/gabNum; % fi steps: 360/12 = 30 degrees
    
    % Cycling through Gabors to adjust the settings:
    dstRects = zeros(gabNum,4);
    gabCol = zeros(gabNum,3);
    for curGabI = 1:gabNum
        % The location for each Gabor (to the non-domEye):
        dstRects(curGabI,:) = CenterRectOnPoint(texrect, ...
            disp.centX(1+domEye) + radius*cosd(fiStep/2+curGabI*fiStep), ...
            disp.centY + radius*sind(fiStep/2+curGabI*fiStep));
        % Colours:
        if singlLoc==curGabI
            gabCol(curGabI,primCol) = 255;
        else
            gabCol(curGabI,1:2) = 255;
            gabCol(curGabI,primCol) = 0;
        end
    end
    
    for i=1:numFrames
        %% Animation loop for the priming stage.
        % The contrast varies through the frames:
        curContr = ( exp(-(-1+(2*i/numFrames))^2*4) )*maxContr*100;
        % Cycling through the Gabors:
        for curGabI = 1:gabNum
            Screen('DrawTexture', wPtr, gab, [], dstRects(curGabI,:), ...
                gabOri(curGabI),[], [], gabCol(curGabI,:), [], ...
                kPsychDontDoRotation, [gabPhase, gabSf, sc, curContr, 1, 0, 0, 0]);
        end
        % Fixation boxes:
        drawFixationBox(wPtr, disp.centX(1), disp.centY, disp.boxSize, ...
            disp.boxColour); % left fixation box
        drawFixationBox(wPtr, disp.centX(2), disp.centY, disp.boxSize, ...
            disp.boxColour); % right fixation box
        % Mask:
        drawMondrians(mondImg{rem(i,10)}, wPtr, disp.centX(2-domEye), disp.centY, disp.boxSize);
        Screen('Flip', wPtr);
        % Monitoring for keypresses.
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