function prelimThresh(subj, domEye)
%% A function to evaluate preliminary visibility threshold for a pop-out/no-pop-out displays.
%subj = 1001;
%domEye = 0; % 0=right, 1=left

%% Preparing the variables.

% Some default variables:
condFileName = 'po-prelimThresh';
gabPhase = 0; % pase of underlying sine grating in degrees
sc = 5; % spatial constant of the exponential "hull"
radius = 70; % for Gabor arrangement
feedbackFC = 0; % the number of frames for feedback
greenRgb = 165; % find this value through equiluminance testing
redRgb = 255;

% Keyboard:
KbName('UnifyKeyNames');
quitkey = 'c';
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
sessionName = strcat(condFileName, '_s', mat2str(subj), ...
    '_d', mat2str(domEye)', '_', dateNtime);
outFileName = strcat('../po-thresh/', sessionName);

% Instructions text:
% textTilt = 'Please\nindicate\nthe tilt:\n<- (left)\n-> (right)'; % 4 lines
% textTiltFdbLeft = ' \n \n \n<- (left)\n ';
% textTiltFdbRight = ' \n \n \n \n-> (right)';
textTilt = ' ';
textTiltFdbLeft = ' ';
textTiltFdbRight = ' ';
textObjVis = 'The "circles"\nwere:\npresent (up)\nabsent (down)';
textObjVisFdbPres = ' \n \npresent (up)\n ';
textObjVisFdbAbs = ' \n \n \nabsent (down)';
textSubjVis = '1. no experience\n2. brief glimpse\n3. almost clear\n4. clear experience';
textSubjVisFdb1 = '1. no experience\n \n \n ';
textSubjVisFdb2 = ' \n2. brief glimpse\n \n ';
textSubjVisFdb3 = ' \n \n3. almost clear\n ';
textSubjVisFdb4 = ' \n \n \n4. clear experience';
textInstr = 'Press\nany button\nto start';
textNextTrial = 'Press spacebar\nto continue';

%% Preparing PsychToolBox and screen.

% Prepping PsychToolBox:
% Screen('Preference', 'SkipSyncTests', 1); % a necessary evil - only for Yosemite
% Screen('Preference', 'SuppressAllWarnings', 1);

AssertOpenGL; % for 3D rendering 
screenid = max(Screen('Screens'));
% screenid = min(Screen('Screens')); %TEMP
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
alphas = -1:.1:2; % assuming the most likely threshold is .5, it gives 1.5 on both sides
for i=1:numofConds
    % Priors and alphas for the best PEST staircases:
    priorMean = c2n(condTable,'priorMean',i);
    priorSD = c2n(condTable,'priorSD',i);
    prior = PAL_pdfNormal(alphas, priorMean, priorSD);
    staircs(i) = PAL_AMRF_setupRF('xMin',0,'xMax',2,...
       'priorAlphaRange', alphas, 'prior', prior,...
       'stopCriterion','reversals','stopRule',c2n(condTable,'nRevs',i)); %#ok<*AGROW>
end
% plot(alphas, prior);

%% Presenting the instructins window.
Screen('TextFont', wPtr, 'Cambria');
Screen('TextSize', wPtr, 26);
Screen('TextColor', wPtr, [255 255 50]);
% Defining the edges of the left and right box:
mult = 11; % for some weird reason, on some computers, the bounding boxes
% don't work correctly - you need much larger bounding boxes than
% necessitated by the font size to even display the text (2015-10-26)
boxL = [disp.centX(1)-disp.boxSize*mult disp.centY-disp.boxSize*mult ...
    disp.centX(1)+disp.boxSize*mult disp.centY+disp.boxSize*mult];
boxR = [disp.centX(2)-disp.boxSize*mult disp.centY-disp.boxSize*mult ...
    disp.centX(2)+disp.boxSize*mult disp.centY+disp.boxSize*mult];
% Drawing the fixation box on the left:
drawFixationBox(wPtr, disp.centX(1), disp.centY, disp.boxSize, disp.boxColour);
% Drawing the text box centered in the left box:
DrawFormattedText(wPtr, textInstr, 'center', 'center', [], [], ...
    [], [], [], [], boxL);
% Drawing the fixation and text boxes on the right:
drawFixationBox(wPtr, disp.centX(2), disp.centY, disp.boxSize, disp.boxColour);
DrawFormattedText(wPtr, textInstr, 'center', 'center', [], [], ...
    [], [], [], [], boxR);
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

% Logging the output of the Matlab window:
diary(strcat(outFileName, '_log.txt'));

% Run the movie animation for a fixed period.
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
numofStops = 0; % the starting number of staircase stop = 0
c = 0; % initiating trial counter:
fc = 0; % frame counter for feedback
bc = 0; % block counter
%while numofStops<numofConds % while there are still staircs to complete
while numofStaircs > 0

bc = bc + 1; % block counter

% Random-shuffling the active (remaining) staircs:
% activeStaircIndcs = find([staircs.reversal]<nRevs);
activeStaircIndcs = find([staircs.stop]==0);
trialCond = activeStaircIndcs(randperm(numofStaircs));

% Some info on this iteration of staircs:
display('=====#=====');
display(sprintf('number of remaining staircases: %i', numofStaircs));

% Screen('CloseAll'); %temp

for blockTrial=1:numofStaircs
    %% Drawing the gratings for this trial.
    c = c + 1; % trial counter
    
    % All of the trial-specific data are recorded in structure d.
    display('=====');
    d.block(c) = bc;
    display(sprintf('current block count: %i', bc));
    d.trial(c) = c;
    display(sprintf('current trial count: %i', c));
    d.curStairc(c) = trialCond(blockTrial);
    d.blockTrial(c) = blockTrial;
    display(sprintf('current block trial: %i', blockTrial));
    display(sprintf('current staircase: %i', d.curStairc(c)));
    
    % Singleton config:
    d.singlCont(c) = c2n(condTable,'singlCont',d.curStairc(c)); % singlton contrast
    display(sprintf('singleton contrast: %.2f', d.singlCont(c)));
        
    % Duration and mask settings:
    d.jitTmax(c) = c2n(condTable,'jitTmax',d.curStairc(c));
    d.stimT(c) = c2n(condTable,'stimT',d.curStairc(c));
    d.postStimBlankT(c) = c2n(condTable,'postStimBlankT',d.curStairc(c));
    d.maskRR(c) = c2n(condTable,'maskRR',d.curStairc(c));
    d.maskOnOff(c) = c2n(condTable,'maskOnOff',d.curStairc(c));
    
    % Gabor settings:
    d.gabNum(c) = c2n(condTable,'gabNum',d.curStairc(c)); % number of Gabors
    d.gabSize(c) = c2n(condTable,'gabSize',d.curStairc(c));
    d.gabSize(c) = round(cm2px(d.gabSize(c), sdims));
    d.gabSf(c) = c2n(condTable,'gabSf',d.curStairc(c));
    
    %% Other variables.
    
    % Randomly assigned vars:
    d.gabOri(c) = (randi(2)-1)*90; % 0=hori, 1=vert
    d.singlLoc(c) = randi(d.gabNum(c));
    display(sprintf('singleton ID: %i', d.singlLoc(c)));
    d.primCol(c) = randi(2); % 1=red, 2=green
    if d.primCol(c)==1, display('primary colour: red');
    else (display('primary colour: green')); end
    
    % Timing vars:
    d.jitFC(c) = randi(round(d.jitTmax(c)*frameRate/1000)); % random jitter
    d.stimFC(c) = round(d.stimT(c)*frameRate/1000);
    d.postStimBlankFC(c) = round(d.postStimBlankT(c)*frameRate/1000);
    d.trialFC(c) = d.jitFC(c)+d.stimFC(c)+d.postStimBlankFC(c);
    
    % Resetting the response variables:
    respMadeObjVis = false;
    respMadeSubjVis = false;
    respMadeCont = false;
    respFeedbackObjVis = false;
    respFeedbackSubjVis = false;
    
    %% Rendering the gratings:

    d.staircVal(c) = staircs(d.curStairc(c)).xCurrent; % current staircase value
    display(sprintf('stairscase value: %.3f', d.staircVal(c)));
    % based on the above staircase value, setting the mask contrast:
    if d.staircVal(c) <= 1
        d.maskContr(c) = 1;
        if d.staircVal(c) < 0
            d.maxContr(c) = 0;
        else
            d.maxContr(c) = d.staircVal(c);
        end
    else
        d.maskContr(c) = 2-d.staircVal(c);
        d.maxContr(c) = 1;
    end
    % If this is a blank trial:
    display(sprintf('stimulus visual contrast: %.3f', d.maxContr(c)));
    display(sprintf('mask visual contrast: %.2f', d.maskContr(c)));
    
    % scale this patch up and down to draw individual patches of the different
    % wanted sizes:
    d.si(c) = d.gabSize(c); % 200;
    % Size of support in pixels, derived from d.si(c):
    tw = 2*d.si(c)+1;
    th = 2*d.si(c)+1;

    % Build a procedural gabor texture for a gabor with a support of tw ~ th
    % pixels and the 'nonsymetric' flag set to 1 == Gabor shall allow runtime
    % change of aspect-ratio:
    gab = CreateProceduralGabor(wPtr, tw, th, 1);
    % The rectangles in which the Gabors will be drawn
    texrect = Screen('Rect', gab);
    
    % Locations:
    d.fiStep(c) = 360/d.gabNum(c); % fi steps: 360/12 = 30 degrees
    
    % Cycling through Gabors to adjust the settings:
    dstRects = zeros(d.gabNum(c),4);
    gabCol = zeros(d.gabNum(c),3);
    % RGB value for primary and pop-out colours:
    if d.primCol(c)==1, primRgb=redRgb; singlRgb=greenRgb; singlCol=2;
    else primRgb=greenRgb; singlRgb=redRgb; singlCol=1; end
    for curGabI = 1:d.gabNum(c)
        % The location for each Gabor (to the non-domEye):
        dstRects(curGabI,:) = CenterRectOnPoint(texrect, ...
            disp.centX(1+domEye) + radius*cosd(d.fiStep(c)/2+curGabI*d.fiStep(c)), ...
            disp.centY + radius*sind(d.fiStep(c)/2+curGabI*d.fiStep(c)));
        % Colours:
        if d.singlLoc(c)==curGabI % for the singleton/pop-out
            gabCol(curGabI,singlCol) = singlRgb*d.singlCont(c);
            gabCol(curGabI,d.primCol(c)) = primRgb*(1-d.singlCont(c));
        else % for Gabors other than the outlier
            gabCol(curGabI,d.primCol(c)) = primRgb;
        end
    end
    
    curFrame = 0; % (re)setting the current frame count
    mondIndx = 1; % (re)assigning the mondrean index
    while ~respMadeCont
        %% Animation loop.
        curFrame = curFrame + 1;
        % 1. The priming stage.
        if curFrame>d.jitFC(c) && curFrame<=d.stimFC(c)+d.jitFC(c)
            % The contrast varies through the frames:
            curContr = ( exp(-(-1+(2*curFrame/d.stimFC(c)))^2*4) )*d.maxContr(c)*100;
            % Cycling through the Gabors:
            for curGabI = 1:d.gabNum(c)
                Screen('DrawTexture', wPtr, gab, [], dstRects(curGabI,:), ...
                    d.gabOri(c),[], [], gabCol(curGabI,:), [], ...
                    kPsychDontDoRotation, [gabPhase, d.gabSf(c), sc, curContr, 1, 0, 0, 0]);
            end
        end
        % 2. Post-stimulus blank.
        % Do nothing. The mask is still shown.
        %% 3. Response: subjective visibility scale.
        if curFrame>d.trialFC(c) && ~respMadeSubjVis
            % Displaying the "subjective visibility" text:
            Screen('TextSize', wPtr, 21);
            DrawFormattedText(wPtr, textSubjVis, 'center', 'center', ...
                [255 255 50], [], [], [], [], [], boxL);
            DrawFormattedText(wPtr, textSubjVis, 'center', 'center', ...
                [255 255 50], [], [], [], [], [], boxR);
            % Monitoring for up/down objective visibility responses.
            [keyIsDown, ~, keyCode] = KbCheck(deviceIndex);
            if keyIsDown,
                if keyCode(KbName('1!'))
                    d.respSubjVis(c) = 1;
                    respMadeSubjVis = true;
                    display('Response made: subjective visibility score 1');
                    % Only update the response here if this is a non-blank staircase:
                    display('Non-blank trial: Recording "unseen" response');
                    staircs(d.curStairc(c)) = PAL_AMRF_updateRF(staircs(d.curStairc(c)), ...
                        d.staircVal(c), 0);
                elseif keyCode(KbName('2@'))
                    d.respSubjVis(c) = 2;
                    respMadeSubjVis = true;
                    display('Response made: subjective visibility score 2');
                    display('Non-blank trial: Recording "seen" response');
                    staircs(d.curStairc(c)) = PAL_AMRF_updateRF(staircs(d.curStairc(c)), ...
                        d.staircVal(c), 1);
                elseif keyCode(KbName('3#'))
                    d.respSubjVis(c) = 3;
                    respMadeSubjVis = true;
                    display('Response made: subjective visibility score 3');
                    display('Non-blank trial: Recording "seen" response');
                    staircs(d.curStairc(c)) = PAL_AMRF_updateRF(staircs(d.curStairc(c)), ...
                        d.staircVal(c), 1);
                elseif keyCode(KbName('4$'))
                    d.respSubjVis(c) = 4;
                    respMadeSubjVis = true;
                    display('Response made: subjective visibility score 4');
                    display('Non-blank trial: Recording "seen" response');
                    staircs(d.curStairc(c)) = PAL_AMRF_updateRF(staircs(d.curStairc(c)), ...
                        d.staircVal(c), 1);
                end
            end
        end
        % Response feedback.
        if respMadeSubjVis && ~respFeedbackSubjVis
            fc = fc + 1;
            if d.respSubjVis(c) == 1
                DrawFormattedText(wPtr, textSubjVisFdb1, 'center', 'center', ...
                    [255 255 50], 20, [], [], [], [], boxL);
                DrawFormattedText(wPtr, textSubjVisFdb1, 'center', 'center', ...
                    [255 255 50], 20, [], [], [], [], boxR);
            elseif d.respSubjVis(c) == 2
                DrawFormattedText(wPtr, textSubjVisFdb2, 'center', 'center', ...
                    [255 255 50], 20, [], [], [], [], boxL);
                DrawFormattedText(wPtr, textSubjVisFdb2, 'center', 'center', ...
                    [255 255 50], 20, [], [], [], [], boxR);
            elseif d.respSubjVis(c) == 3
                DrawFormattedText(wPtr, textSubjVisFdb3, 'center', 'center', ...
                    [255 255 50], 20, [], [], [], [], boxL);
                DrawFormattedText(wPtr, textSubjVisFdb3, 'center', 'center', ...
                    [255 255 50], 20, [], [], [], [], boxR);
            else
                DrawFormattedText(wPtr, textSubjVisFdb4, 'center', 'center', ...
                    [255 255 50], 20, [], [], [], [], boxL);
                DrawFormattedText(wPtr, textSubjVisFdb4, 'center', 'center', ...
                    [255 255 50], 20, [], [], [], [], boxR);
            end
            if fc >= feedbackFC
                respFeedbackSubjVis = true;
                fc=0;
                % Writing into d. the current reversals and mean/sd estimates:
                curRevs = staircs(d.curStairc(c)).reversal;
                d.revs(c) = curRevs(end);
                if ~isempty(staircs(d.curStairc(c)).meanUniformPrior)
                    d.meanUniformPrior(c) = staircs(d.curStairc(c)).meanUniformPrior;
                    d.sdUniformPrior(c) = staircs(d.curStairc(c)).sdUniformPrior;
                    d.mean(c) = staircs(d.curStairc(c)).mean;
                    d.sd(c) = staircs(d.curStairc(c)).sd;
                else
                    d.meanUniformPrior(c) = 0;
                    d.sdUniformPrior(c) = 0;
                    d.mean(c) = 0;
                    d.sd(c) = 0;
                end
            end
        end
        %% 4. Response: continue.
        if respFeedbackSubjVis && ~respMadeCont
            % Displaying the "continue" text:
            Screen('TextSize', wPtr, 32);
            DrawFormattedText(wPtr, textNextTrial, 'center', 'center', ...
                [255 255 50], 14, [], [], [], [], boxL);
            DrawFormattedText(wPtr, textNextTrial, 'center', 'center', ...
                [255 255 50], 14, [], [], [], [], boxR);
            % Monitoring for a "space" press:
            [keyIsDown, ~, keyCode] = KbCheck(deviceIndex);
            if keyIsDown,
                if keyCode(KbName('space'))
                    Screen('Close'); % for optimization
                    respMadeCont = true;
                    display('Response made: continuing to next trial');
                    % Data recording.
                    % The data in the .mat format:
                    save([outFileName '.mat'], 'd'); % the data structure
                    save([outFileName '_staircs.mat'], 'staircs'); % the staircases
                    % The data in the csv format:
                    dColNames = fieldnames(d)';
                    dMat = cell2mat(struct2cell(d))';
                    csvwrite([outFileName '.csv'], dMat);
                    cell2csv([outFileName '_colNames.csv'], dColNames);
                    % The output of the Matlab output window:
                    diary(strcat(outFileName, '_log.txt'));
                    display('Data recorded');
                end
            end
        end
        %% Tail
        % Mask:
        if curFrame<=(d.jitFC(c)+d.stimFC(c)+d.postStimBlankFC(c)) & d.maskOnOff(c)
            % Alternating frequency is set by maskRR (refresh rate):
            if rem(curFrame,d.maskRR(c))==0
                mondIndx = mondIndx + 1;
                if rem(curFrame,d.maskRR(c)*length(mondFileIds2))==0
                    mondIndx = 1;
                end
            end
            % Drawing the mask based on the Mondrean index:
            drawMondrians(mondImg{mondIndx}, wPtr, ...
                disp.centX(2-domEye), disp.centY, disp.boxSize, d.maskContr(c));
        end
        % Fixation boxes:
        drawFixationBox(wPtr, disp.centX(1), disp.centY, disp.boxSize, ...
            disp.boxColour); % left fixation box
        drawFixationBox(wPtr, disp.centX(2), disp.centY, disp.boxSize, ...
            disp.boxColour); % right fixation box
        Screen('DrawingFinished', wPtr); % an optimization bit
        Screen('Flip', wPtr);
        % Monitoring for the "quit key" press.
        [keyIsDown, ~, keyCode] = KbCheck(deviceIndex);
        if keyIsDown,
            if keyCode(KbName(quitkey)),
                display('Quit key pressed');
                Screen('CloseAll');
                ShowCursor;
                return;
            end
        end
    end
end

% Updating the number of remaining staircs for the while loop:
numofStops = sum([staircs.stop]);
numofStaircs = numofConds - numofStops;

% Displaying the number of reversals:
revs = zeros(1, size(staircs,2));
for i = 1:size(staircs,2)
    revs(i) = max(staircs(i).reversal);
end
display(revs);

end

Screen('CloseAll');

catch %#ok<CTCH>
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    Priority(0);
    psychrethrow(psychlasterror);
end
