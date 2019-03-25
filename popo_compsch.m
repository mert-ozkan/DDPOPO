clc; clear all; sca;
 
AssertOpenGL;
 
path.exp = '/Users/mertozkan/Documents/MATLAB/DD/DDPOPO';
path.mov = './movies';
path.dat = './data';
addpath(genpath(path.exp));
cd(path.exp)
%% Open Reaction Log
rxnlog.ExpID = 'POPO_CompSch';
rxnlog.createfolder = 0;
subinit = input('Initials: ','s');
if strcmp(subinit,'try')
    subinit = 'test';
end
sxn_no = input('Session number: ');
rxnlog.ExpID = sprintf('%s_S%d',rxnlog.ExpID,sxn_no);
rxnlog.dir = sprintf('%s/%s',path.dat,subinit);
while ~isfolder(rxnlog.dir)
    disp('The participants folder does not exist.');
    subinit = input('Initials: ','s');
    rxnlog.dir = sprintf('%s/%s',path.dat,subinit);
end

rxnlog = recrxn(rxnlog,'open',subinit);
%%
kb.esc = KbName('escape');
kb.space = KbName('space');
%% Open Window
scr_cfg.blendfunction = 'yes';
scr_cfg.sourcefactor = 'GL_SRC_ALPHA';
scr_cfg.destinationfactor = 'GL_ONE_MINUS_SRC_ALPHA';
scr_cfg.mode = 'openwindow';
scr_cfg.skipsynctests = 1;
scr_cfg.debugrect = 0;
scr_cfg.viewingdistance = 57;
scr_cfg.backgroundcolor = [127 127 127 255];

[scr]=openexperimentwindow(scr_cfg);
%% Parameters
ppd_h = scr.pixperdeg_h;
ppd_v = scr.pixperdeg_v;
cent_x = scr.xcenter;
cent_y = scr.ycenter;

env_sz = 1;
ext_spd = 3;
int_spd = 0;
distdefix = 8;
trl_dur = .5;
isReversal = 1;
intv_pre = [0.2 0.7];
intv_post = 0.3;

qFrm = round(trl_dur * scr.framerate);
intv.prestim = [0.3 0.5];
intv.poststim = 0.3;
fix_col = [255,255,255,255];


%% Pink noise texture
pink_cfg.mode = 'pinknoise';
pink_cfg.pixelsperdegree = ppd_h;
pink_cfg.framerate = scr.framerate;
pink_cfg.internalspeed = int_spd;
pink_cfg.trialduration = 1/scr.framerate; %(seconds)
pink_cfg.envelopesize = env_sz; %(degrees)
pink_cfg.backgroundcolour = 127;
pink_cfg.tar_frequency = .2;
pink_cfg.numberofoctaves = 3;
pink_cfg.randomnoise = 0;% = 1 only for the control stimulus without illusion
pink_cfg.reversal = isReversal; %(if 1, the stimulus ends at the initial location)

pink = createdoubledriftstimulus(pink_cfg);
pink_tex = Screen('MakeTexture', scr.windowptr,squeeze(pink.mat(:,:,1)));

qTrl = 1024;
setsize = 8;
trl_ord = settrialorder(qTrl,setsize);
%% Rect of the stimuli
initpt_angle = [pi/8:pi/4:2*pi-pi/8]+pi;

%% register parameters
recrxn(rxnlog,'line');
recrxn(rxnlog,'tab','Duration of one full path:'); recrxn(rxnlog,'line',trl_dur);
recrxn(rxnlog,'tab','External speed:'); recrxn(rxnlog,'line',ext_spd);
recrxn(rxnlog,'tab','Distance from foveal focus:'); recrxn(rxnlog,'line',distdefix);
recrxn(rxnlog,'tab','Set size:'); recrxn(rxnlog,'line',setsize);
recrxn(rxnlog,'line');

recrxn(rxnlog,'line');
recrxn(rxnlog,'tab','whBlock');
recrxn(rxnlog,'tab','whTrl');
recrxn(rxnlog,'tab','po_pos');
recrxn(rxnlog,'tab','po_dir');
recrxn(rxnlog,'tab','whCol');
recrxn(rxnlog,'tab','rxn_key');
recrxn(rxnlog,'line','rxn_t');
%% Start session
po(ceil(rand*8)) = 1;
po = po*pi;

dot_cols = [1,0,0,1;...
    .3,1,0,1]*180;

qBlk = 32;
whBlk = 0;
endsxn = 0;
whTrl = 1;
PsychHID('KbQueueCreate');
while ~endsxn
    
    %% Start the block
    if ~mod(whTrl-1,length(trl_ord)/qBlk)
        whBlk = whBlk+1;
        text = {sprintf('%d of %d blocks',whBlk,qBlk),'Press the spacebar to start the block.'};
        drawtexts(scr.windowptr,text,[cent_x-8*ppd_h;cent_x-8*ppd_h],[cent_y-5*ppd_v,cent_y],'flip',24);
        
        go = 0;
        PsychHID('KbQueueStart');
        PsychHID('KbQueueFlush');
        while ~go
            [keyisdown,keycode] = PsychHID('KbQueueCheck');
            WaitSecs(0.01);
            if keyisdown
                if keycode(kb.esc)
                    sca;
                    endsxn = 1;
                    break;
                elseif keycode(kb.space)
                    go = 1;
                end
            end
        end
        PsychHID('KbQueueStop');
        if endsxn, break; end
        
        text = 'Press the spacebar after you obtain fixation!';
        drawtexts(scr.windowptr,text,cent_x-8*ppd_h,cent_y-5*ppd_v,[],24);
        
        fix_col = [255 255 255 255];
        drawfixation(scr.windowptr,cent_x,cent_y,10,fix_col,'flip');

        go = 0;
        PsychHID('KbQueueStart');
        PsychHID('KbQueueFlush');
        while ~go
            [keyisdown,keycode] = PsychHID('KbQueueCheck');
            WaitSecs(0.01);
            if keyisdown
                if keycode(kb.esc)
                    sca;
                    endsxn = 1;
                    break;
                elseif keycode(kb.space)
                    go = 1;
                end
            end
        end
        PsychHID('KbQueueStop');
        if endsxn, break; end
    end
    
    %% Start the trial    
    % Pre stimulus interval
    secs = rand*(intv_pre(2)-intv_pre(1))+intv_pre(1);
    WaitSecs(secs)
    %% Place stimuli
    po = ones(1,8)*trl_ord(1,whTrl)/2;
    po(trl_ord(2,whTrl)) = (trl_ord(1,whTrl) == 1)/2 + 1/2;
    po = po*pi;
    
    colors = repmat([1,2],1,setsize/2);
    colors = colors(randperm(setsize));
    whCol = colors(trl_ord(2,whTrl));

    for whPatch = 1:8
        fname = sprintf('Patch%d',whPatch);

        rect_cfg.(fname).mode = 'linearrect';
        rect_cfg.(fname).pixelsperdegree = ppd_h;
        rect_cfg.(fname).framerate = scr.framerate;
        rect_cfg.(fname).distancefromthefixation = distdefix; %(degrees)
        rect_cfg.(fname).xcenter = scr.xcenter;
        rect_cfg.(fname).ycenter = scr.ycenter;
        rect_cfg.(fname).angleoftheinitialpointvector = initpt_angle(whPatch); %the angle of the vector drawn from
                                       %the origin to the initial position of
                                       %the stimulus.
        rect_cfg.(fname).quadrantlocation = 4; %(i.e. in which quadrant of the coordinate system
                           %the stimulus will appear)
        rect_cfg.(fname).trialduration = trl_dur; %(seconds)
        rect_cfg.(fname).reversal =  isReversal; %(if 1, the rect ends at the initial location)
        rect_cfg.(fname).horizontaldirection = 'none'; %('left','right','none')
        rect_cfg.(fname).verticaldirection = 'none'; %('up','down','none')
        rect_cfg.(fname).envelopesize = env_sz;%(degrees)
        rect_cfg.(fname).externalspeed = ext_spd;
        rect_cfg.(fname).externalmotionorientation = po(whPatch);
        rect_stim.(fname) = createdoubledriftstimulus(rect_cfg.(fname));
        
        for whFrm = 1:qFrm 
            rect_dot.(fname)(whFrm,:) = [ceil((rect_stim.(fname)(whFrm,3)-rect_stim.(fname)(whFrm,1))/2)+rect_stim.(fname)(whFrm,1),...
                ceil((rect_stim.(fname)(whFrm,4)-rect_stim.(fname)(whFrm,2))/2)+rect_stim.(fname)(whFrm,2)];
        end
        txt_x(whPatch) = rect_dot.(fname)(1,1);
        txt_y(whPatch) = rect_dot.(fname)(1,2);
    end
    

    
    
    %% Stimulus presentation
    isRxn = 0;
    whFrm = 1;
    onset = 0;
    PsychHID('KbQueueStart');
    PsychHID('KbQueueFlush');
    while ~isRxn
        
        whFrm = mod(whFrm-1,qFrm)+1;
        
        fix_col = [255 255 255 255];
        drawfixation(scr.windowptr,scr.xcenter,scr.ycenter,10,fix_col);
        for whPatch = 1:8
            fname = sprintf('Patch%d',whPatch);
            
            
            Screen('DrawTexture',scr.windowptr,...
                pink_tex,...
                [],...
                rect_stim.(fname)(whFrm,:));
            Screen('DrawDots', scr.windowptr,rect_dot.(fname)(whFrm,:), 3, dot_cols(colors(whPatch),:), [], 2);
        end
        t = Screen('Flip',scr.windowptr);
        
        if whFrm == 1 && onset == 0
            onset = t;
        end
        
        [keyisdown,keycode] = PsychHID('KbQueueCheck');
        if keyisdown
            if keycode(kb.esc)
                sca;
                endsxn = 1;
                break;
            elseif sum(keycode~=0)~=1
                rxn_key = '-1';
                rxn_t = 0;
                isRxn = 1;
                fix_col = [255 0 0 255];
            else
                rxn_key = KbName(find(keycode));
                if strcmp(rxn_key,'r')
                    rxn_key = '1'; %green
                elseif strcmp(rxn_key,'g')
                    rxn_key = '2';
                end
                rxn_t = GetSecs - onset;
                isRxn = 1;
                fix_col = [0 255 0 255];
            end                
        end    
        whFrm = whFrm + 1;
    end
    PsychHID('KbQueueStop');
    if endsxn, break; end
    
    drawfixation(scr.windowptr,cent_x,cent_y,10,fix_col,'flip')
   
    recrxn(rxnlog,'tab',whBlk);
    recrxn(rxnlog,'tab',whTrl);
    recrxn(rxnlog,'tab',trl_ord(2,whTrl));
    recrxn(rxnlog,'tab',trl_ord(1,whTrl));
    recrxn(rxnlog,'tab',whCol);
    recrxn(rxnlog,'tab',rxn_key);
    recrxn(rxnlog,'line',rxn_t);

    whTrl = whTrl + 1;
    
    WaitSecs(intv_post)
end


%% Complementary Functions
%%
function rxnlog = recrxn(rxnlog,command,var1)
 
if ~isfield(rxnlog,'dir')
    rxnlog.dir = '';
elseif ~strcmp(rxnlog.dir(end),'/')
    rxnlog.dir = [rxnlog.dir,'/'];
end
 
if ~isfield(rxnlog,'createfolder')
    rxnlog.createfolder = 0;
end
 
switch command
    case 'open'
        
        if ~isfield(rxnlog,'ExpID')
            rxnlog.ExpID = input('Provide an experiment identifier: ','s');
        end
        
        % Open a LogFile
        if nargin < 3
            filename = input('Name: ','s');
        else
            filename = var1;
        end
        
        if strcmp('try',filename)
                filename = 'test';
        end
        
        if rxnlog.createfolder
            if isfolder([rxnlog.dir,filename])
                disp('The folder already exists!');
            else
                mkdir([rxnlog.dir,filename]);
            end
            rxnlog.dir = [rxnlog.dir,filename,'/'];
        end
        
        if strcmp(filename,'test')
            filename = 'test';
        else
            adddate=clock;
            filename = [filename,'_',num2str(adddate(2)),'-',num2str(adddate(3)),'_',num2str(adddate(4)),'-',num2str(adddate(5))];
        end
        
        rxnlog.filename = [filename,'.txt'];
        rxnlog.id = fopen([rxnlog.dir,rxnlog.ExpID,'_',filename],'w');
        
        fprintf(rxnlog.id,'Experiment Info\n');
        fprintf(rxnlog.id,'%s\n',rxnlog.ExpID);
        fprintf(rxnlog.id,'File Name:\t%s\n',rxnlog.filename);
    otherwise
        switch command
            case 'line'
                writein = '\n';
            case 'tab'
                writein = '\t';
        end
        if nargin<3
            fprintf(rxnlog.id,writein);
        else
            var_class = class(var1);
            switch var_class
                case 'char'
                    notation = '%s';
                case 'double'
                    notation = '%d';
            end
            fprintf(rxnlog.id,[notation,writein],var1);
        end
end
end
%%
function drawtexts(w,text,x,y,flip,sz,col)
if nargin<7
    col = [255 255 255 255];
end
if nargin<6
    sz=100;
end
if nargin<5
    flip='no';
end
 
Screen('TextSize',w,sz)
 
if length(x)>1
    for whText = 1:length(x)
        Screen('DrawText',w,text{whText},x(whText),y(whText),col);
    end
else
    Screen('DrawText',w,text,x,y,col);
end
 
if strcmp(flip,'flip')
    Screen('Flip',w);
end
end
%% Set trial order
function trl_ord = settrialorder(qTrl,setsize)

while mod(qTrl,setsize)
    qTrl = qTrl + 1;
end

if mod(qTrl,2)
    qTrl = qTrl + 1;
end



cond1 = [ones(1,qTrl/2), ones(1,qTrl/2)*2];
cond2 = repmat(1:setsize,1,qTrl/setsize);

trl_ord = [cond1; cond2];
trl_ord = trl_ord(:,randperm(length(trl_ord)));

end
%% Open window
function [op]=openexperimentwindow(cfg)
%% cfg.mode = 'openwindow';
%   cfg.skipsynctests = 1;
%   cfg.debugrect = 0;
%   cfg.backgroundcolor = [127 127 127 255];
%   cfg.blendfunction = 'no';
%   cfg.sourcefactor
%   cfg.destinationfactor
%       % e.g. % Screen('BlendFunction', window, 'GL_ONE', 'GL_ONE');
%              % Screen('BlendFunction', window, 'GL_DST_ALPHA', 'GL_ONE_MINUS_DST_ALPHA');
%              % Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
%   cfg.viewingdistance
%     op.windowptr = window;
%     op.windowrect = windowRect;
%     op.widthinpix = wid_inpix;
%     op.heightinpix = height_inpix;
%     op.widthincm = wid_incm;
%     op.heightincm = height_incm;
%     op.xcenter = xCenter;
%     op.ycenter = yCenter;
%     op.framerate = framerate;
%     op.pixperdeg_h
%     op.pixperdeg_v
%%
switch cfg.mode
    case 'openwindow'
        PsychDefaultSetup(1);
        commandwindow;
        Screen('Preference', 'SkipSyncTests',cfg.skipsynctests);
        
        
        screens = Screen('Screens');
        screenNumber = max(screens);
        
        if cfg.debugrect == 1
            rect =  [800 400 1440 850];
        else
            rect = [];
        end
        
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
        
        
        [window, windowRect] = PsychImaging('OpenWindow',...
            screenNumber, cfg.backgroundcolor,rect);
        [wid_inpix, height_inpix] = Screen('WindowSize', window); % sz in px
        [xCenter, yCenter] = RectCenter(windowRect);
        [wid_incm, height_incm] = Screen('DisplaySize', screenNumber); % phys sz in mm
        
        framerate = FrameRate(window);
        
        switch cfg.blendfunction
            case 'yes'
                try
                    Screen('BlendFunction', window, cfg.sourcefactor, cfg.destinationfactor);
                catch
                    error('Define source and destination factors for the blending function');
                end
            otherwise
                disp('No blending function is being used.');
        end
        
        [ppd_h,ppd_v] = deg2pix(1,cfg.viewingdistance,window);
        
        op.windowptr = window;
        op.windowrect = windowRect;
        op.widthinpix = wid_inpix;
        op.heightinpix = height_inpix;
        op.widthincm = wid_incm;
        op.heightincm = height_incm;
        op.xcenter = xCenter;
        op.ycenter = yCenter;
        op.framerate = framerate;
        op.pixperdeg_h = ppd_h;
        op.pixperdeg_v = ppd_v;
    otherwise
        disp('Not developed yet!');
end
end
%% Function that converts degrees to pixels
function  [qPix_h,qPix_v] = deg2pix(deg,dist,windowpointer,roundthenumber)
if nargin < 3
    windowpointer = 0;
    warning('IF THE WINDOW POINTER IS OTHER THAN "0" PLEASE SPECIFY!');
end
if nargin < 4
    roundthenumber = 1;
    warning('The number of pixels are rounded to an integer. If you want otherwise the 4th input should be "0"');
end

screens = Screen('Screens');
screen_no = max(screens);

[wid_inpix,height_inpix] = Screen('WindowSize', windowpointer);
[wid_incm, height_incm] = Screen('DisplaySize', screen_no);

qPix_h = dist*tan(deg*pi/180)/...
    (wid_incm/(10*wid_inpix));

qPix_v = dist*tan(deg*pi/180)/...
    (height_incm/(10*height_inpix));

if roundthenumber
    qPix_h = round(qPix_h);
    qPix_v = round(qPix_v);
end
end
%%
function op = createdoubledriftstimulus(cfg)
%%
%% cfg.mode = 'linearrect'
%   cfg.pixelsperdegree
%   cfg.framerate
%   cfg.distancefromthefixation %(degrees)
%   cfg.xcenter
%   cfg.ycenter
%   cfg.angleoftheinitialpointvector %the angle of the vector drawn from
%                                    %the origin to the initial position of
%                                    %the stimulus.
%   cfg.quadrantlocation %(i.e. in which quadrant of the coordinate system
%                        %the stimulus will appear)
%   cfg.trialduration %(seconds)
%   cfg.reversal %(if 1, the rect ends at the initial location)
%   cfg.horizontaldirection %('left','right','none')
%   cfg.verticaldirection %('up','down','none')
%   cfg.envelopesize %(degrees)
%   cfg.externalspeed
%   cfg.externalmotionorientation;
%
%   "op" is the (number of frames x 4) destination rect for each frame;
%% cfg.mode = 'pinknoise' or '1/f'
%   cfg.pixelsperdegree
%   cfg.framerate
%   cfg.internalspeed
%   cfg.trialduration %(seconds)
%   cfg.envelopesize %(degrees)
%   cfg.backgroundcolour
%   cfg.tar_frequency
%   cfg.numberofoctaves
%   cfg.randomnoise % = 1 only for the control stimulus without illusion
%   cfg.reversal %(if 1, the stimulus ends at the initial location)
%
%   "op.mat" is the stimulus matrix (pixels x pixels x frames)
%   "op.rect' is the stimulus rect

%% cfg.mode = 'gabor'
%  cfg.function = 'create'
%   cfg.pixelsperdegree
%   cfg.spatialfrequency
% 	cfg.contrast % = .5;
% 	cfg.aspectratio % = 1;
% 	cfg.bgoffset % = [0 0 0 1];
% 	cfg.gaborwidth
% 	cfg.gaborheight
% 	cfg.disablenorm % = 1;
% 	cfg.precontrasmultiplier % = 1;
% 	cfg.windowptr
% 	cfg.phase
% 	cfg.sigma
%
%   "op.ptr" gabor texture
%   "op.rect"
%   "op.propmat" properties matrix
%
% cfg.mode = 'gabor';
% cfg.function = 'updatepropmat'
%   cfg.propmat = op.propmat
%   cfg.phasefactor
% 	cfg.framerate
% 	cfg.internalcyclespersec
%
%   "op" propmat
%
% cfg.mode = 'gabor';
% cfg.function = 'creatematrix'
%     cfg.gaborsize;
%     cfg.pixelsperdegree;
%     cfg.trialduration;
%     cfg.framerate;
%     cfg.angle;
%     cfg.sigma
%     cfg.cyclesperdeg;

%% cfg.mode = 'gratingmask';
%  cfg.function = 'create';
%     cfg.pixelsperdegree 
%     cfg.gratingwidth
%     cfg.gratingheight
%     cfg.spatialfrequency
%     cfg.cyclesize
%     cfg.bgcolor
%     cfg.contrast

%  cfg.mode = 'gratingmask';
%  cfg.function = 'updaterect';
%     cfg.pixelsperdegree 
%     cfg.framenumber
%     cfg.cyclesize
%     cfg.cyclespersecond
%     cfg.framerate
%     cfg.gratingwidth
%     cfg.gratingheight

%% cfg.mode = 'perlinnoise'
%     cfg.size;
%     cfg.pixelsperdegree

%% cfg.mode = 'gaussianmask';
%     cfg.pixelsperdegree
%     cfg.stimulusmatrix
%     cfg.sigma
%%
switch cfg.mode
    case 'linearrect'
        
        ppd = cfg.pixelsperdegree;
        frmrate = cfg.framerate;
        dist_decent = cfg.distancefromthefixation;
        th_decent = cfg.angleoftheinitialpointvector;
        whQuad = cfg.quadrantlocation;
        qFrm = round(cfg.trialduration*cfg.framerate);
        reversal = cfg.reversal;
        dir_h = cfg.horizontaldirection;
        dir_v = cfg.verticaldirection;
        env_sz = round(cfg.envelopesize*cfg.pixelsperdegree);
        env_sz = env_sz+(mod(env_sz,2)==0);
        rect_env = [0 0 env_sz env_sz];
        ext_spd = cfg.externalspeed;
        ext_or = cfg.externalmotionorientation;
        cent_x = cfg.xcenter;
        cent_y = cfg.ycenter;
        
        dplc_perfrm = displaceperfrm(ext_spd,frmrate);
        [dplc_perfrm_h,dplc_perfrm_v] = length_inxy(dplc_perfrm,ext_or,ppd,ppd);
        
        
        
        % rows: 1st is left; 2nd is right; 1st is up; 2nd is down
        if reversal
            drift_h(1,:) = -1*[0:dplc_perfrm_h:(qFrm/2)*dplc_perfrm_h (qFrm/2)*dplc_perfrm_h-1:-1*dplc_perfrm_h:1];
            drift_v(1,:) = -1*[0:dplc_perfrm_v:(qFrm/2)*dplc_perfrm_v (qFrm/2)*dplc_perfrm_v-1:-1*dplc_perfrm_v:1];
        else
            drift_h(1,:) = -1*[0:dplc_perfrm_h:(qFrm)*dplc_perfrm_h];
            drift_v(1,:) = -1*[0:dplc_perfrm_v:(qFrm)*dplc_perfrm_v];
        end
        
        drift_h(2,:) = -1*drift_h(1,:);
        drift_v(2,:) = -1*drift_v(1,:);
        
        if isempty(drift_v)
            drift_v = zeros(size(drift_h));
        end
        
        if isempty(drift_h)
            drift_h = zeros(size(drift_v));
        end

        [x,y]= length_inxy(dist_decent,th_decent,ppd,ppd);
        dist_decent_x = x + max(abs(drift_h(1,:)));
        dist_decent_y = y + max(abs(drift_v(1,:)));
        
        signx = 2*(whQuad==2 || whQuad==4) - 1;
        signy = 2*(whQuad==3 || whQuad==4) - 1;
        
        whH = 1+strcmp(dir_h,'right');
        whV = 1+strcmp(dir_v,'down');
        
        for whFrm = 1:qFrm
            pos_x = cent_x +...
                signx * dist_decent_x +...
                drift_h(whH,whFrm);
            pos_y = cent_y +...
                signy * dist_decent_y +...
                drift_v(whV,whFrm);
            rect_dst(whFrm,:) = CenterRectOnPointd(...
                rect_env, ...
                pos_x,...
                pos_y);
        end
        
        op = rect_dst;
        
    case {'pinknoise','1/f'}
        ppd = cfg.pixelsperdegree;
        frmrate = cfg.framerate;
        int_spd = cfg.internalspeed;
        qFrm = round(cfg.trialduration*cfg.framerate);
        env_sz = round(cfg.envelopesize*cfg.pixelsperdegree);
        env_sz = env_sz+(mod(env_sz,2)==0);
        gray = cfg.backgroundcolour;
        tar_fq = cfg.tar_frequency*cfg.pixelsperdegree;
        qOct = cfg.numberofoctaves;
        randomnoise = cfg.randomnoise;
        isReversal = cfg.reversal;
        
        stim_mat = pinknoisetexture(ppd,int_spd,1/frmrate,qFrm,env_sz,gray(1),tar_fq,qOct,isReversal,randomnoise);
        [sz,~,~] = size(stim_mat);
        op.mat = stim_mat;
        op.rect = [0 0 sz sz];
    case 'gabor'
        switch cfg.function
            case 'create'
                ppd = cfg.pixelsperdegree;
                spat_fq = cfg.spatialfrequency/cfg.pixelsperdegree;
                cont = cfg.contrast;
                asp_rat = cfg.aspectratio;
                bg_offset = cfg.bgoffset;
                wid = cfg.gaborwidth;
                wid = round(wid*ppd);
                wid = wid + (mod(wid,2)==0);
                height = cfg.gaborheight;
                height = round(height*ppd);
                height = height + (mod(height,2)==0);
                dblnorm = cfg.disablenorm;
                precont_multplr = cfg.precontrasmultiplier;
                w = cfg.windowptr;
                phase = cfg.phase;
                sig = cfg.sigma*cfg.pixelsperdegree;

                [gabor.ptr, gabor.rect] = CreateProceduralGabor(w, ...
                    wid,...
                    height,...
                    [],...
                    bg_offset,...
                    dblnorm,...
                    precont_multplr);

                gabor.propmat = [phase,...
                    spat_fq,...
                    sig,...
                    cont,...
                    asp_rat,...
                    0, 0, 0];

                op = gabor;
            case 'creatematrix'
                % Marvin's code
                env_sz = cfg.gaborsize;
                ppd = cfg.pixelsperdegree;
                trl_dur = cfg.trialduration;
                temporal_freq = 1/trl_dur;
                RefreshRate = cfg.framerate;
                sigma = cfg.sigma*ppd;
                cycperdeg = cfg.cyclesperdeg;
                cycperdeg2 = cycperdeg+4;

                gaborSizePix = env_sz *ppd;
                if ~mod(gaborSizePix,2)
                    gaborSizePix = gaborSizePix + 1;
                end
                
                halfgaborsize=floor(gaborSizePix/2);
                spatfreq = cycperdeg/ppd; % cycles per deg in cycles per pix
                pixelsPerCycle = ceil(1/spatfreq);
                intDriftCycPerFrame = temporal_freq/RefreshRate; % this fraction of a cycle per frame 
                intDriftPixPerFrame = ceil(pixelsPerCycle * intDriftCycPerFrame);
                framesPerCyc = ceil(RefreshRate/temporal_freq);
                % Compute actual cosine grating:
                gratingNEW = repmat(sin(linspace(0,2*pi,pixelsPerCycle)),gaborSizePix,ceil(gaborSizePix/pixelsPerCycle));
                gratingNEW = [gratingNEW,gratingNEW(:,2:end),gratingNEW(:,2:end);
                              gratingNEW,gratingNEW(:,2:end),gratingNEW(:,2:end)]; % makes grating from -1 to 1
%                 gratingNEW = gratingNEW + 1; % now grating from 0 to 2
%                 gratingNEW = gratingNEW/2; % grating from 0 to 1
%                 gratingNEW = gratingNEW * (255); % grating from 0 to 255

                spatfreq2 = cycperdeg2/ppd; % cycles per deg in cycles per pix
                pixelsPerCycle2 = ceil(1/spatfreq2);
                intDriftCycPerFrame = temporal_freq/RefreshRate; % this fraction of a cycle per frame 
                intDriftPixPerFrame2 = ceil(pixelsPerCycle2 * intDriftCycPerFrame);
                framesPerCyc = ceil(RefreshRate/temporal_freq);
                % Compute actual cosine grating:
                gratingNEW2 = repmat(sin(linspace(0,2*pi,pixelsPerCycle2)),gaborSizePix,ceil(gaborSizePix/pixelsPerCycle2));
                gratingNEW2 = [gratingNEW2,gratingNEW2(:,2:end),gratingNEW2(:,2:end);
                              gratingNEW2,gratingNEW2(:,2:end),gratingNEW2(:,2:end)]; % makes grating from -1 to 1
%                 gratingNEW2 = gratingNEW2 + 1; % now grating from 0 to 2
%                 gratingNEW2 = gratingNEW2/2; % grating from 0 to 1
%                 gratingNEW = gratingNEW * (255); % grating from 0 to 255
                
                % make gabors 
                % gabor = nan(2*sourcerect_size+1, 2*sourcerect_size+1,4);
                [x,y] = meshgrid(-halfgaborsize:halfgaborsize,-halfgaborsize:halfgaborsize);
                gabor = nan(gaborSizePix, gaborSizePix,3);
%                 gaussian =round(255 * (1 - exp(-(x.^2)/(2*sigma.^2)-(y.^2)/(2*sigma.^2))));
                gaussian =exp(-(x.^2)/(2*sigma.^2)-(y.^2)/(2*sigma.^2));               
                ngabors = ceil(framesPerCyc);
                % gabortex = nan(1,ngabors);
                for texi = 1:ngabors
                    startpoint = mod(texi*intDriftPixPerFrame,pixelsPerCycle)+ 1;
%                     startpoint = (texi-1)*intDriftPixPerFrame + 1;
                    grating = gratingNEW(1: size(gaussian,1),...
                        startpoint:startpoint -1 + size(gaussian,1));
                    
                    startpoint2 = mod(texi*intDriftPixPerFrame2,pixelsPerCycle2)+ 1;
%                     startpoint = (texi-1)*intDriftPixPerFrame + 1;
                    grating2 = gratingNEW2(1: size(gaussian,1),...
                        startpoint2:startpoint2 -1 + size(gaussian,1));
                    
                    for layeri = 1:3
                        gabor(:,:,layeri) = rescale(grating + transpose(grating2));
                    end
                    
                    for layeri = 1:3
                        gabor(:,:,layeri) = gabor(:,:,layeri)*255;%.*  (gaussian * 255);%((-1*gaussian) + 255);
                    end
%                     opaqueThresh = 0;
%                     gaussian2 = ((-1*gaussian) + 255);
%                     gaussian2 = 0.*(gaussian2 < opaqueThresh) + gaussian2.*(gaussian2 >=opaqueThresh);
                    gabor(:,:,4) = gaussian*255;%gaussian2;

                %     gabors(:,:,:,texi) = gabor;
                %     gaborins(:,:,texi) = grating .*  ((-1*gaussian) + 255);
                    gaborstim(texi,:,:,:) = gabor;%Screen('MakeTexture', w, gabor);
                end
                
                op = gaborstim;
                
            case 'updatepropmat'
                propmat = cfg.propmat;
                phasefactor = cfg.phasefactor;
                cyc = cfg.internalcyclespersec;
                incperfrm = cyc*360/cfg.framerate;
                
                propmat(1) = propmat(1)...
                    + phasefactor ...
                    * incperfrm;
                
                op = propmat;
        end
    case 'gratingmask'
        switch cfg.function
            case 'create'
                ppd = cfg.pixelsperdegree;
                wid = round(cfg.gratingwidth*cfg.pixelsperdegree);
                wid = wid + (mod(wid,2)==0);
                height = round(cfg.gratingheight*cfg.pixelsperdegree);
                height = height + (mod(height,2)==0);
                spat_fq = cfg.spatialfrequency;
                pixpercyc = round(cfg.cyclesize*cfg.pixelsperdegree);
                fq = (1/pixpercyc)*2*pi;
                bgcol = cfg.bgcolor;
                cont = cfg.contrast;

                texsize = wid/2;
                x = meshgrid(-texsize:texsize + pixpercyc, 1);

                grating = bgcol * cos(fq*x) + bgcol;
                mask = ones(1, numel(x), 2) * bgcol;

                mask(:, :, 2)= grating .* cont;

                op = mask;
            case 'updaterect'
                ppd = cfg.pixelsperdegree;
                frameCounter = cfg.framenumber-1;
                pixpercyc = round(cfg.cyclesize*cfg.pixelsperdegree);
                shiftPerFrame = (cfg.cyclespersecond*pixpercyc)/cfg.framerate;
                wid = round(cfg.gratingwidth*cfg.pixelsperdegree);
                wid = wid + (mod(wid,2)==0);
                height = round(cfg.gratingheight*cfg.pixelsperdegree);
                height = height + (mod(height,2)==0);
                
                
                xoffset = mod(frameCounter * shiftPerFrame, pixpercyc);
                srcRect = [xoffset 0 xoffset+wid height];
                
                op = srcRect;
        end
    case 'perlinnoise'
        %from Liwei's ddspiral_nulling2
        
        m = round(cfg.size*cfg.pixelsperdegree);
        m = m + (mod(m,2)==0);
        
        
        s = zeros([m,m]);     % Prepare output image (size: m x m)
        w = m;
        i = 0;
        while w > 8
            i = i + 1;
            d = interp2(randn([m,m]), i-1, 'spline');
            s = s + i * d(1:m, 1:m);
            w = w - ceil(w/2 - 1);
        end
        s = (s - min(min(s(:,:)))) ./ ((max(max(s(:,:))) - min(min(s(:,:))))*4)+.375;%range .375 to .625
        
        % want a circle
        [x,y]=meshgrid(1:m,1:m);
        s((x-m/2).^2 + (y-m/2).^2 > (m/2)^2) = .5;
        
        op = s*255;
    case 'gaussianmask'
        ppd = cfg.pixelsperdegree;
        stim = cfg.stimulusmatrix;
        sig = cfg.sigma*ppd;
        [env_sz,~] = size(stim);
        
        imWidth = floor(env_sz/2);
        [gx,gy]=meshgrid(-imWidth:imWidth, -imWidth:imWidth);
        env = exp( -((gx.^2)+(gy.^2)) /(2*(sig)^2));
        
        op = stim.*env+127;
end
end
%% Draw Fixation
function drawfixation(w,x,y,sz,fix_col,flip,wait)
if nargin<6
    flip = 'dont';
    wait = 0;
end
if nargin<7
    wait = 0;
end
Screen('DrawDots', w,[x y], sz, fix_col, [], 2);
if strcmp(flip,'flip')
    Screen('Flip',w);
end
if wait
    WaitSecs(wait);
end
end
%%
function mov = recdisplay(rec,mov,command,var1,var2)
% mov.name
% mov.framerate
% mov.dir
% 1. command='create'
%       rec==true = 1
%       
% 2. command='record'
%       rec==true = 1
%       var1: window pointer
%       var2: duration of the frame to stay in the display in sec
%
% 3. command='finalize'
%       rec==true = 1

if rec
    if ~isfield(mov,'name')
        mov.name = ['genericmoviename',num2str(randi(10000)),'.avi'];
    else
        mov.name = [mov.name,'.avi'];
    end
    if ~isfield(mov,'framerate'),   mov.framerate = 60;    end
    
    switch command
        case 'create'
            if isfield(mov,'dir')
             cd(mov.dir)
            end
            mov.obj = VideoWriter(mov.name);
            mov.obj.FrameRate = mov.framerate;
        case 'record'            
            currFrame = Screen('GetImage', var1);
            if nargin<5
                var2 = 0;
                repeatframe = 1;
            else
                repeatframe = round(var2*mov.framerate);
            end
            for whFrm = 1:repeatframe
                open(mov.obj);
                writeVideo(mov.obj,currFrame);
            end
        case 'finalize'
            close(mov.obj)        
    end
end
end
%%
%%
function [stim] = pinknoisetexture(ppd,int_spd,flp_intv,qFrm,env_sz,bg_col,tar_fq,qOct,isReversal,isControl)

if nargin<9
    isControl = 0;
end

if isControl
    stepf = round(ppd*(int_spd*flp_intv)/2);
    img = zeros(env_sz, env_sz, qFrm+2);
else
    stepf = round(ppd*(int_spd*flp_intv));
    img = zeros(env_sz*2+stepf*qFrm, env_sz);
end

pink = (255 * makepinknoise(img, tar_fq, qOct,isControl,stepf))  - bg_col;
[stim] = framesIllusion(stepf, env_sz, bg_col, tar_fq, qFrm, pink,isReversal,isControl);
end
%%
function img = makepinknoise(img, w, qOct,isControl,step,persistence, lacunarity)
% creates and sum successive noise functions, each with higher
% frequency and lower amplitude. nomalize output within [0 1]
%
% input:
% - im: initial matrix
% - wl: grid size in pixels of the first noise octave (lowest spatial frequency)
% - octaves: number of noisy octaves with increasing spatial frequency to
%            be added
%
% optional:
% - lacunarity: frequency multiplier for each octave (usually set to 2 so
%               spatial frequency doubles each octave)
% - persistence: amplitude gain (usually set to 1/lacunarity)
%
% When lacunarity=2 and persistence=0.5 you get ~ 1/f noise
%
if nargin <= 5
    lacunarity = 2;
    persistence = 0.5;
end
if isControl
    [n, m, v] = size(img);
    a = 1;
    for oct = 1:qOct
        rndim = -1 +2*rand(ceil(n/w),ceil(m/w),step*ceil(v/w));   % uniform
        [Xq,Yq,Zq] = ndgrid(linspace(1,size(rndim,2),m),linspace(1,size(rndim,1),n),linspace(1,size(rndim,3),v));
        d = interp3(rndim,Xq,Yq,Zq, 'cubic');
        img = img + a*d(1:n, 1:m, 1:v);
        a = a*persistence;
        w = w/lacunarity;
    end
    img = (img - min(min(min(img(:,:,:))))) ./ (max(max(max(img(:,:,:)))) - min(min(min(img(:,:,:)))));
else
    
    [n, m] = size(img);
    a = 1;
    for oct = 1:qOct
        rndim = -1 + 2*rand(ceil(n/w),ceil(m/w)); % uniform
        [Xq,Yq] = meshgrid(linspace(1,size(rndim,2),m),linspace(1,size(rndim,1),n));
        d = interp2(rndim,Xq,Yq, 'cubic');
        img = img + a*d(1:n, 1:m);
        a = a*persistence;
        w = w/lacunarity;
    end
    img = (img - min(min(img(:,:)))) ./ (max(max(img(:,:))) - min(min(img(:,:))));
end
end
%%
function[m] = framesIllusion(stepf,env_sz,bg_col, tar_fq, nFrames, noiseimg,isReversal, control)

% here we prepare one set texture with drifting internal motion (2D noise)
% - this version works with the perceptual task left/right tilt -

if isReversal
    reversal = nFrames/2;
else
    reversal = nFrames;
end

% gaussian envelope
imWidth = floor(env_sz/2);
[gx,gy]=meshgrid(-imWidth:imWidth, -imWidth:imWidth);
env = exp( -((gx.^2)+(gy.^2)) /(2*(tar_fq)^2));
m = zeros(env_sz, env_sz, nFrames);

if control
    c = 3;
    for fi=1:nFrames
        if fi>1
            if fi<=reversal
                c = c+1; %step;
            else
                c = c-1; %step;
            end
        end
        noisePatt = noiseimg(:,:,c);
        m(:,:,fi) = uint8(bg_col + noisePatt.*env);
    end
else
    segBeg = 1;
    segEnd = env_sz;
    segBeg2 = 1 + stepf*reversal;
    segEnd2 = env_sz + stepf*reversal;
    
    % compute textures for individual frames
    cf = 0; cb = 0; fi = 0;
    for i=1:nFrames
        if i<=reversal
            aBeg = segBeg + (cf*stepf);
            aEnd = segEnd + (cf*stepf);
            cf = cf+1;
        else
            aBeg = segBeg2 - (cb*stepf);
            aEnd = segEnd2 - (cb*stepf);
            cb = cb+1;
        end
        fi = fi + 1;
        noisePatt = noiseimg(aBeg:aEnd,:);
        m(:,:,fi) = uint8(bg_col + noisePatt.*env);
    end
end
end
%% Function that computes displacement units per frame:
function displace_perfrm = displaceperfrm(displace_persec,frmrate)
if nargin < 2
    frmrate = FrameRate(0);
    warning('If the window pointer is different than "0" introduce frame rate as the second input! Or else, the function FrameRate might interfere with your display')
end
displace_perfrm = displace_persec/frmrate;
end
%% Function that gives pixel size in x and y coordinates separately for any vector with a certain angle:
function [qPix_h,qPix_v] = length_inxy(length_indeg,theta,h_pix_perdeg,v_pix_perdeg)

% theta = mod(theta,pi/2);
qPix_h = round(h_pix_perdeg * length_indeg * cos(theta));
qPix_v = round(v_pix_perdeg * length_indeg * sin(theta));

if theta == 0
    qPix_v = 0;
elseif theta == pi/2
    qPix_h = 0;
end
end