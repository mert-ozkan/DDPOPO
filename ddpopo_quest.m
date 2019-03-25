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
hor_ppd = scr.pixperdeg_h;
ver_ppd = scr.pixperdeg_v;

qPatch = 8;
env_sz = 1;
int_spd = 4;
dist_deFix = 8;
trl_dur = .31;

ext_th = [135, 45];
int_th = [315, 135;...
    225, 45]; % Rows correspond to ext_th; Columns correspond to horizontal and vertical illusion, respectively.
int_th = convert_toPtbTh(int_th);

intv_pre = [0.2 0.7];

qFrm = round(trl_dur * scr.framerate);
%% Stimulus trajectories
rect_cfg.numberoftrajectories = qPatch;
rect_cfg.trialduration = trl_dur;
rect_cfg.framerate = scr.framerate;
rect_cfg.degrees_persecond = 3;
rect_cfg.horizontalpixels_perdegree = hor_ppd;
rect_cfg.verticalpixels_perdegree = ver_ppd;
rect_cfg.equidistant_point = 'middle';
rect_cfg.distancefromcenter = dist_deFix;
rect_cfg.xcenter = scr.xcenter;
rect_cfg.ycenter = scr.ycenter;
rect_cfg.angleoftrajectory = ext_th;
coordX = equidistant_trajectory_coordinates(rect_cfg);

%% Pink noise texture
pink_cfg.mode = 'pinknoise';
pink_cfg.pixelsperdegree = hor_ppd;
pink_cfg.framerate = scr.framerate;
% pink_cfg.internalspeed = int_spd;
pink_cfg.trialduration = trl_dur; %(seconds)
pink_cfg.envelopesize = env_sz; %(degrees)
pink_cfg.backgroundcolour = 127;
pink_cfg.tar_frequency = .2;
pink_cfg.numberofoctaves = 3;
pink_cfg.randomnoise = 0;% = 1 only for the control stimulus without illusion
pink_cfg.reversal = 0; %(if 1, the stimulus ends at the initial location)

%% Display
whPos = 1;
whPhysTrj = 1;
whIllTrj = 1;
for whTrl = 1:40    
    %% Initial screen
    if whTrl == 1
        isEndSxn = initial_screen(scr,whBlk,qTrl);
        if isEndSxn, break; end
    end
    %% Create stimulus
    
    %% Jittered Screen
    draw_dots(scr.windowptr,scr.xcenter,scr.ycenter,'Flip','WaitSecs',unifrnd(intv_pre(1),intv_pre(2)));
    %% Start trial
    currFrm = ceil(qFrm/2);
    isNxt = false;
    while ~isNxt        
        whFrm = get_currFrmNo(currFrm,qFrm,true);
        
        draw_dots(scr.windowptr,scr.xcenter,scr.ycenter);
        drawDDpatches(scr,coordX,whFrm,pink,int_thX,ext_thX);
        Screen('Flip',scr.windowptr);
    end
end
%%
%%
function whFrm = get_currFrmNo(currFrm,qFrm,isRev,onsetX)
if nargin < 3
    isRev = false;
end

if nargin < 4
    onsetX = 1;
end

if isRev    
    frm_noX = [1:qFrm qFrm-1:-1:2];
    mod_qFrm = qFrm*2-2;
else
    frm_noX = 1:qFrm;
    mod_qFrm = qFrm;
end

for idx = 1:length(onsetX)
    if currFrm >= onsetX(idx)
        idx_frm(idx) = mod(currFrm-onsetX(idx),mod_qFrm)+1;
    else
        idx_frm(idx) = 0;
    end
    
end

for idx = 1:length(onsetX)
    if idx_frm(idx)
        whFrm(idx) = frm_noX(idx_frm(idx));
    else
        whFrm(idx) = 0;
    end
end
end
%%
function isEndSxn = initial_screen(scr,whBlk,qTrl)
hor_ppd = scr.pixperdeg_h;
txtX = {...
    sprintf('This is the Block No%d.',whBlk),...
    'Press the SPACE BAR in order to PROCEED',...
    'If you want to QUIT the session, press the ESCAPE key. You will continue from this block in the next session.'};
txt_x = ones(1,length(txtX))*scr.xcenter-(5*hor_ppd);
txt_y = ones(1,length(txtX))*scr.ycenter-([5:-1:5-length(txtX)+1]*hor_ppd);
draw_texts(scr.windowptr,txtX,txt_x,txt_y);
draw_dots(scr.windowptr,scr.xcenter,scr.ycenter,'Flip');

[isEndSxn, ~, ~, ~] = wait_forReaction('space');

if isEndSxn, return; end
PsychHID('KbQueueStop');

txtX = {...
    sprintf('There will be %d trials in this block.',qTrl),...
    'You need to complete each trial without giving a break in between.',...
    'Press the SPACE BAR in order to PROCEED',...
    'In the NEXT SCREEN, press the SPACE BAR again after you make sure that you OBTAINED FIXATION.',...
    'If you want to QUIT the session, press the ESCAPE key. You will continue from this block in the next session.'};
txt_x = ones(1,length(txtX))*scr.xcenter-(5*hor_ppd);
txt_y = ones(1,length(txtX))*scr.ycenter-([5:-1:5-length(txtX)+1]*hor_ppd);
draw_texts(scr.windowptr,txtX,txt_x,txt_y);
draw_dots(scr.windowptr,scr.xcenter,scr.ycenter,'Flip');
[isEndSxn, ~, ~, ~] = wait_forReaction('space');
if isEndSxn, return; end

draw_dots(scr.windowptr,scr.xcenter,scr.ycenter,'Flip');
[isEndSxn, ~, ~, ~] = wait_forReaction('space');
if isEndSxn, return; end
end