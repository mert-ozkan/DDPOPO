function [op, trl_sq, init_trl] = discriminate_ddpopo(op, trl_sq, init_trl, int_spdX)
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
dist_deFix = 8;
trl_dur = .31;

ext_th = [45, 135];
int_th = [315, 135;...
    225, 45]; % Rows correspond to ext_th; Columns correspond to horizontal and vertical illusion, respectively.
int_th = convert_toPtbTh(int_th);

intv_pre = [0.2 0.7];
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
pink_cfg.trialduration = trl_dur; %(seconds)
pink_cfg.envelopesize = env_sz; %(degrees)
pink_cfg.backgroundcolour = 127;
pink_cfg.tar_frequency = .2;
pink_cfg.numberofoctaves = 3;
pink_cfg.randomnoise = 0;% = 1 only for the control stimulus without illusion
pink_cfg.reversal = 0; %(if 1, the stimulus ends at the initial location)

%% Display
qTrl = length(trl_sq);
whTrl = init_trl;
isEndSxn = false;
PsychHID('KbQueueCreate');
while whTrl <= qTrl && ~isEndSxn
    %% Parameters
    stim_posX = trl_sq(whTrl,[1, 3]);
    trj_tipX = trl_sq(whTrl,[2, 4]);
    ill_tipX = ones(1,2)*trl_sq(whTrl,5);
    if sum(ill_tipX) == false
        ill_tipX = mod([1,0]-randi(2),2)+1;
    end
    int_thX = [int_th(trj_tipX(1),ill_tipX(1)), int_th(trj_tipX(2),ill_tipX(2))];
    
    coordN = squeeze(coordX(stim_posX,trj_tipX,:,:)); 
    %% Create stimulus
    pink_cfg.internalspeed = int_spdX(stim_posX(1),trj_tipX(1),ill_tipX(1));
    pink1 = createdoubledriftstimulus(pink_cfg);
    [~,~,qFrm] = size(pink1.mat);
    for currFrm = 1:qFrm
        pink1.tex(currFrm) = Screen('MakeTexture', scr.windowptr,squeeze(pink1.mat(:,:,currFrm)));
    end

    pink_cfg.internalspeed = int_spdX(stim_posX(2),trj_tipX(2),ill_tipX(2));
    pink2 = createdoubledriftstimulus(pink_cfg);
    [~,~,qFrm] = size(pink2.mat);
    for currFrm = 1:qFrm
        pink2.tex(currFrm) = Screen('MakeTexture', scr.windowptr,squeeze(pink2.mat(:,:,currFrm)));
    end
    %% Initial screen
    if whTrl == init_trl
        isEndSxn = initial_screen(scr,whTrl,qTrl);
        if isEndSxn, break; end
    end
    %% Jittered Screen
    draw_dots(scr.windowptr,scr.xcenter,scr.ycenter,'Flip','WaitSecs',unifrnd(intv_pre(1),intv_pre(2)));
    %% Start trial
    currFrm = ceil(qFrm/2);
    
    PsychHID('KbQueueStart'); PsychHID('KbQueueFlush');
    t0 = GetSecs;
    while GetSecs-t0 <= .5   
        whFrm = get_currFrmNo(currFrm,qFrm,true);
        draw_dots(scr.windowptr,scr.xcenter,scr.ycenter);
        draw_DDpatches(scr,coordN(1,:,:,:),whFrm,pink1,int_thX(1),trj_tipX(1));
        Screen('Flip',scr.windowptr);
        
        [isEndSxn, ~, ~, ~] = reaction;
        if isEndSxn, break; end
                  
        currFrm = currFrm + 1;
    end
    if isEndSxn, break; end
    
    draw_dots(scr.windowptr,scr.xcenter,scr.ycenter,'Flip','WaitSecs',0.2);
    
    t0 = GetSecs;
    while GetSecs-t0 <= .5   
        whFrm = get_currFrmNo(currFrm,qFrm,true);
        draw_dots(scr.windowptr,scr.xcenter,scr.ycenter);
        draw_DDpatches(scr,coordN(2,:,:,:),whFrm,pink2,int_thX(2),trj_tipX(2));
        Screen('Flip',scr.windowptr);
        
        [isEndSxn, ~, ~, ~] = reaction;
        if isEndSxn, break; end
                  
        currFrm = currFrm + 1;
    end
    if isEndSxn, break; end
    PsychHID('KbQueueStop');
    
    draw_dots(scr.windowptr,scr.xcenter,scr.ycenter,'Flip');
    [isEndSxn, isRxn, rxn_key, ~] = wait_forReaction({'f','j'});
    if isEndSxn, break; end
    
    isSame = false;
    if strcmp(rxn_key,'f')
        isSame = true;
    end
    
    op(whTrl,:) = [stim_posX, trj_tipX, ill_tipX, isSame];
    whTrl = whTrl + 1;
end
sca
init_trl = whTrl + isRxn; % if quit before completing the next trial, next session will start from the init_trl
end
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
function isEndSxn = initial_screen(scr,whTrl,qTrl)
hor_ppd = scr.pixperdeg_h;
txtX = {...
    sprintf('You will start this session by Trial %d/%d.',whTrl,qTrl),...
    'Press the SPACE BAR in order to PROCEED',...
    'If you want to QUIT the session, press the ESCAPE key. You will continue from this block in the next session.'};
txt_x = ones(1,length(txtX))*scr.xcenter-(5*hor_ppd);
txt_y = ones(1,length(txtX))*scr.ycenter-([5:-1:5-length(txtX)+1]*hor_ppd);
draw_texts(scr.windowptr,txtX,txt_x,txt_y);
draw_dots(scr.windowptr,scr.xcenter,scr.ycenter,'Flip');

[isEndSxn, ~, ~, ~] = wait_forReaction('space');
if isEndSxn, return; end

txtX = {...
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