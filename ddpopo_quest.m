function [qst, trl_sq, init_trl, int_spdX] = ddpopo_quest(qst, trl_sq, init_trl, int_spdX)
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
    %% Create stimulus
    int_spd = easy_quest('Quantile',qst,'Condition',trl_sq(whTrl,:));
    pink_cfg.internalspeed = int_spd;
    pink = createdoubledriftstimulus(pink_cfg);
    [~,~,qFrm] = size(pink.mat);
    for currFrm = 1:qFrm
        pink.tex(currFrm) = Screen('MakeTexture', scr.windowptr,squeeze(pink.mat(:,:,currFrm)));
    end
    %% Parameters
    stim_pos = trl_sq(whTrl,1);
    trj_tipN = trl_sq(whTrl,2);
    ill_tipN = trl_sq(whTrl,3);
    int_thN = int_th(trl_sq(whTrl,2),ill_tipN);
    
    coord = squeeze(coordX(:,trj_tipN,:,:));
    coord = coord(stim_pos,:,:);
    
    %% Initial screen
    if whTrl == init_trl
        isEndSxn = initial_screen(scr,whTrl,qTrl);
        if isEndSxn, break; end
    end
    %% Jittered Screen
    draw_dots(scr.windowptr,scr.xcenter,scr.ycenter,'Flip','WaitSecs',unifrnd(intv_pre(1),intv_pre(2)));
    %% Start trial
    currFrm = ceil(qFrm/2);
    isNxt = false;
    PsychHID('KbQueueStart'); PsychHID('KbQueueFlush');
    th = randi(180)+.5;
    lil_step_th = 1;
    big_step_th = 10;
    while ~isNxt        
        whFrm = get_currFrmNo(currFrm,qFrm,true);
        draw_line(scr,'Slope',th);
        draw_dots(scr.windowptr,scr.xcenter,scr.ycenter);
        draw_DDpatches(scr,coord,whFrm,pink,int_thN,trj_tipN);
        Screen('Flip',scr.windowptr);
        
        [isEndSxn, ~, rxn_key, ~] = reaction({'UpArrow','DownArrow','LeftArrow','RightArrow','space'});
        if isEndSxn, break; end
        switch rxn_key
            case 'UpArrow'
                th = th + lil_step_th;
            case 'DownArrow'
                th = th - lil_step_th;
            case 'LeftArrow'
                th = th - big_step_th;
            case 'RightArrow'
                th = th + big_step_th;
            case 'space'
                rxn = th_to2AFC(th, trj_tipN, ill_tipN);
                isNxt = true;
        end
                
        currFrm = currFrm + 1;
    end
    if isEndSxn, break; end
    PsychHID('KbQueueStop');
    
    int_spdX(whTrl,:) = [stim_pos, trj_tipN, ill_tipN, int_spd, th, rxn];
    qst = easy_quest('Update',qst,int_spd,rxn,'Condition',trl_sq(whTrl,:));
    
    whTrl = whTrl + 1;
end

init_trl = whTrl+1; % if quit before completing the next trial, next session will start from the init_trl
end
%%
function rxn = th_to2AFC(th, trj_tip, ill_tip)

th = mod(mod(th,360),180);
doe = 10; % degree of error

if trj_tip == 1
    if ill_tip == 1
        rxn = th > 0 && th < (45 + doe);
    else
        rxn = th > (45 - doe) && th < 90;       
    end
else
    if ill_tip == 1
        rxn = th > 90 && th < (135 + doe);
    else
        rxn = th > (135 - doe) && th < 180;
    end
end

rxn = ~rxn;

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