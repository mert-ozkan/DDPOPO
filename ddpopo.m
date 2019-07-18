function [isEndSxn, isComplete, scr] = ddpopo(scr, f, whBlk, trl_sq, int_spdX)
%% Priming of Pop-Out in Double-Drift
%% Open Window
scr_cfg.blendfunction = 'yes';
scr_cfg.sourcefactor = 'GL_SRC_ALPHA';
scr_cfg.destinationfactor = 'GL_ONE_MINUS_SRC_ALPHA';
scr_cfg.mode = 'openwindow';
scr_cfg.skipsynctests = 1;
scr_cfg.debugrect = 0;
scr_cfg.viewingdistance = 57;
scr_cfg.backgroundcolor = [127 127 127 255];

if isempty(fieldnames(scr))
    [scr]=openexperimentwindow(scr_cfg);
end
%% Parameters
hor_ppd = scr.pixperdeg_h;
ver_ppd = scr.pixperdeg_v;

qPatch = 8;
env_sz = 1;
dist_deFix = 6;
trl_dur = .31;

ext_th = [45, 135];
int_th = [315, 135;...
    225, 45]; % Rows correspond to ext_th; Columns correspond to horizontal and vertical illusion, respectively.
int_th = convert_toPtbTh(int_th);

intv_pre = [0.2 0.7];

qFrm = round(trl_dur * scr.framerate);
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
%% Write variable names
write_inFile(f,',','whTrl');
write_inFile(f,',','cryp');
write_inFile(f,',','repCryp');
write_inFile(f,',','tar_clr');
write_inFile(f,',','rxn_key');
write_inFile(f,',','rxn_t');
write_inFile(f,',','prop1');
write_inFile(f,',','prop2');
write_inFile(f,',','prop3');
write_inFile(f,',','prop4');
write_inFile(f,',','prop5');
write_inFile(f,',','prop6');
write_inFile(f,',','prop7');
write_inFile(f,'line','prop8');

%% Display
isEndSxn = false;
whTrl = 1;
qTrl = length(trl_sq);
prevCryp = 99999999;
PsychHID('KbQueueCreate');
while ~isEndSxn && whTrl <= qTrl
    %% Trial parameters
    cryp = trl_sq(whTrl);
    tar_pos = return_digit(cryp,'last');
    
    repCryp = compare_digits(cryp,prevCryp);
    
    [ext_thX, int_thX, ill_tipX] = dcryp_patch_properties(cryp,qPatch,int_th);
    crypX = ill_tipX*10 + ext_thX;
    
    dot_clr = repmat([1;2],qPatch/2,1);
    dot_clr = dot_clr(randperm(qPatch));
    
    onsetX = asynchronous_onsets(scr.framerate,.5,qPatch);
    %% Create stimulus
    for whPatch = 1:qPatch
        pink_cfg.internalspeed = int_spdX(tar_pos, ext_thX(whPatch), ill_tipX(whPatch));
        pink = createdoubledriftstimulus(pink_cfg);
        for currFrm = 1:qFrm
            pink_tex(whPatch, currFrm) = Screen('MakeTexture', scr.windowptr,squeeze(pink.mat(:,:,currFrm)));
        end
    end
    pink.tex = pink_tex;
    %% Initial screen
    if whTrl == 1
        isEndSxn = initial_screen(scr,whBlk,qTrl);
        if isEndSxn, break; end
    end   
    %% Jittered Screen
    draw_dots(scr.windowptr,scr.xcenter,scr.ycenter,'Flip','WaitSecs',unifrnd(intv_pre(1),intv_pre(2)));
    %% Start trial
    currFrm = ceil(qFrm/2);
    isRxn = false;
    isFrm1 = true;
    PsychHID('KbQueueStart'); PsychHID('KbQueueFlush');
    while ~isRxn        
        whFrm = get_currFrmNo(currFrm,qFrm,true,onsetX);
        draw_dots(scr.windowptr,scr.xcenter,scr.ycenter);
        draw_DDpatches(scr,coordX,whFrm,pink,int_thX,ext_thX,'SuperimposeDots','AssignColor',dot_clr);
        flip_t = Screen('Flip',scr.windowptr);
        
        if isFrm1
            onset_t = flip_t;
        end
        isFrm1 = false;

        [isEndSxn, isRxn, rxn_key, rxn_t] = reaction({'r','g'});
        if isEndSxn, break; end
        
        currFrm = currFrm + 1;
        
    end
    if isEndSxn, break; end
    PsychHID('KbQueueStop');
    
    if strcmpi(rxn_key,'r')
        rxn_key = 1;
    else
        rxn_key = 2;
    end
    
    rxn_t = rxn_t - onset_t;
    
    write_inFile(f,',',whTrl);
    write_inFile(f,',',cryp);
    write_inFile(f,',',repCryp);
    write_inFile(f,',',dot_clr(tar_pos));
    write_inFile(f,',',rxn_key);
    write_inFile(f,',',rxn_t);
    write_inFile(f,',',crypX(1));
    write_inFile(f,',',crypX(2));
    write_inFile(f,',',crypX(3));
    write_inFile(f,',',crypX(4));
    write_inFile(f,',',crypX(5));
    write_inFile(f,',',crypX(6));
    write_inFile(f,',',crypX(7));
    write_inFile(f,'line',crypX(8));
    
    prevCryp = cryp;
    whTrl = whTrl + 1;
end
isComplete = true;
if whTrl-1 ~= qTrl
    isComplete = false;
end
write_inFile(f,'close');
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
%%
function comp = compare_digits(num1,num2)
st1 = num2str(num1);
st2 = num2str(num2);

if length(st1) ~= length(st2)
    error('The number of digits do not match.')
end

comp = sprintf('%d',st1==st2);

end
%%
function [ext_thX, int_thX, illX] = dcryp_patch_properties(cryp, qPatch, int_th)

tar_ill = return_digit(cryp,1);
tar_phys = return_digit(cryp,2);
tar_pos = return_digit(cryp,'last');

ext_thX = [ones(1,floor((qPatch-1)/2)+1).*((tar_phys==1)+1),ones(1,floor((qPatch-1)/2)).*tar_phys];
ext_thX = ext_thX(randperm(qPatch-1));
ext_thX = insert(ext_thX,tar_phys,tar_pos);

illX = ones(1,qPatch-1)*(tar_ill==1)+1;
illX = insert(illX,tar_ill,tar_pos);

int_thX = diag(int_th(ext_thX,illX));
end
%%
function vec = insert(vec, el, idx)
if idx > length(vec)+1
    error('Index is greater than the vector length.')
end

if idx == 0
    idx = 1;
end

if idx ~= 1
    vec = [vec(1:idx-1) el vec(idx:end)];
else
    vec = [el vec];
end

end
%%
function dig = return_digit(num,whDig)
st = num2str(num);
if strcmp(whDig,'last')
    whDig = length(st);
end
dig = str2double(st(whDig));
end
%%
function onsetX = asynchronous_onsets(frm_rate,secF,q)
if secF
    frmF = round(secF*frm_rate);
    onsetX = ceil(unifrnd(0,frmF,1,q));
    while ~ismember(1,onsetX) || ~ismember(frmF,onsetX)
        onsetX = randi(frmF,1,q);
    end
else
    onsetX = ones(1,q);
end
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
