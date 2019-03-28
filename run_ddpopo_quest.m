function run_ddpopo_quest(varargin)
close all;
sca;
%%
addpath('/Users/mertozkan/Documents/MATLAB/PTB');
addpath(genpath('/Users/mertozkan/Documents/MATLAB/DD/DDPOPO'));
dr_dat = '/Users/mertozkan/Documents/MATLAB/DD/DDPOPO/data';
%%
init = input('Initials: ','s');
f_nm = sprintf('DDPOPO_Quest_%s',init);
if strcmp('test',init) || strcmp('try',init)
    f_nm = init;    
end
f_nm = sprintf('%s/%s.mat',dr_dat,f_nm);
if isfile(f_nm)
    load(f_nm);
else
    trl_sq = setTrialOrder_inDDPOPOQuest(40);
    init_trl = 1;
    qst = easy_quest('Create',4,.1,'StimulusTypes',[8,2,2]);
end

[qst, trl_sq, init_trl] = ddpopo_quest(qst, trl_sq, init_trl);
sca
end
%%
function trlX = setTrialOrder_inDDPOPOQuest(qTrl_perStim)
trlX = repmat(combvec(1:8,1:2,1:2),1,qTrl_perStim);
trlX = trlX(:,randperm(length(trlX)))';
end