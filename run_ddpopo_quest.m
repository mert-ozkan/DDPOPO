function run_ddpopo_quest(varargin)
close all;
sca;
%%
addpath('/Users/mertozkan/Documents/MATLAB/PTB');
addpath(genpath('/Users/mertozkan/Documents/MATLAB/DD/DDPOPO'));
dr_dat = '/Users/mertozkan/Documents/MATLAB/DD/DDPOPO/data';
%%
init = input('Initials: ');
f_nm = sprintf('DDPOPO_StimParam_%s',init);
if strcmp('test',init) || strcmp('try',init)
    f_nm = init;    
end
%%
qPos = 8;
qPhysTrj = 2;
qIllTrj = 2;
qTrl_perStim = 40;

end