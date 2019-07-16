function run_discriminate_ddpopo(varargin)
clear all;
close all;
sca;
%%
addpath('/Users/mertozkan/Documents/MATLAB/PTB');
addpath(genpath('/Users/mertozkan/Documents/MATLAB/DD/DDPOPO'));
dr_spd = '/Users/mertozkan/Documents/MATLAB/DD/DDPOPO/data/internal_speeds';
dr_dat = '/Users/mertozkan/Documents/MATLAB/DD/DDPOPO/data';
%%
init = input('Initials: ','s');
f_nm = sprintf('DDPOPO_Dcr_%s',init);
if strcmp('test',init) || strcmp('try',init)
    f_nm = 'test';    
end
f_nm = sprintf('%s/%s.mat',dr_dat,f_nm);

if isfile(f_nm) && ~strcmp(f_nm,'test.mat')
    load(f_nm);
else
    trl_sq = setTrialOrder(12);
    init_trl = 1;
    op = [];
end

int_spdX = return_speeds(init,dr_spd);

if ~(init_trl > length(trl_sq))
    [op, trl_sq, init_trl] = discriminate_ddpopo(op, trl_sq, init_trl, int_spdX);
else
    error('Data collection for this participant seem to have completed.')
end
sca;

save(f_nm, 'op', 'trl_sq', 'init_trl')
if init_trl > length(trl_sq)
    disp('Quest session finalised');
end
end
%%
function int_spdX = return_speeds(init,dr)
if strcmp(init,'test') || strcmp(init,'try')
    int_spdX = ones(8,2,2)*4;
else
    f = sprintf('%s/DDPOPO_Quest_%s.mat',dr,init);
    if isfile(f)
        load(f,'op')
        for whPos = 1:8
            for whTrj = 1:2
                for whIll = 1:2
                    field_nm = sprintf('Stim%d%d%d',whPos,whTrj,whIll);
                    int_spdX(whPos,whTrj,whIll) = op.(field_nm).avg;
                end
            end
        end
    else
        error('There is no file under this name.');
    end
end
end
%%
function trlX = setTrialOrder(qTrl)
trlX = repmat(combvec(1:8,1:2,1:2,1:2),1,qTrl)';
trlX = trlX(randperm(length(trlX)),:);
trlX = [trlX(:,1:2), repmat(1:8,1,length(trlX)/8)', trlX(:,3:4)];
trlX = trlX(randperm(length(trlX)),:);

qTestTrlX = length(trlX);
ctch_trlX = randperm(qTestTrlX,qTestTrlX/3);
trlX = [trlX; trlX(ctch_trlX,:)];
trlX(qTestTrlX+1:end,5) = 0;
trlX = trlX(randperm(length(trlX)),:);
end