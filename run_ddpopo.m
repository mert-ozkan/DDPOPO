function run_ddpopo(varargin)
close all;
sca;
%%
addpath('/Users/mertozkan/Documents/MATLAB/PTB');
addpath(genpath('/Users/mertozkan/Documents/MATLAB/DD/DDPOPO'));
dr_dat = '/Users/mertozkan/Documents/MATLAB/DD/DDPOPO/data';
dr_incomplete = '/Users/mertozkan/Documents/MATLAB/DD/DDPOPO/data/incomplete';
%%
whBlk = 1;
qBlk = 1;
isEndSxn = false;   

while whBlk <= qBlk && ~isEndSxn
    if length(varargin) == 1
        [f_nm, init, trl_sq_no, whBlk] = create_filename(dr_dat,'debug');
    elseif whBlk == 1
        [f_nm, init, trl_sq_no, whBlk] = create_filename(dr_dat);
    else
        [f_nm, init, trl_sq_no, whBlk] = create_filename(dr_dat, 'SessionInfo', init, trl_sq_no, whBlk);
    end
    f.directory = dr_dat;
    f = write_inFile(f,'open',f_nm);
    blkX = imp_blkX(sprintf('Trial_Sequence_No%d.txt',trl_sq_no));
    [qBlk,~] = size(blkX);
    trl_sq = blkX(whBlk,:);    
    [isEndSxn, isComplete] = ddpopo(f, whBlk, trl_sq);
    
    whBlk = whBlk + 1;
end
sca;

if ~isComplete && length(varargin) ~= 1
    movefile(f.name,dr_incomplete);
end
end
%%
function [f_nm, init, trl_sq_no, whBlk] = create_filename(dr,varargin)

debug = false;
isAskSxnInfo = true;

for idx = 1:length(varargin)
    argN = varargin{idx};
    switch argN
        case 'debug'
            debug = true;
        case 'SessionInfo'
            isAskSxnInfo = false;
            init = varargin{idx+1};
            trl_sq_no = varargin{idx+2};
            whBlk = varargin{idx+3};
    end
end  

if debug
    init = 'test';
    trl_sq_no = 1;
    whBlk = 1;
elseif isAskSxnInfo
    init = input('Initials: ','s');
    trl_sq_no = input('The number of the trial sequence: ');
    whBlk = input('The number of the block to start from: ');
end

f_nm = sprintf('DDPOPO_%s_TrlSqNo%d_BlkNo%d',init,trl_sq_no,whBlk);
prevFNm = sprintf('DDPOPO_%s_TrlSqNo%d_BlkNo%d',init,trl_sq_no,whBlk-1);

isVld = false;
if strcmp('test', init) || strcmp('try', init)
    f_nm = 'test';
    isVld = true;
end

if ~isVld
    fX = dir(dr);
    for idx = 1:length(fX)
        f_nmN = fX(idx).name;
        if contains(f_nmN,f_nm)
            error('The file already exists.')
        end
        
        if whBlk-1 && contains(f_nmN,prevFNm)
            isVld = true;
        elseif ~(whBlk-1) && contains(f_nmN,init)
            error('A file with the same initials is detected.')
        elseif ~(whBlk-1) && ~contains(f_nmN,init)
            isVld = true;
        end 
    end
end

if ~isVld
    error('The file of the previous trial could not be found.')
end
end