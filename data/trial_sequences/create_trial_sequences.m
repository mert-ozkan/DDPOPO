%%
addpath(genpath('/Users/mertozkan/Documents/MATLAB/DD/DDPOPO'))
addpath('/Users/mertozkan/Documents/MATLAB/PTB')
dr = '/Users/mertozkan/Documents/MATLAB/DD/DDPOPO/data/trial_sequences';
for no = 1:20
    print_trial_order(dr,no)
end
%%
[a,b] = check_sq_similarity(20);
%%
function f = print_trial_order(dr,no)
if nargin < 2
    no = input('Initials: ','s');
end
no = string(no);
f_nm = sprintf('Trial_Sequence_No%s.txt',no);
if strcmp(no,'try') || strcmp(no, 'test')
    f_nm = 'test.txt';
end
f.directory = dr;
f.uniquefilename = false;
f = write_inFile(f,'open',f_nm);

blkX = setTrialOrder_inDDPOPO;
write_inFile(f,'matrix',blkX);
end
%%
function [comp_mat,similarity_count] = check_sq_similarity(qSq)
for whSq = 1:qSq
    f_nm = sprintf('Trial_Sequence_No%d.txt',whSq);
    blkX = imp_blkX(f_nm);
    
    [qBlk, qTrl] = size(blkX);
    for whBlk = 1:qBlk
        st = num2str(blkX(whBlk,:));
        st(st==' ') = '';
        comp_mat(whBlk,whSq) = string(st);        
    end
end
qDig = length(num2str(blkX(1)));
similarity_count = zeros(1,qSq);
seg_lng = 5;

for whSq = 1:qSq
    for whBlk = 1:qBlk
        blkN = char(comp_mat(whBlk,whSq));
        for whSeg = 1:qDig*seg_lng:qTrl*qDig-qDig*seg_lng
            segN = blkN(whSeg:whSeg+qDig*seg_lng);
            for whCompSq = 1:qSq
                for whCompBlk = 1:qBlk
                    if whCompSq ~= whSq
                        comp_blk = char(comp_mat(whCompBlk,whCompSq));
                        for whCompSeg = 1:qDig*seg_lng:qTrl*qDig-qDig*seg_lng
                            comp_seg = comp_blk(whCompSeg:whCompSeg+qDig*seg_lng);
                            similarity_count(whSq) = similarity_count(whSq) + strcmp(segN,comp_seg);
                        end
                        
                    end
                end
            end
        end
    end
end
        
    

end