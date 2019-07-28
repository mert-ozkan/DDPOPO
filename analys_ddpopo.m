%%
clear all
close all
dr = '/Users/mertozkan/Documents/MATLAB/DD/DDPOPO/data/MOz';

%function analys_ddpopo(dr)
op = concatenate_files(dr);

rxn_t = str2num(char(op(:,6)));

trl = label_trials(op(:,3), op(:,2), rxn_t);

[isOk, acc] = accuracies(op,trl);

rt = reaction_times(trl, rxn_t, isOk);


% end
%%
function rt = reaction_times(trl, rxn_t, isOk)
%%
y_lim = [2 3];
%%
isVld = isOk & trl.vld_rt;

rt.ill = rxn_t(trl.rep.ill & isVld);
rt.phys = rxn_t(trl.rep.phys & isVld);
rt.both = rxn_t(trl.rep.both & isVld);
rt.pos = rxn_t(trl.rep.pos & isVld);
rt.none = rxn_t(trl.rep.none & isVld);
    

condX = {'ill', 'phys', 'both', 'none', 'pos'};

for idx = 1:length(condX)
    rt.stats.err(idx,1) = mean(rt.(condX{idx})(:,1));
    rt.stats.sample_size.(condX{idx}) = length(rt.(condX{idx})(:,1));
    rt.stats.err(idx,2) = std(rt.(condX{idx})(:,1))/sqrt(length(rt.(condX{idx})(:,1)));
end


errorbar(1:5,rt.stats.err(:,1), rt.stats.err(:,2));

ax = gca;
ax.XLim = [.5, 5.5];
ax.YLim = y_lim;
ax.XTick = 1:5;
ax.XTickLabel = {'Illusory Feature Only', 'Physical Feature Only', 'Both Features', 'No Feature Repetitions', 'Position'};
xlabel('Repeated Feature');
ylabel('Reaction Times (s)');
title('Priming of Pop-Out');

%% Spatial Proximity
trl.pos.dist(:,2:5) = NaN;

trl.pos.dist(trl.rep.ill & isVld, 2) = rxn_t(trl.rep.ill & isVld);
trl.pos.dist(trl.rep.phys & isVld, 3) = rxn_t(trl.rep.phys & isVld);
trl.pos.dist(trl.rep.both & isVld, 4) = rxn_t(trl.rep.both & isVld);
trl.pos.dist(trl.rep.none & isVld, 5) = rxn_t(trl.rep.none & isVld);

condX = {'ill', 'phys', 'both','none'};
trl.pos.dist = array2table(trl.pos.dist,'VariableNames',{'dist','ill', 'phys', 'both','none'});

for whDist = 1:5
    for whCond = 1:4
        condN = condX{whCond};
        distN = whDist-1;
        
        datN = rmmissing(trl.pos.dist.(condN)(trl.pos.dist.dist==distN));
        rt.stats.sample_size.spat_prox(whDist,whCond) = length(datN);
        dist(whCond,whDist,1) = mean(datN);
        dist(whCond,whDist,2) = std(datN)/sqrt(length(datN));
    end
end

figure
for whCond = 1:4
    
%     errorbar(0:4, dist(whCond,:,1),dist(whCond,:,2));
    hold on
end

bar(0:4,dist(:,:,1)')

ax = gca;
ax.XLim = [-.5, 4.5];
ax.XTick = 0:4;
ax.YLim = y_lim;
legend(condX);
xlabel('Target Distances between Trials by Item');
ylabel('Reaction Times (s)');
title('Spatial Proximity of Succesive Targets');

%% Hemifield Repetitions

hemX = {'v_hem', 'h_hem', 'diff_quad'};
condX = {'ill', 'phys', 'both','none'};

for whCond = 1:length(condX)
    
    for whHem = 1:length(hemX)
        hemN = hemX{whHem};
        condN = condX{whCond};
        datN = rxn_t(trl.pos.loc.(hemN) & trl.rep.(condN) & isVld);
        
        rt_hem(whCond,whHem,1) = mean(datN);
        rt_hem(whCond,whHem,2) = std(datN)/sqrt(length(datN));
        rt.stats.sample_size.hemXcond(whHem,whCond) = length(datN);
    end
end

condX = {'Illusory Feature Only', 'Physical Feature Only', 'Both Features', 'No Feature Repetitions'};
figure
for whCond = 1:length(condX)
    subplot(1,length(condX),whCond)
    errorbar(1:3, rt_hem(whCond,:,1),rt_hem(whCond,:,2));
    
    ax = gca;
    ax.XLim = [.5, 3.5];
    ax.YLim = y_lim;
    ax.XTick = 1:3;
    ax.XTickLabels = {'Top/Bottom','Left/Right','Different Quadrant'};
    title(condX{whCond});
    xlabel('Location of Repetition')
    ylabel('Reaction Times (s)')
end




end
%%
function trl = label_trials(rep, cryp, rxn_t)
rep = char(rep);
cryp = char(cryp);

isRepIll = logical(str2num(rep(:,1)));
isRepPhys = logical(str2num(rep(:,2)));
isRepPos = logical(str2num(rep(:,8)));

isRepBoth = isRepIll & isRepPhys;
isRepIll = isRepIll & ~isRepBoth;
isRepPhys = isRepPhys & ~isRepBoth;

isRepTrl = table(isRepIll, isRepPhys, isRepBoth, isRepPos,'VariableNames',{'ill','phys', 'both','pos'});
isPrevTrl = ones(length(isRepIll),3)&false;
isPrevTrl(find(isRepIll)-1,1) = true;
isPrevTrl(find(isRepPhys)-1,2) = true;
isPrevTrl(find(isRepBoth)-1,3) = true;
isPrevTrl(find(isRepPos)-1,4) = true;
isPrevTrl = array2table(isPrevTrl,'VariableNames',{'ill','phys', 'both', 'pos'});

trl.rep = isRepTrl;
trl.rep.none = ~(trl.rep.ill|trl.rep.phys|trl.rep.both);

trl.prev = isPrevTrl;

dist = str2num(cryp(:,7));
loc = str2num(cryp(:,8));

dic_v_hem = [1, 2, 7, 8;
    3, 4, 5, 6];
dic_h_hem = [1:4; 5:8];

relv_loc = [loc, [loc(2:end); 0]];

for idx = 1:length(relv_loc)
    check_v = ismember(dic_v_hem,relv_loc(idx,:));
    check_h = ismember(dic_h_hem,relv_loc(idx,:));
    
    isRepVHem(idx) = ismember(2,sum(check_v,2));
    isRepHHem(idx) = ismember(2,sum(check_h,2));
end

trl.pos.dist = dist;

trl.pos.loc = loc;

isNoRepHem = ~(isRepVHem | isRepHHem);
trl.pos.loc(:,2) = isRepVHem;
trl.pos.loc(:,3) = isRepHHem;
trl.pos.loc(:,4) = isNoRepHem;

trl.pos.loc = array2table(trl.pos.loc,'VariableNames',{'loc','v_hem','h_hem','diff_quad'});
trl.vld_rt = rxn_t<(mean(rxn_t)+2*std(rxn_t)) & rxn_t>(mean(rxn_t)-2*std(rxn_t));
end
%%
function op = concatenate_files(dr)
f_nmX = struct2table(dir(dr));
f_nmX = f_nmX.name(3:end);

op = [];
for idx = 2:length(f_nmX)
    f_nm = f_nmX{idx};
    opN = imp_ddpopo(f_nm);
    op = [op; opN];
end
end
%%
function op = imp_ddpopo(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   DDPOPOMOZTRLSQNO16BLKNO717JUL2019142222 = IMPORTFILE(FILENAME) Reads
%   data from text file FILENAME for the default selection.
%
%   DDPOPOMOZTRLSQNO16BLKNO717JUL2019142222 = IMPORTFILE(FILENAME,
%   STARTROW, ENDROW) Reads data from rows STARTROW through ENDROW of text
%   file FILENAME.
%
% Example:
%   DDPOPOMOzTrlSqNo16BlkNo717Jul2019142222 = importfile('DDPOPO_MOz_TrlSqNo16_BlkNo7_17Jul2019_142222', 2, 32);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/07/17 16:34:40

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
%	column4: text (%s)
%   column5: text (%s)
%	column6: text (%s)
%   column7: text (%s)
%	column8: text (%s)
%   column9: text (%s)
%	column10: text (%s)
%   column11: text (%s)
%	column12: text (%s)
%   column13: text (%s)
%	column14: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
op = [dataArray{1:end-1}];

end
%%
function [isOk, acc] = accuracies(op,trl)
isOk = strcmp(op(:,4),op(:,5));

pc_ok = sum(isOk)*100/length(isOk);

acc = [sum(isOk(trl.rep.ill))/sum(trl.rep.ill);
    sum(isOk(trl.rep.phys))/sum(trl.rep.phys);
    sum(isOk(trl.rep.both))/sum(trl.rep.both);
    sum(isOk(trl.rep.pos))/sum(trl.rep.pos);
    sum(isOk(trl.rep.none))/sum(trl.rep.none)]*100;
    
end
