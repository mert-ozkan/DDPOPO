function blkX = setTrialOrder_inDDPOPO
dispX = [...
    112122,...
    122122,...
    211112,...
    221211];
dispX = combvec(0:4,dispX);
dispX = dispX(2,:)*10 + dispX(1,:);

qRep_perDisp = 2;
rng_qTrl_inBlk = [30 35];
rep_cond = combvec(dispX,dispX);
test_trial_randomization(rep_cond)

trl_sq = generate_randomly_ordered_trials(rep_cond,qRep_perDisp);
blkX = divide_trials_inblocks(trl_sq,rng_qTrl_inBlk);
blkX = insert_random_sequences(blkX, trl_sq, 7, 27);
blkX = assign_position(blkX);
q_rep = number_of_repetitions(blkX);
end
%%
function q_rep = number_of_repetitions(blkX)

[qBlk, qTrl] = size(blkX);
digX = zeros(8, qBlk, qTrl);

for whDig = 1:8
    for whBlk = 1:qBlk
        for whTrl = 1:qTrl
            digX(whDig,whBlk,whTrl) = return_digit(blkX(whBlk,whTrl),whDig);
        end
    end
end
for whDig = 1:8
    addxn = 0;
    for whBlk = 1:qBlk
        histo = count_repeating_elements(squeeze(digX(whDig,whBlk,:))');
        addxn = addxn + sum(histo(:,2));
    end
    q_rep(whDig,1) = addxn;
end
end
%%
function histo_rep = count_repeating_elements(vec)

qDig = length(num2str(vec(1)));
vec2 = [vec(2:end),0];
rep_id = unique(vec)+unique(vec)*10^qDig;
vec = vec*10^qDig+vec2;
vec(end) = [];
[q,id] = hist(vec,unique(vec));

idx_rep_id = ismember(id, rep_id);

histo_rep = [id(idx_rep_id)', q(idx_rep_id)'];
end
%%
function dig = return_digit(num,whDig)
st = num2str(num);
dig = str2double(st(whDig));
end
%%
function valid = test_random_positions(blkX)

valid = true;

[qBlk, qTrl] = size(blkX);
posX = zeros(qBlk, qTrl);
for whBlk = 1:qBlk
    for whTrl = 1:qTrl
        currPos = num2str(blkX(whBlk,whTrl));
        currPos = str2double(currPos(end));
        posX(whBlk,whTrl) = currPos;
    end
end

[q, ~] = hist(posX,unique(posX));
q = sum(q,2);
q = q/(qBlk*qTrl);
if sum(q<.11) || sum(q>.14)
    valid = false;
end

end
%%
function blkX = assign_position(blkX)
org_blkX = blkX;
valid = false;
while ~valid
    
    blkX = org_blkX;
    
    [qBlk, qTrl] = size(blkX);
    blkX(:,1) = blkX(:,1)*10 + floor(unifrnd(1,9,qBlk,1));

    for whBlk = 1:qBlk
        for whTrl = 2:qTrl
            prevPos = num2str(blkX(whBlk,whTrl-1));
            prevPos = str2double(prevPos(end));

            dist = num2str(blkX(whBlk,whTrl));
            dist = str2double(dist(end));

            sign = mod(randi(1000),2)*-2+1;

            currPos = mod(prevPos+sign*dist-1,8)+1;
            blkX(whBlk,whTrl) = blkX(whBlk,whTrl)*10+currPos;
        end
    end
    
    valid = test_random_positions(blkX);
    
end

end
%%
function op = divide_trials_inblocks(trl_sq,rng_qTrl_inBlk, recursive, test)

if nargin < 3
    recursive = false;  
end

if nargin < 4
    test = false;
end

if recursive
    qTrl_inBlk = rng_qTrl_inBlk;
else
    qTrl_inBlk = find_optimal_qTrl_inBlk(trl_sq,rng_qTrl_inBlk);
end

qTrl = length(trl_sq);
idx = qTrl_inBlk+1;
while idx <= qTrl
    if test
        trl_atIdx = trl_sq(idx-1);
    else
        trl_atIdx = floor(trl_sq(idx-1)/10)*10;
    end
    trl_sq = insert(trl_sq, trl_atIdx, idx);
    idx = idx + qTrl_inBlk;
    qTrl = length(trl_sq);
end

if recursive
    op = trl_sq;
else
    blkX = [];
    for idx = 1:qTrl_inBlk:length(trl_sq)
        blkX = [blkX; trl_sq((1:qTrl_inBlk)+idx-1)];
    end
    
    blkX = blkX(randperm(length(trl_sq)/qTrl_inBlk),:);
    
    op = blkX;
end
end

%%
function best_qTrl =find_optimal_qTrl_inBlk(trl_sq,rng_qTrl)

org_trl_sq = trl_sq;
qTrlX = [];
for whQTrl = min(rng_qTrl):max(rng_qTrl)
    trl_sq = divide_trials_inblocks(org_trl_sq,whQTrl,true);
    if mod(length(trl_sq),whQTrl) == 0
        qTrlX = [qTrlX, whQTrl];
    end
end

if isempty(qTrlX)
    error('Equal number of trials in each block is not possible in this range.')
else
    best_qTrl = qTrlX(ceil(length(qTrlX)/2));
end

end
%%
function trl_sq = generate_randomly_ordered_trials(rep_cond,qRep_perDisp)

lng = 0;
while lng ~= length(rep_cond)*qRep_perDisp+1    
    trl_sq =[];
    rep_cond(3,:) = ones(length(rep_cond),1)*qRep_perDisp;
    while sum(rep_cond(3,:)) ~= 0
        
        prevTrl = length(trl_sq);
        currTrl = prevTrl + 1;
        
        if currTrl == 1
            pick_cond = randi(length(rep_cond));
            digitF = num2str(rep_cond(1,pick_cond));
            digitF = str2double(digitF(end));
            while digitF ~= 0
                pick_cond = randi(length(rep_cond));
                digitF = num2str(rep_cond(1,pick_cond));
                digitF = str2double(digitF(end));
            end
            trl_sq(currTrl+[0,1]) = [rep_cond(1,pick_cond), rep_cond(2,pick_cond)];
        else
            lst_cond = find(trl_sq(prevTrl)==rep_cond(1,:));
            if sum(rep_cond(3,lst_cond)) == 0
                break
            else
                pick_cond = lst_cond(randi(length(lst_cond)));
                while rep_cond(3,pick_cond) == 0
                    pick_cond = lst_cond(randi(length(lst_cond)));
                end
                trl_sq(currTrl) = rep_cond(2,pick_cond);
            end
        end
        rep_cond(3,pick_cond) = rep_cond(3,pick_cond)-1;
    end
    lng = length(trl_sq);
end

end

%%
function vec = insert(vec, el, idx)
if idx > length(vec)
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
function test_trial_randomization(rep_cond)

trl_sq = generate_randomly_ordered_trials(rep_cond,3);
blkX = divide_trials_inblocks(trl_sq,[30 35],false,true);

[qBlk, ~] = size(blkX);
qDig = length(num2str(blkX(1)));
repX = [];
dispX = [];
for whBlk = 1:qBlk
    currBlk = blkX(whBlk,:);
    
    rep = currBlk*10^qDig + [currBlk(2:end), 0];
    rep(end) = [];
    
    [qRep_inBlk, dispX_inBlk] = hist(rep, unique(rep));
    
    for idx = 1:length(dispX_inBlk)
        whDisp = dispX_inBlk(idx);
        whRep = qRep_inBlk(idx);
        
        idx_disp = ismember(dispX, whDisp);
        if sum(idx_disp) == 1
            repX(idx_disp) = repX(idx_disp) + whRep;
        else
            dispX = [dispX, whDisp];
            repX = [repX, whRep];
        end
    end
end

dispX(2,:) = dispX-floor(dispX(1,:)/10^qDig)*10^qDig;
dispX(1,:) = floor(dispX(1,:)/10^qDig);

if ~isequal(sortrows(dispX')',sortrows(rep_cond')')
    error('Test failed');
else
    disp('Test passed.');
end

if sum(repX~=3)
    error('Test failed.')
else
    disp('Test passed.');
end
end
%%
function new_blkX = insert_random_sequences(blkX, trl_sq, q_sq, q_el_inSq)
[q_blk, q_trl_inBlk] = size(blkX);
for whBlk = 1:q_blk
    blkN = blkX(whBlk,:);
    
    [q_trl,~] = randfixedsum(q_sq,1,q_el_inSq-q_sq,1,q_trl_inBlk);
    q_trl = round(q_trl);
    while sum(q_trl) ~= q_el_inSq-q_sq
        [q_trl,~] = randfixedsum(q_sq,1,q_el_inSq-q_sq,1,q_trl_inBlk);
        q_trl = round(q_trl);
    end
    
    ord_trl = sort(ceil(rand(1,q_sq)*(q_trl_inBlk-2)))+1;
    for idx = q_sq:-1:1 %descending order
        rnd_sqN = [trl_sq(ceil(rand(1,q_trl(idx))*length(trl_sq))), blkN(ord_trl(idx)-1)];
        blkN = insert(blkN, rnd_sqN, ord_trl(idx));
    end
    
    new_blkX(whBlk,:) = blkN;
end
end