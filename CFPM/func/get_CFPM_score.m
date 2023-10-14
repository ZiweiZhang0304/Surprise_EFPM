function [CFPM_score, CFPM_score_pos, CFPM_score_neg] = get_CFPM_score(brain_dat, subj_num_external, pos_flat, neg_flat)
% Apply CFPM and obtain time resolved edge score in the CFPM.

% Inputs:
%    brain_dat: edge time series
%    subs: number of subjects
%    pos_flat: flat vector of the positive CFPM
%    neg_flat: flat vector of the negative CFPM

% Outputs:
%    CFPM_score: trial to trial CFPM score
%    CFPM_score_pos: trial to trial CFPM score, positive network only
%    CFPM_score_neg: trial to trial CFPM score, negative network only

%keep all score
CFPM_score = [];
CFPM_score_pos = [];
CFPM_score_neg = [];

for s = 1:subj_num_external 
    %which subject is this?
    disp(s)
    
    %calculate CFPM score
    sub_ets =[];
        
    temp_ets = brain_dat(:,:,s);
    sub_ets = cat(1, sub_ets, temp_ets);
    
    ets_pos = sub_ets .* pos_flat;
    ets_neg = sub_ets .* neg_flat;

    ets_neg_TR =  mean(ets_neg(:,neg_flat==1),2,'omitnan');
    ets_pos_TR =  mean(ets_pos(:,pos_flat==1),2,'omitnan');
    
    ets_TR = ets_pos_TR - ets_neg_TR;
 
    CFPM_score = [CFPM_score ets_TR];
    CFPM_score_pos = [CFPM_score_pos ets_pos_TR];
    CFPM_score_neg = [CFPM_score_neg ets_neg_TR];
    disp(size(CFPM_score))

end

% CFPM_score= reshape(CFPM_score,[size(CFPM_score,1)*size(subs,2),1]);
% disp(size(CFPM_score))


