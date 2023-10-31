function [EFPM_score, EFPM_score_pos, EFPM_score_neg] = get_EFPM_score(brain_dat, subj_num_external, pos_flat, neg_flat)
% Apply EFPM and obtain time resolved edge score in the EFPM.

% Inputs:
%    brain_dat: edge time series
%    subs: number of subjects
%    pos_flat: flat vector of the positive EFPM
%    neg_flat: flat vector of the negative EFPM

% Outputs:
%    EFPM_score: trial to trial EFPM score
%    EFPM_score_pos: trial to trial EFPM score, positive network only
%    EFPM_score_neg: trial to trial EFPM score, negative network only

%keep all score
EFPM_score = [];
EFPM_score_pos = [];
EFPM_score_neg = [];

for s = 1:subj_num_external 
    %which subject is this?
    disp(s)
    
    %calculate EFPM score
    sub_ets =[];
        
    temp_ets = brain_dat(:,:,s);
    sub_ets = cat(1, sub_ets, temp_ets);
    
    ets_pos = sub_ets .* pos_flat;
    ets_neg = sub_ets .* neg_flat;

    ets_neg_TR =  mean(ets_neg(:,neg_flat==1),2,'omitnan');
    ets_pos_TR =  mean(ets_pos(:,pos_flat==1),2,'omitnan');
    
    ets_TR = ets_pos_TR - ets_neg_TR;
 
    EFPM_score = [EFPM_score ets_TR];
    EFPM_score_pos = [EFPM_score_pos ets_pos_TR];
    EFPM_score_neg = [EFPM_score_neg ets_neg_TR];
    disp(size(EFPM_score))

end

% EFPM_score= reshape(EFPM_score,[size(EFPM_score,1)*size(subs,2),1]);
% disp(size(EFPM_score))


