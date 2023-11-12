function [r_pos, r_neg, pos_overlap, neg_overlap] = internal_train(behav_dat, ets_dat, r_value_observed, kfolds, shuffle)
% Internal training (cross validation for edge feature selection)

% INPUT
%          ets_dat: brain edge time series data. Dim: #time * #edges * #subject.
%          r_value_observed: observed correlation value between trial-to-trial brain network
%          score and trial-to-trial behavioral variable.
%          kfold: number of fold for cross validation. If using leave-one-out,
%          kfold=1.

% OUTPUT
%          r_pos, observed partial Rho values between brain scores from the selected feature
%          and behavioral variable (positive edges).
%          r_neg, observed partial Rho values between brain scores from the selected feature
%          and behavioral variable (negative edges).
%          pos_overlap, selected positive edges. Surviror of every fold of CV.
%          neg_overlap, selected negative edges. Surviror of every fold of CV.

%%
% Split data into training and testing.
time        = size(ets_dat,1);
nsubs       = size(ets_dat,3);
randinds    = randperm(nsubs); % Randomize subj number.
ksample     = floor(nsubs/kfolds);
node        = 268;
behav_var   = behav_dat.behav_var_y;
ets_FD      = behav_dat.ets_FD;

% Output is the partial Rho values between brain and behavior.
r_pos = NaN(ksample, kfolds);
r_neg = NaN(ksample, kfolds);

for leftout = 1:kfolds
    
    % ------------------ SEPARATE TRAIN AND TEST ------------------ %
    %fprintf('%1.0f ',leftout);
    if kfolds == nsubs % leave-one-out
        testinds=randinds(leftout);
        traininds=setdiff(randinds,testinds);
    else
        si=1+((leftout-1)*ksample);
        fi=si+ksample-1;
        
        testinds=randinds(si:fi);
        traininds=setdiff(randinds,testinds);
    end
    
    
    y_train_r = r_value_observed(:,traininds);
    x_test = ets_dat(:,:,testinds); 
    
    
    % --------------------------- TRAIN --------------------------- %
    
    % Obtain group level t-statistics testing correlation Rho against zero.
    stats_train_for_fdr = NaN(size(y_train_r,1), 1);
    p_train_for_fdr = NaN(size(y_train_r,1), 1);
    
    for i = 1:size(y_train_r,1)
        [h,p,ci,stats] = ttest(y_train_r(i,:));
        stats_train_for_fdr(i) = stats.tstat;
        p_train_for_fdr(i) = p;
    end

    
    % Identify postitive or negative edges.
    idx_t_pos = find(stats_train_for_fdr >0);
    idx_t_neg = find(stats_train_for_fdr <0);
    
    % Identify significant edges.
    idx_sig = find(p_train_for_fdr < 0.01);
    
    % Identify the intersection between significant and positive/negative edges.
    idx_pos = intersect(idx_t_pos,idx_sig);
    idx_neg = intersect(idx_t_neg,idx_sig);
    
    
    % Map the index back into the node * node matrix.
    % Positive network
    pos = zeros(268,268);
    triup = zeros(1,size(y_train_r,1));
    triup(idx_pos) = 1; 
    pos(logical(triu(ones(268),1))) = triup;
    param_pos_mask = pos + pos';
    
    % negative network
    neg = zeros(268,268);
    triup = zeros(1,size(y_train_r,1));
    triup(idx_neg) = 1;
    neg(logical(triu(ones(268),1))) = triup;
    param_neg_mask = neg + neg';
    
 
    % Organize obtained masks from training into flat vectors.
    % Positive network
    At = pos_mask.';
    m  = triu(true(size(At)),1);
    pos_flat  = At(m).';
    
    % Negative network
    At = neg_mask.';
    m  = triu(true(size(At)),1);
    neg_flat  = At(m).';
    
    % Store the obtained mask from this round of training.
    pos_mask_all(:,:,leftout) = pos_mask;
    neg_mask_all(:,:,leftout) = neg_mask;
    
    
    % --------------------------- TEST --------------------------- %
    within_r_pos = NaN(ksample,length(testinds));
    within_r_neg = NaN(ksample,length(testinds));
    
    for s = 1:length(testinds)
        s_idx = testinds(s);
        sub_behav = behav_dat.subjNum==s_idx;
        behav_var_sub = behav_var(sub_behav,:);
        ets_FD_sub= ets_FD( (1+time*(s_idx-1)):(time+time*(s_idx-1)), :);

        sub_ets = x_test(:,:,s);

        ets_pos = sub_ets .* pos_flat;
        ets_neg = sub_ets .* neg_flat;
%         disp("this is the size of brain data before calculating network score")
%         disp(size(ets_pos))
        
        ets_pos_trial = mean(ets_pos(:,pos_flat==1),2,'omitnan');
        ets_neg_trial = mean(ets_neg(:,neg_flat==1),2,'omitnan');
%         disp("this is the size of network score")
%         disp(size(ets_pos_trial))
        
        % Concatenate brain, model and motion separately for pos and neg networks.
        corr_data_pos = [ets_pos_trial';behav_var_sub';ets_FD_sub']'; % size should be time * (feature_num + 2)
        corr_data_neg = [ets_neg_trial';behav_var_sub';ets_FD_sub']';
%         disp("this is the size of correlation data")
%         disp(size(corr_data_pos))
        
        % Positive network
        % Keep time points where none of the three vars (brain, behavior, and motion) is NaN.
        corr_data_pos_dropna = corr_data_pos(sum(isnan(corr_data_pos),2)==0,:);
        corr_data_neg_dropna = corr_data_neg(sum(isnan(corr_data_neg),2)==0,:);
        disp("this is the size of correlation data drop na")
        disp(size(corr_data_pos_dropna));
        
        score_outcome_trial_pos_trim_dropna = corr_data_pos_dropna(:,1);
        behav_var_sub_trim_pos_dropna = corr_data_pos_dropna(:,end-1);
        ets_FD_sub_trim_pos_dropna = corr_data_pos_dropna(:,end);
        
        % If generate null with phase randomized brain data, shuffle is set to 1.
        if shuffle == 1
            if rem(size(corr_data_pos_dropna,1), 2) == 0
                %disp(size(score_outcome_trial_pos_trim_dropna));
                score_outcome_trial_pos_trim_dropna(size(corr_data_pos_dropna,1)+1,1) = score_outcome_trial_pos_trim_dropna(size(corr_data_pos_dropna,1),1);
            end
            disp(size(score_outcome_trial_pos_trim_dropna));
            score_outcome_trial_pos_trim_dropna_surr = phaseran(score_outcome_trial_pos_trim_dropna,1);
            score_outcome_trial_pos_trim_dropna_surr = score_outcome_trial_pos_trim_dropna_surr(1:size(corr_data_pos_dropna,1),1);
            %disp('pos edge surr size');
            %disp(size(score_outcome_trial_pos_trim_dropna_surr));
        end
        
        % Negative network
        score_outcome_trial_neg_trim_dropna = corr_data_neg_dropna(:,1);
        behav_var_sub_trim_neg_dropna = corr_data_neg_dropna(:,end-1);
        ets_FD_sub_trim_neg_dropna = corr_data_neg_dropna(:,end);
        
        % If generate null with phase randomized brain data, shuffle is set to 1.
        if shuffle == 1
            if rem(size(corr_data_neg_dropna,1), 2) == 0
                score_outcome_trial_neg_trim_dropna(size(corr_data_neg_dropna,1)+1,1) = score_outcome_trial_neg_trim_dropna(size(corr_data_neg_dropna,1),1);
            end
            disp(size(score_outcome_trial_neg_trim_dropna));
            score_outcome_trial_neg_trim_dropna_surr = phaseran(score_outcome_trial_neg_trim_dropna,1);
            score_outcome_trial_neg_trim_dropna_surr = score_outcome_trial_neg_trim_dropna_surr(1:size(corr_data_neg_dropna,1),1);
        end
        
        
        % Correlate trial-to-trial behavior with edge time series.
        
        % No positive / negative edges remain after dropna.
        if size(corr_data_pos_dropna,1) == 0 && size(corr_data_neg_dropna,1) == 0
            within_r_pos(:, s) = NaN;  
            within_r_neg(:, s) = NaN;

        % No positive edges remain after dropna, only negative edges remained.
        elseif size(corr_data_pos_dropna,1) == 0 && size(corr_data_neg_dropna,1) ~= 0
            within_r_pos(:, s) = NaN;
            if shuffle == 1
                [r_val_neg, ~] = partialcorr(score_outcome_trial_neg_trim_dropna_surr, behav_var_sub_trim_neg_dropna, ets_FD_sub_trim_neg_dropna, 'Type','Spearman');
            else
                [r_val_neg, ~]= partialcorr(score_outcome_trial_neg_trim_dropna, behav_var_sub_trim_neg_dropna, ets_FD_sub_trim_neg_dropna,'Type','Spearman');
                disp(r_val_neg);
            end
            within_r_neg(:, s) = r_val_neg;
           
            
        % No negative edges remain after dropna, only positive edges remained.
        elseif size(corr_data_pos_dropna,1) ~= 0 && size(corr_data_neg_dropna,1) == 0
            within_r_neg(:, s) = NaN;
            if shuffle == 1
                [r_val_pos, ~] = partialcorr(score_outcome_trial_pos_trim_dropna_surr, behav_var_sub_trim_pos_dropna_surr, ets_FD_sub_trim_pos_dropna, 'Type','Spearman');
            else
                [r_val_pos, ~]= partialcorr(score_outcome_trial_pos_trim_dropna, behav_var_sub_trim_pos_dropna, ets_FD_sub_trim_pos_dropna,'Type','Spearman');
                disp(r_val_pos);
            end
            within_r_pos(:, s) = r_val_pos;            
        
        
        % Both positive and negative edges remain after dropna.
        else
            if shuffle == 1
                [r_val_pos, ~] = partialcorr(score_outcome_trial_pos_trim_dropna_surr, behav_var_sub_trim_pos_dropna, ets_FD_sub_trim_pos_dropna, 'Type','Spearman');%,'rows','complete' somehow doesn't work with a covariate vector with NaN
                [r_val_neg, ~] = partialcorr(score_outcome_trial_neg_trim_dropna_surr, behav_var_sub_trim_neg_dropna, ets_FD_sub_trim_neg_dropna, 'Type','Spearman');
            else
                [r_val_pos, ~]= partialcorr(score_outcome_trial_pos_trim_dropna, behav_var_sub_trim_pos_dropna, ets_FD_sub_trim_pos_dropna,'Type','Spearman');
                [r_val_neg, ~]= partialcorr(score_outcome_trial_neg_trim_dropna, behav_var_sub_trim_neg_dropna, ets_FD_sub_trim_neg_dropna,'Type','Spearman');
            end
            within_r_pos(:, s) = r_val_pos;
            within_r_neg(:, s) = r_val_neg;
        end    
    end
    
    % Store the r for each round of test subjects.
    r_pos(:,leftout) = within_r_pos;
    r_neg(:,leftout) = within_r_neg;
end

% Find egdes appearing in every round of leave-one-out cross-validation.
pos_overlap = zeros(node,node);
neg_overlap = zeros(node,node);
pos_overlap(sum(pos_mask_all,3)==nsubs) = 1;
neg_overlap(sum(neg_mask_all,3)==nsubs) = 1;


