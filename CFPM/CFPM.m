%% CFPM pipeline

% clean up workspace
clearvars; close all; clc

% add helpful functions
addpath(genpath('func'));

%% Generate example edge time series and behavioral data
% example edge time series dim: #time * #edges * #subject
% example behavioral data dim: 1*1 struct with fields. One field should be
% the framewise displacement; the others should be the behavioral variables
% interested. Dim of each field should be #time * #subject * 1.
time                      = 480;
subj_num                  = 32;
max_node                  = 268;
feat_num                  = ((max_node * max_node) - max_node)/2;
ets_FD                    = rand(time * subj_num, 1);
behav_var_y               = rand(time * subj_num, 1);
subjNum                   = repelem([1:subj_num]', time);

example_behav.subjNum     = subjNum;
example_behav.behav_var_y = behav_var_y;
example_behav.ets_FD      = ets_FD;
example_ets               = rand(time, feat_num, subj_num);
%% Load example edge time series and behavioral data
%load example_ets example_behav

%% Calculate observed brain-behavior correlation
[r_matrix, p_matrix,~] = gen_r_matrix(example_behav, example_ets);
r_observed_fisherz = atanh(r_matrix);

% save true_r_matrix
% dlmwrite('./output/r_observed_fisherz.txt', r_observed_fisherz, 'delimiter','\t');

%% Internal training (cross validation for edge feature selection)
[r_pos,    ... % get internal cross-validation Rho values for positive edges
 r_neg,    ... % get internal cross-validation Rho values for negative edges
 mask_pos, ... % get the positive CFPM mask
 mask_neg, ... % get the negative CFPM mask
 corr_data_pos, ...
 behav_var_sub ...
] = internal_train(example_behav, example_ets, r_observed_fisherz, subj_num, 0);

% Save internal CV Rho values.
% dlmwrite('./output/masks/r_CV_pos.txt', r_pos', 'delimiter','\t');
% dlmwrite('./output/masks/r_CV_neg.txt', r_neg', 'delimiter','\t');

% Save identified masks.
% save('./output/masks/neg_mask.mat', 'neg_mask')
% save('./output/masks/pos_mask.mat', 'pos_mask')

%% External testing

% Read in edge time series data from external dataset
time_external        = 340;
subj_num_external    = 20;
example_ets_external = rand(time_external, feat_num, subj_num_external, 1);

%%
% Load the above masks from internal CV
% load('./output/masks/pos_mask.mat');
% load('./output/masks/neg_mask.mat');

At = mask_pos.';
m  = triu(true(size(At)),1);
pos_flat  = At(m).';

At = mask_neg.';
m  = triu(true(size(At)),1);
neg_flat  = At(m).';
%%
% Apply and obtain time resolved edge score in the CFPM [avg(high) - avg(low)]
[CFPM_score, CFPM_score_pos, CFPM_score_neg] = get_CFPM_score(example_ets_external, subj_num_external, pos_flat, neg_flat);
% writetable(array2table(CFPM_score),'./output/CFPM_score.csv');

% Calculate simple correlation between brain CFPM score and behavior or
% build regression models using brain CFPM score to predict behavior



