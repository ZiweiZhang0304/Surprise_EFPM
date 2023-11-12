function [r_matrix, p_matrix,behav_var_sub,ets_FD_sub] = gen_r_matrix(behav_dat, brain_dat)

% Obtain trial-to-trial difference in human & model behavior.

% INPUT
%          brain_dat: brain edge time series data. Dim: #time * #edges * #subject.
%          behav_dat: 

% OUTPUT
%          r_matrix: 
%          p_matrix:
%%
behav_var = behav_dat.behav_var_y;
ets_FD = behav_dat.ets_FD;

time           = 480;
max_node       = 268;
lower_tri_size = (max_node * max_node - max_node)/2;

r_matrix       = NaN(lower_tri_size,size(unique(behav_dat.('subjNum')),1));
p_matrix       = NaN(lower_tri_size,size(unique(behav_dat.('subjNum')),1));

for s = 1: max(behav_dat.subjNum)   
    
    % find the x variable for this subject and run
    sub_behav = behav_dat.subjNum==s;
    behav_var_sub = behav_var(sub_behav,:);
    ets_FD_sub = ets_FD( (1+time*(s-1)):(time+time*(s-1)), :);
    disp(size(ets_FD_sub));
   
    sub_ets = brain_dat(:,:,s);

    corr_data = [sub_ets';behav_var_sub';ets_FD_sub']'; % size should be 480 * (35778 + 2)

    % correlate trial-to-trial model parameter with edge time series
    for edge = 1:size(sub_ets,2)
        corr_data_edge = corr_data(:, [edge, end-1, end]);
        corr_data_edge_dropna = corr_data_edge(sum(isnan(corr_data_edge),2)==0,:);
        
        disp("this is the size of correlation data drop na")
        disp(size(corr_data_edge_dropna));
        
        if size(corr_data_edge_dropna,1) == 0
            r_matrix(edge, s) = NaN; 
            p_matrix(edge, s) = NaN;
        else
            sub_ets_trim_edge_dropna = corr_data_edge_dropna(:,1);
            behav_var_sub_trim_edge_dropna = corr_data_edge_dropna(:,end-1);
            ets_FD_sub_trim_edge_dropna = corr_data_edge_dropna(:,end); %size(corr_data_edge_dropna,2)
            
            [r_val, p_val] = partialcorr(sub_ets_trim_edge_dropna, behav_var_sub_trim_edge_dropna, ets_FD_sub_trim_edge_dropna, 'Type','Spearman');%,'rows','complete' somehow doesn't work with a covariate vector with NaN

            p_matrix(edge, s) = p_val;
            r_matrix(edge, s) = r_val;    
        end
    end
end

