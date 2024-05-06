function proj_auc = cal_proj_auc_random_dist_2sides(env1_data, env2_data, num_permutation)

% envX_data: features X cells/voxels

n_features = length(env1_data(:, 1));
num_cells = length(env1_data(1, :));

% env1
env1_data = env1_data - repmat(mean(env1_data),n_features,1); % substract the mean over condition
[U1,S1,~] = svd(env1_data'*env1_data,'econ');
% env2
env2_data = env2_data - repmat(mean(env2_data),n_features,1); % substract the mean over condition
[U2,S2,~] = svd(env2_data'*env2_data,'econ');

proj_auc = -ones(num_permutation, 1);
parfor (n = 1:num_permutation, 3)
%for n = 1:num_permutation
    % permute and project:
    proj21 = (env2_data(:, randperm(num_cells))*U1).^2;
    proj12 = (env1_data(:, randperm(num_cells))*U2).^2;
    sum_var21 = sum(proj21)/trace(S2);
    sum_var12 = sum(proj12)/trace(S1);
    sum_var = 0.5 * (sum_var21 + sum_var12);
    cumsum_var = cumsum(sum_var,2)/num_cells;

    proj_auc(n) = sum(cumsum_var);

end
