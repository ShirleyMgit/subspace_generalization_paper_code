function proj_auc_gp = cal_proj_auc_random_gp_dist2sides(env1_data, env2_data, num_permutation)

% envX_data: features X cells/voxels, gird and place cells together
% num_cells should be even

n_features = length(env1_data(:, 1));
num_cells = length(env1_data(1, :));

% env1
env1_data = env1_data - repmat(mean(env1_data),n_features,1); % substract the mean over condition

% env2
env2_data = env2_data - repmat(mean(env2_data),n_features,1); % substract the mean over condition

proj_auc_gp = -ones(num_permutation, 1);
parfor (n = 1:num_permutation, 3)
%for n = 1:num_permutation
    % select cells and calculate PCs
    v_rand = randperm(num_cells);
    d_env1 = env1_data(:, v_rand);
    d_env1 = d_env1(:, num_cells/2)
    [U1,S1,~] = svd(d_env1'*d_env1, 'econ');
    d_env2 = env2_data(:, v_rand);
    d_env2 = d_env2(:, num_cells/2)
    [U2,S2,~] = svd(d_env2'*d_env2,'econ');
    % project:
    proj21 = (d_env2*U1).^2;
    proj12 = (d_env1*U2).^2;
    sum_var21 = sum(proj21)/trace(S2);
    cumsum_var21 = cumsum(sum_var21,2);
    sum_var12 = sum(proj12)/trace(S1);
    cumsum_var12 = cumsum(sum_var12,2);
    cumsum_var = 0.5*(cumsum_var12 + cumsum_var21)
    proj_auc_gp(n) = sum(cumsum_var);

end
