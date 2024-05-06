clear all
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
more_code = fullfile(exp_root, 'cells_data_muller');
data_folder = fullfile(more_code, 'wetransfer_neil-s-data_2023-09-25_1057');
addpath(genpath(more_code));
addpath(genpath(data_folder));

grid_rats = [1, 3];
num_grid_rats = length(grid_rats);

for n=1:num_grid_rats
    load(['grid_animal',num2str(grid_rats(n)),'.mat'])
    num_grids = length(cells_array(1, :, 1));
    env1_data = squeeze(cells_array(1, :, :));
    env1_data(isnan(env1_data)) = 0;
    env2_data = squeeze(cells_array(2, :, :));
    env2_data(isnan(env2_data)) = 0;

    cumsum_var_grid11 = cal_projection_plot(transpose(env1_data), transpose(env1_data));
    cumsum_var_grid12 = cal_projection_plot(transpose(env1_data), transpose(env2_data));

    cumsum_var_grid22 = cal_projection_plot(transpose(env2_data), transpose(env2_data));
    cumsum_var_grid21 = cal_projection_plot(transpose(env1_data), transpose(env1_data));

    cumsum_var_grid_same = 0.5 * (cumsum_var_grid11 + cumsum_var_grid22)/ num_grids;
    cumsum_var_grid_dif = 0.5 * (cumsum_var_grid12 + cumsum_var_grid21)/ num_grids;

    num_grid_cells = length(cumsum_var_grid_same);
    pc_fraction_grid = linspace(0,1,length(cumsum_var_grid_same)+1);
    figure(1)
    subplot(2,1,n)
    plot(pc_fraction_grid, [0, cumsum_var_grid_same], 'k', 'LineWidth', 2)
    hold on
    plot(pc_fraction_grid, [0, cumsum_var_grid_dif], 'k:', 'LineWidth', 2)

    num_permutation = 5000;
    proj_auc_grid = cal_proj_auc_random_dist_2sides(transpose(env1_data), transpose(env2_data), num_permutation);

    % calculate stats:
    grid_effect = (sum(cumsum_var_grid_same) - sum(cumsum_var_grid_dif));
    null_dist = sum(cumsum_var_grid_same) - proj_auc_grid;
    [count_null_grid, auc_dif_nul] = hist(null_dist, 100);
    p_null = count_null_grid / num_permutation;
    cdf_null = cumsum(p_null);

    n_effect = sum(auc_dif_nul < grid_effect);
    if n_effect > 0
        p_val_grid = cdf_null(n_effect);
    else
        p_val_grid = 0;
    end

    figure(10)
    subplot(2,1,n)
    plot(auc_dif_nul, cdf_null)
    title(['grid permutation cdf, p = ', num2str(p_val_grid)])

end
