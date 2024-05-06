clear all
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
more_code = fullfile(exp_root, 'cells_data_muller');
data_folder = fullfile(more_code, 'wetransfer_neil-s-data_2023-09-25_1057');
addpath(genpath(more_code));
addpath(genpath(data_folder));

place_cells_rats = [4, 6, 7, 8];
num_place_rats = length(place_cells_rats);

for n=1:num_place_rats
    load(['place_cells_animal',num2str(place_cells_rats(n)),'.mat'])

    env1_data = squeeze(cells_array(1, :, :));
    env1_data(isnan(env1_data)) = 0;
    env2_data = squeeze(cells_array(2, :, :));
    env2_data(isnan(env2_data)) = 0;

    cumsum_var_place11 = cal_projection_plot(transpose(env1_data), transpose(env1_data));
    cumsum_var_place12 = cal_projection_plot(transpose(env1_data), transpose(env2_data));
    cumsum_var_place21 = cal_projection_plot(transpose(env2_data), transpose(env1_data));
    cumsum_var_place22 = cal_projection_plot(transpose(env2_data), transpose(env2_data));
    
    num_place_cells = length(cumsum_var_place11);
    cumsum_var_place_same = 0.5 * (cumsum_var_place11 + cumsum_var_place22) / num_place_cells;
    cumsum_var_place_dif = 0.5 * (cumsum_var_place12 + cumsum_var_place21) / num_place_cells;

    disp(num_place_cells)
    pc_fraction_place = linspace(0, 1, num_place_cells+1);
    figure(1)
    subplot(num_place_rats,1,n)
    plot(pc_fraction_place, [0, cumsum_var_place11], 'g', 'LineWidth', 2)
    hold on
    plot(pc_fraction_place, [0, cumsum_var_place12], 'g:', 'LineWidth', 2)

    num_permutation = 8000;
    proj_auc_place = cal_proj_auc_random_dist_2sides(transpose(env1_data), transpose(env2_data), num_permutation);

    % calculate stats:
    place_effect = sum(cumsum_var_place_same) - sum(cumsum_var_place_dif);
    null_dist = sum(cumsum_var_place_same) - proj_auc_place;
    [count_null_place, auc_dif_nul] = hist(null_dist,100);
    p_null = count_null_place / num_permutation;
    cdf_null = cumsum(p_null);

    n_effect = sum(auc_dif_nul < place_effect);
    if n_effect > 0
        p_val_place = cdf_null(n_effect);
    else
        p_val_place = 0;
    end

    figure(10)
    subplot(num_place_rats,1,n)
    plot(auc_dif_nul, cdf_null)
    title(['place cells permutation cdf, p = ', num2str(p_val_place)])

end
