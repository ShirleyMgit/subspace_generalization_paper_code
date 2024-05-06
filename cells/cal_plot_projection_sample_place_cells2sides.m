clear all
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
more_code = fullfile(exp_root, 'cells_data_muller');
data_folder = fullfile(more_code, 'wetransfer_neil-s-data_2023-09-25_1057');
addpath(genpath(more_code));
addpath(genpath(data_folder));


load('grid_animal4.mat')

env1_data = squeeze(cells_array(1, :, :));
env1_data(isnan(env1_data)) = 0;
env2_data = squeeze(cells_array(2, :, :));
env2_data(isnan(env2_data)) = 0;

cumsum_var_grid11 = cal_projection_plot(transpose(env1_data), transpose(env1_data));
cumsum_var_grid12 = cal_projection_plot(transpose(env1_data), transpose(env2_data));
cumsum_var_grid22 = cal_projection_plot(transpose(env2_data), transpose(env2_data));
cumsum_var_grid21 = cal_projection_plot(transpose(env2_data), transpose(env1_data));

cumsum_var_grid_same = 0.5 * (cumsum_var_grid11 + cumsum_var_grid22);
cumsum_var_grid_dif = 0.5 * (cumsum_var_grid12 + cumsum_var_grid21);

grid_effect12 = sum(cumsum_var_grid11) - sum(cumsum_var_grid12);
grid_effect21 = sum(cumsum_var_grid22) - sum(cumsum_var_grid21);
grid_effect = 0.5*(grid_effect12 + grid_effect21);

num_grids = length(cumsum_var_grid11);
pc_fraction_grid = linspace(0,1,length(cumsum_var_grid11)+1);
figure(1)
plot(pc_fraction_grid, [0, cumsum_var_grid_same], 'k', 'LineWidth', 2)
hold on
plot(pc_fraction_grid, [0, cumsum_var_grid_dif], 'k:', 'LineWidth', 2)


load('place_cells_animal4.mat')
num_selections = 5000;
num_place_cells = length(cells_array(1, : , 1));

env1_data_all = squeeze(cells_array(1, :, :));
env1_data_all(isnan(env1_data_all)) = 0;
env2_data_all = squeeze(cells_array(2, :, :));
env2_data_all(isnan(env2_data_all)) = 0;

cumsum_var_p11_samples = -ones(num_selections, num_grids);
cumsum_var_p12_samples = -ones(num_selections, num_grids);
cumsum_var_p22_samples = -ones(num_selections, num_grids);
cumsum_var_p21_samples = -ones(num_selections, num_grids);

cumsum_var_same_samples = -ones(num_selections, num_grids);
cumsum_var_dif_samples = -ones(num_selections, num_grids);
for n = 1:num_selections
        v_rand = randperm(num_place_cells);
        env1_data = env1_data_all(v_rand(1:num_grids),:);
        env2_data = env2_data_all(v_rand(1:num_grids),:);

        cumsum_var_p11 = cal_projection_plot(transpose(env1_data), transpose(env1_data));
        cumsum_var_p12 = cal_projection_plot(transpose(env1_data), transpose(env2_data));
        cumsum_var_p11_samples(n, :) = cumsum_var_p11;
        cumsum_var_p12_samples(n, :) = cumsum_var_p12;

        cumsum_var_p22 = cal_projection_plot(transpose(env2_data), transpose(env2_data));
        cumsum_var_p21 = cal_projection_plot(transpose(env2_data), transpose(env1_data));
        cumsum_var_p22_samples(n, :) = cumsum_var_p22;
        cumsum_var_p21_samples(n, :) = cumsum_var_p21;

        cumsum_var_same_samples(n, :) = 0.5 * (cumsum_var_p11 + cumsum_var_p22);
        cumsum_var_dif_samples(n, :) = 0.5 * (cumsum_var_p12 + cumsum_var_p21);

end


cumsum_var_same_mean = mean(cumsum_var_same_samples);
cumsum_var_dif_mean = mean(cumsum_var_dif_samples);
cumsum_var_same_q025 = [0, quantile(cumsum_var_same_samples, 0.025)];
cumsum_var_dif_q025 = [0, quantile(cumsum_var_dif_samples, 0.025)];

cumsum_var_same_q975 = [0, quantile(cumsum_var_same_samples, 0.975)];
cumsum_var_dif_q975 = [0, quantile(cumsum_var_dif_samples, 0.975)];

pc_fraction_place = linspace(0,1,length(cumsum_var_same_mean)+1);
figure(1)
hold on
plot(pc_fraction_place, [0, cumsum_var_same_mean], 'g', 'LineWidth', 2)
hold on
plot(pc_fraction_place, [0, cumsum_var_dif_mean], 'g:', 'LineWidth', 2)
legend('grid within-env', 'grid across-env','place cells within-env', 'place cells across-env')
% hold on
% patch([pc_fraction_place, fliplr(pc_fraction_place)], [cumsum_var_p11_q025 fliplr(cumsum_var_p11_q975)],'g','facealpha',0.2,'edgecolor','w')
% hold on
% patch([pc_fraction_place, fliplr(pc_fraction_place)], [cumsum_var_p12_q025 fliplr(cumsum_var_p12_q975)],'g','facealpha',0.2,'edgecolor','w')

% calculate statistics:
dif_auc_place_cells_sample = sum(cumsum_var_same_samples, 2) - sum(cumsum_var_dif_samples, 2);
[count_place, auc_dif_place] = hist(dif_auc_place_cells_sample(:), 200);
p_place_cells = count_place / num_selections;
cdf_place_cells = cumsum(p_place_cells);

figure(10)
plot(auc_dif_place,cdf_place_cells, 'k')
title('cdf place cells dif(suc), p_{val}<0.001')
n_effect = sum(auc_dif_place < grid_effect);
p_val_grid = cdf_place_cells(n_effect);
