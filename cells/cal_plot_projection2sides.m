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

num_grid_cells = length(cumsum_var_grid11);
pc_fraction_grid11 = linspace(0,1,length(cumsum_var_grid11) + 1);
pc_fraction_grid22 = linspace(0,1,length(cumsum_var_grid22) + 1);
figure(1)
plot(pc_fraction_grid11, [0, cumsum_var_grid11], 'k', 'LineWidth', 2)
hold on
plot(pc_fraction_grid11, [0, cumsum_var_grid12], 'k:', 'LineWidth', 2)

figure(2)
plot(pc_fraction_grid22, [0, cumsum_var_grid22], 'k', 'LineWidth', 2)
hold on
plot(pc_fraction_grid22, [0, cumsum_var_grid21], 'k:', 'LineWidth', 2)

figure(3)
plot(pc_fraction_grid22, [0, cumsum_var_grid22], 'k', 'LineWidth', 2)
hold on
plot(pc_fraction_grid22, [0, cumsum_var_grid21], 'k:', 'LineWidth', 2)

num_permutation = 5000;
proj_auc_grid12 = cal_proj_auc_random_dist(transpose(env1_data), transpose(env2_data), num_permutation);
proj_auc_grid21 = cal_proj_auc_random_dist(transpose(env2_data), transpose(env1_data), num_permutation);
%save(['proj_auc_grid_random_n',num2str(num_permutation),'.mat'],'proj_auc_grid');
save(['proj_auc_grid21_random_n',num2str(num_permutation),'.mat'],'proj_auc_grid21');

load('place_cells_animal4.mat')
all_cells = true;
num_place_cells = length(cells_array(1, : , 1));
if all_cells
    env1_data = squeeze(cells_array(1,: , :));
    env1_data(isnan(env1_data)) = 0;
    env2_data = squeeze(cells_array(2, : , :));
    env2_data(isnan(env2_data)) = 0;
else
    v_rand = randperm(num_place_cells);
    env1_data = squeeze(cells_array(1, v_rand(1:num_grid_cells) , :));
    env1_data(isnan(env1_data)) = 0;
    env2_data = squeeze(cells_array(2, v_rand(1:num_grid_cells) , :));
    env2_data(isnan(env2_data)) = 0;
end
cumsum_var_place11 = cal_projection_plot(transpose(env1_data), transpose(env1_data));
cumsum_var_place12 = cal_projection_plot(transpose(env1_data), transpose(env2_data));

cumsum_var_place22 = cal_projection_plot(transpose(env2_data), transpose(env2_data));
cumsum_var_place21 = cal_projection_plot(transpose(env2_data), transpose(env1_data));

proj_auc_place12 = cal_proj_auc_random_dist(transpose(env1_data), transpose(env2_data), num_permutation);
proj_auc_place21 = cal_proj_auc_random_dist(transpose(env2_data), transpose(env1_data), num_permutation);
%save(['proj_auc_place_random_n',num2str(num_permutation),'.mat'],'proj_auc_place');

pc_fraction_place11 = linspace(0,1,length(cumsum_var_place11)+1);
pc_fraction_place22 = linspace(0,1,length(cumsum_var_place22)+1);

figure(1)
hold on
plot(pc_fraction_place11, [0, cumsum_var_place11], 'g', 'LineWidth', 2)
hold on
plot(pc_fraction_place11, [0, cumsum_var_place12], 'g:', 'LineWidth', 2)

figure(2)
hold on
plot(pc_fraction_place22, [0, cumsum_var_place22], 'g', 'LineWidth', 2)
hold on
plot(pc_fraction_place22, [0, cumsum_var_place21], 'g:', 'LineWidth', 2)

% calculate stats:
grid_effect12 = sum(cumsum_var_grid11) - sum(cumsum_var_grid12);
grid_effect21 = sum(cumsum_var_grid22) - sum(cumsum_var_grid21);

null_dist12grid = sum(cumsum_var_grid11) - proj_auc_grid12;
null_dist21grid  = sum(cumsum_var_grid22) - proj_auc_grid21;
[count_null_grid12, auc_dif_nul12grid ] = hist(null_dist12grid , 200);
[count_null_grid21, auc_dif_nul21grid ] = hist(null_dist21grid , 200);

p_null12 = count_null_grid12 / num_permutation;
cdf_null12grid  = cumsum(p_null12);
figure(10)
plot(auc_dif_nul12grid, cdf_null12grid)

p_null21 = count_null_grid21 / num_permutation;
cdf_null21grid  = cumsum(p_null21);
figure(11)
plot(auc_dif_nul21grid , cdf_null21grid )


n_effect12grid  = sum(auc_dif_nul12grid  < grid_effect12);
if n_effect12grid  > 0
    p_val_grid12 = cdf_null(n_effect12grid );
end

n_effect21grid  = sum(auc_dif_nul21grid  < grid_effect21);
if n_effect21grid  > 0
    p_val_grid21 = cdf_null(n_effect21grid);
end

place_effect12 = sum(cumsum_var_place11) - sum(cumsum_var_place12);
place_effect21 = sum(cumsum_var_place22) - sum(cumsum_var_place21);
null_dist_place12 = sum(cumsum_var_place11) - proj_auc_place12;
null_dist_place21 = sum(cumsum_var_place22) - proj_auc_place21;

[count_null_place12, auc_dif_nul_place12] = hist(null_dist_place12, 200);
[count_null_place21, auc_dif_nul_place21] = hist(null_dist_place21, 200);

p_null_place12 = count_null_place12 / num_permutation;
cdf_null_place12 = cumsum(p_null_place12);
figure(12)
plot(auc_dif_nul_place12, cdf_null_place12)

p_null_place21 = count_null_place21 / num_permutation;
cdf_null_place21 = cumsum(p_null_place21);
figure(13)
plot(auc_dif_nul_place21, cdf_null_place21)

n_effect_place12 = sum(auc_dif_nul_place12 < place_effect12);
p_val_place12 = cdf_null_place12(n_effect_place12);
disp(p_val_place12)

n_effect_place21 = sum(auc_dif_nul_place21 < place_effect21);
p_val_place21 = cdf_null_place21(n_effect_place21);
disp(p_val_place21)