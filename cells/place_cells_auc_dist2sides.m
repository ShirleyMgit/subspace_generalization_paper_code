clear all
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
more_code = fullfile(exp_root, 'cells_data_muller');
data_folder = fullfile(more_code, 'wetransfer_neil-s-data_2023-09-25_1057');
addpath(genpath(more_code));
addpath(genpath(data_folder));

%load('grid_animal4.mat')
load('grid_animal3.mat')
num_grids = length(cells_array(1, :, 1));
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
grid_effect = (sum(cumsum_var_grid_same) - sum(cumsum_var_grid_dif))/num_grids ;

place_cells_animals = [4, 6 , 8]; %all 3 animals have 25 place cells
num_place_cells = 25;
num_rats = 3;

num_selections = 1000;
auc_place_cells = -ones(3, num_selections);
for rat = 1:num_rats
    % for each animal choose num_grids cells and calculate the difference
    % in auc
    load(['place_cells_animal',num2str(place_cells_animals(rat)),'.mat'])
    env1_data_all = squeeze(cells_array(1, : , :));
    env1_data_all(isnan(env1_data_all)) = 0;
    env2_data_all = squeeze(cells_array(2, : , :));
    env2_data_all(isnan(env2_data_all)) = 0;

    for n = 1:num_selections
        v_rand = randperm(num_place_cells);
        env1_data = env1_data_all(v_rand(1:num_grids),:);
        env2_data = env2_data_all(v_rand(1:num_grids),:);
        cumsum_var_p11 = cal_projection_plot(transpose(env1_data), transpose(env1_data));
        cumsum_var_p12 = cal_projection_plot(transpose(env1_data), transpose(env2_data));
        cumsum_var_p22 = cal_projection_plot(transpose(env2_data), transpose(env2_data));
        cumsum_var_p21 = cal_projection_plot(transpose(env2_data), transpose(env1_data));

        cumsum_var_same = 0.5 * (cumsum_var_p11 + cumsum_var_p22);
        cumsum_var_dif = 0.5 * (cumsum_var_p12 + cumsum_var_p21);
        auc_place_cells(rat, n) = (sum(cumsum_var_same) - sum(cumsum_var_dif))/num_grids;
    end
end

save(['place_cells_auc_diff_num_cell',num2str(num_grids),'n_selctions',num2str(num_selections),'.mat'],'auc_place_cells');

% calculate statistics:
[count_place, auc_dif_place] = hist(auc_place_cells(:),100);
p_place_cells = count_place / (num_selections*num_rats);
cdf_place_cells = cumsum(p_place_cells);
figure(10)
subplot(2,1,1)
plot(auc_dif_place, cdf_place_cells)
subplot(2,1,2)
plot(auc_dif_place, p_place_cells)

n_effect = sum(auc_dif_place < grid_effect);
p_val_grid = cdf_place_cells(n_effect);