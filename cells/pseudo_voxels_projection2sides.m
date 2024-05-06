clear all
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
more_code = fullfile(exp_root, 'cells_data_muller');
data_folder = fullfile(more_code, 'wetransfer_neil-s-data_2023-09-25_1057');
addpath(genpath(more_code));
addpath(genpath(data_folder));

grid_rats = [3, 4, 1];
num_grid_rats = length(grid_rats);

place_cells_rats = [4, 8, 6];
num_place_rats = length(place_cells_rats);

ncells = 7;
psedo_voxels_number = 2;
num_psudo_subjects = 400;
psudo_voxels_grid = zeros(num_psudo_subjects, 2, psedo_voxels_number*num_grid_rats, 4096);

% create voxel like data:
for n=1:num_grid_rats
    load(['grid_animal',num2str(grid_rats(n)),'.mat']);
    mean_cell = mean(cells_array,2);
    std_cell_array = std(cells_array,[],2);
    norm_cells_arr = (cells_array - mean_cell) ./ std_cell_array;
    psudo_voxels1 = cal_psudo_voxels(norm_cells_arr, num_psudo_subjects, psedo_voxels_number, ncells);
    psudo_voxels_grid(:, :, (n - 1)*psedo_voxels_number + 1:n*psedo_voxels_number, :) = psudo_voxels1;
end

psudo_voxels_place = zeros(num_psudo_subjects, 2, psedo_voxels_number*num_place_rats, 4096);
for n=1:num_place_rats
    load(['place_cells_animal',num2str(place_cells_rats(n)),'.mat']);    mean_cell = mean(cells_array,2);
    std_cell_array = std(cells_array,[],2);
    norm_cells_arr = (cells_array - mean_cell) ./ std_cell_array;
    psudo_voxels1 = cal_psudo_voxels(norm_cells_arr, num_psudo_subjects, psedo_voxels_number, ncells);
    psudo_voxels_place(:, :, (n - 1)*psedo_voxels_number + 1:n*psedo_voxels_number, :) = psudo_voxels1;
end
disp('done creating psudo voxels')

grid_effect12 = zeros(num_psudo_subjects, 1);
place_effect12 = zeros(num_psudo_subjects, 1);
grid_effect21 = zeros(num_psudo_subjects, 1);
place_effect21 = zeros(num_psudo_subjects, 1);
for s = 1:num_psudo_subjects

    [place_effect12(s), place_effect21(s)] = cal_proj_effect(squeeze(psudo_voxels_place(s, :, :, :)));
    [grid_effect12(s), grid_effect21(s)] = cal_proj_effect(squeeze(psudo_voxels_grid(s, :, :, :)));

end

place_effect = [place_effect12; place_effect21];
grid_effect = [grid_effect12; grid_effect21];
[h,p] = kstest2(place_effect, grid_effect)

n_smooth = 9;
nbins = 50;
[h_diff_place, dif_place] = hist(place_effect, nbins);
h_diff_place_smooth = smooth(h_diff_place, n_smooth);
d_dif_place = mean(dif_place(2:end) - dif_place(1:end-1));
p_diff_place = h_diff_place_smooth/(d_dif_place * sum(h_diff_place_smooth));

[h_diff_grid, dif_grid] = hist(grid_effect, nbins);
h_diff_grid_smooth = smooth(h_diff_grid, n_smooth);
d_dif_grid = mean(dif_grid(2:end) - dif_grid(1:end-1));
p_diff_grid = h_diff_grid_smooth/(d_dif_grid * sum(h_diff_grid_smooth));

figure(2)
plot(dif_place, p_diff_place, 'g')
hold on
plot(dif_grid, p_diff_grid, 'k')


function psudo_voxels1 = cal_psudo_voxels(cells_array, num_psudo_subjects, psedo_voxels_number, ncells)
    env1_data = squeeze(cells_array(1, :, :));
    env1_data(isnan(env1_data)) = 0;
    env2_data = squeeze(cells_array(2, :, :));
    env2_data(isnan(env2_data)) = 0;
    
    psudo_voxels1 = zeros(num_psudo_subjects, 2, psedo_voxels_number, 4096);
    num_cells = length(env1_data(:, 1));
    for s = 1:num_psudo_subjects
        p = randperm(num_cells);
        %v_num = floor(num_cells / psedo_voxels_number);
        for n = 1:psedo_voxels_number
            env1_data_p = env1_data(p, :);
            env2_data_p = env2_data(p, :);
            psudo_voxels1(s, 1, n, :) = mean(env1_data_p((1+(n-1)*ncells):(n*ncells) , :));
            psudo_voxels1(s, 2, n, :) = mean(env2_data_p((1+(n-1)*ncells):(n*ncells) , :));
        end
    end
end

function [dif_effect12, dif_effect21, dif_effect] = cal_proj_effect(psudo_voxels)

    env1_data = squeeze(psudo_voxels(1, :, :));
    env2_data = squeeze(psudo_voxels(2, :, :));

    num_cells1 = length(env1_data(:, 1));
    cumsum_var_11 = cal_projection_plot(transpose(env1_data), transpose(env1_data));
    cumsum_var_12 = cal_projection_plot(transpose(env1_data), transpose(env2_data));

    num_cells2 = length(env2_data(:, 1));
    cumsum_var_22 = cal_projection_plot(transpose(env2_data), transpose(env2_data));
    cumsum_var_21 = cal_projection_plot(transpose(env2_data), transpose(env1_data));

    % calculate stats:
    dif_effect12 = (sum(cumsum_var_11) - sum(cumsum_var_12))/  num_cells1;
    dif_effect21 = (sum(cumsum_var_22) - sum(cumsum_var_21))/  num_cells2;

    dif_effect = 0.5*(dif_effect12 + dif_effect21);
end
