clear all
close all
clc

% Calculates subspace generalization between simulated voxels.
% Each voxel is an average over simulated grid cells. Cells are sampled
% randomly. This script examine the signal as a function of the numer of
% voxels per module.
% noise is added to the voxels' activity map.

global cells_in_voxels

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
more_code = fullfile(exp_root, 'cells_data_muller');

saving_folder = 'C:\Users\User\Documents\fMRI_EXP\cells_data_muller\wetransfer_neil-s-data_2023-09-25_1057';
gridscale_all = [4, 3, 2, 5];
phase_env12 = [pi/2, pi/3, 0, pi/4; pi/4, pi/5, pi/3, pi/6];
shift_x_env12 = [0, 0.1, 0, 0.2; 0.1, 0, 0.1, 0];

addpath(genpath(more_code));

n_modules = 4;

nres = 2500;
num_grids = 13456; % in a module

% load the simulated data:
module_env = zeros(n_modules,2, num_grids, nres);

for module = 1:n_modules
    gridscale = gridscale_all(module);
    file_name1 = ['simulated_grid_cells_env_1_module_s', num2str(gridscale),'.mat'];
    file_name2 = ['simulated_grid_cells_env_2_module_s', num2str(gridscale),'.mat'];
    load(file_name1)
    module_env(module, 1, :, :)  = grid_cells_env2save;
    load(file_name2)
    module_env(module, 2, :, :)  = grid_cells_env2save;
end

grid_dim = sqrt(num_grids);

%%% calculate pseudo voxels %%%
a_noise = 0.05;% 0.1;

v_max = sqrt(num_grids);
dv = 1;
shift_vec_1 = linspace(0, v_max-dv , v_max) ;
[shift_vecX,shift_vecY] = meshgrid(shift_vec_1 ,shift_vec_1 );
m_shift = shift_vecY.*shift_vecX;
median_m_shift = median(m_shift(:));

v_voxels = 2:4:36;

index_all = 1:num_grids;

n_rep = 10;

p_val_effect = zeros(length(v_voxels), n_rep);
mean_auc = zeros(length(v_voxels), n_rep);
auc_within = zeros(length(v_voxels), n_rep);
auc_between = zeros(length(v_voxels), n_rep);
grid_effect = zeros(length(v_voxels), n_rep);

snr = zeros(n_modules, length(v_voxels), n_rep);

is_plot_auc = false;

for rep = 1: n_rep

    for n=1:length(v_voxels)

        n_voxels1 = v_voxels(n);
        cells_in_voxels  = ceil(num_grids/ n_voxels1); % number of cells within a voxels

        n_voxels = length(gridscale_all)* n_voxels1; % for all modules

        env1_voxels_m = [];
        env2_voxels_m = [];

        for module = 1:n_modules

            index_random_v = randperm(num_grids);

            [env1_voxels_m1, env2_voxels_m1, snr(module, n, rep)] = cal_voxels(a_noise, index_random_v, squeeze(module_env(module, 1, :, :)), squeeze(module_env(module, 2, :, :)), n_voxels1, nres);

            env1_voxels_m = [env1_voxels_m; env1_voxels_m1];
            env2_voxels_m = [env2_voxels_m; env2_voxels_m1];

            disp(["done module ", num2str(module)])

        end

        % Calculate subspace generalization:
        cumsum_var_grid11 = cal_projection_plot(transpose(env1_voxels_m), transpose(env1_voxels_m));
        cumsum_var_grid12 = cal_projection_plot(transpose(env1_voxels_m), transpose(env2_voxels_m));

        cumsum_var_grid22 = cal_projection_plot(transpose(env2_voxels_m), transpose(env2_voxels_m));
        cumsum_var_grid21 = cal_projection_plot(transpose(env2_voxels_m), transpose(env1_voxels_m));

        cumsum_var_grid_same = 0.5 * (cumsum_var_grid11 + cumsum_var_grid22);
        cumsum_var_grid_dif = 0.5 * (cumsum_var_grid12 + cumsum_var_grid21);

        num_grid_cells = length(cumsum_var_grid_same);
        pc_fraction_grid = linspace(0,1,length(cumsum_var_grid_same)+1);

        mean_auc(n, rep) = 0.5 * ((sum(cumsum_var_grid_same) + sum(cumsum_var_grid_dif))/n_voxels);
        auc_within(n, rep) = sum(cumsum_var_grid_same) / n_voxels;
        auc_between(n, rep) = sum(cumsum_var_grid_dif) / n_voxels;
        grid_effect1 = (sum(cumsum_var_grid_same) - sum(cumsum_var_grid_dif))/n_voxels;
        grid_effect(n, rep) = grid_effect1;

        num_permutation = 10000;
        proj_auc_grid = cal_proj_auc_random_dist_2sides(transpose(env1_voxels_m), transpose(env2_voxels_m), num_permutation);

        % calculate stats:

        null_dist = sum(cumsum_var_grid_same)/n_voxels - proj_auc_grid;
        [count_null_grid, auc_dif_nul] = hist(null_dist, 100);
        p_null = count_null_grid / num_permutation;
        cdf_null = cumsum(p_null);

        n_effect = sum(auc_dif_nul < grid_effect1);
        if n_effect > 0
            p_val_grid = cdf_null(n_effect);
        else
            p_val_grid = 0;
        end
        disp([p_val_grid, n_effect])
        p_val_effect(n, rep) = p_val_grid;

        if is_plot_auc

            figure(4)
            subplot(1,2,n)

            if a_noise==0
                plot(pc_fraction_grid, [0, cumsum_var_grid_same], 'k', 'LineWidth', 2)
                hold on
                plot(pc_fraction_grid, [0, cumsum_var_grid_dif], 'g', 'LineWidth', 2)
                hold on
                plot(0:0.1:1, 0:0.1:1, ':k')
            else
                plot(pc_fraction_grid, [0, cumsum_var_grid_same], '-.k', 'LineWidth', 2)
                hold on
                plot(pc_fraction_grid, [0, cumsum_var_grid_dif], '-.g', 'LineWidth', 2)
                hold on
                plot(0:0.1:1, 0:0.1:1, ':k')
            end

        end
        disp(["done number of voxel (1D)", num2str(v_voxels(n))])
    end

    disp(["end repetition: ", num2str(rep)])
end
% calculate subject generalization stats:
auc_within_mean = mean(auc_within, 2);
auc_between_mean = mean(auc_between, 2);
grid_effect_mean = mean(grid_effect, 2);
p_val_effect_mean = mean(p_val_effect, 2);

auc_within_std = std(auc_within, 0, 2)/sqrt(n_rep);
auc_between_std = std(auc_between, 0, 2)/sqrt(n_rep);
grid_effect_std = std(grid_effect, 0, 2)/sqrt(n_rep);
p_val_effect_std = std(p_val_effect, 0, 2)/sqrt(n_rep);

v_voxels2 = [v_voxels, fliplr(v_voxels)];

figure(10)
snr_mean = mean(snr,3);
imagesc(snr_mean)
colorbar

mean_module_snr = squeeze(mean(snr, 1));
mean_mean_snr = squeeze(mean(mean_module_snr, 2));
std_mean_snr = squeeze(std(mean_module_snr, [], 2));

figure(11)
subplot(1,4,1)
plot(v_voxels, auc_within_mean, 'k', 'LineWidth', 2)
hold on
fill(v_voxels2, [auc_within_mean' + auc_within_std', fliplr(auc_within_mean' - auc_within_std')], 'k','FaceAlpha',0.3)
hold on
plot(v_voxels, auc_between_mean, 'g', 'LineWidth', 2)
hold on
fill(v_voxels2, [auc_between_mean' + auc_between_std', fliplr(auc_between_mean' - auc_between_std')], 'g','FaceAlpha',0.3)
ylabel('AUC')
xlabel('voxels/module')
ylim([0.5, 1])
subplot(1,4,2)
plot(v_voxels, grid_effect_mean, 'k', 'LineWidth', 2)
hold on
fill(v_voxels2, [grid_effect_mean' + grid_effect_std', fliplr(grid_effect_mean' - grid_effect_std')], 'k','FaceAlpha',0.3)
ylabel('AUC_{within - between}')
xlabel('voxels/module')
ylim([-0.25, 0.25])
subplot(1,4,3)
plot(v_voxels, p_val_effect_mean, 'k', 'LineWidth', 2)
hold on
fill(v_voxels2, [p_val_effect_mean' + p_val_effect_std', fliplr(p_val_effect_mean' - p_val_effect_std')], 'k', 'FaceAlpha',0.3)
xlabel('voxels/module')
ylabel('p-value_{grid effect}')
subplot(1,4,4)
plot(v_voxels, mean_mean_snr, 'k', 'LineWidth', 2)
hold on
fill(v_voxels2, [mean_mean_snr' + std_mean_snr', fliplr(mean_mean_snr' - std_mean_snr')], 'k', 'FaceAlpha',0.3)
ylabel('<SNR>')
xlabel('voxels/module')

figure(1)
subplot(2, 2, 1)
imagesc(reshape(squeeze(env1_voxels_m(1, :)), sqrt(nres), sqrt(nres)))
colorbar
subplot(2, 2, 2)
imagesc(reshape(squeeze(env1_voxels_m(2, :)), sqrt(nres), sqrt(nres)))
colorbar
subplot(2, 2, 3)
imagesc(reshape(squeeze(env1_voxels_m(3, :)), sqrt(nres), sqrt(nres)))
colorbar
subplot(2, 2, 4)
imagesc(reshape(squeeze(env1_voxels_m(4, :)), sqrt(nres), sqrt(nres)))
colorbar


for n = 1:n_modules
    figure(3)
    subplot(2,2,n)
    imagesc(reshape(mean(env1_voxels_m(1+(n-1)*n_voxels1:n*n_voxels1, :)), sqrt(nres), sqrt(nres)))
    colorbar
end

function [env1_voxels_m, env2_voxels_m, snr] = cal_voxels(a_noise, index_random_v, env1_data, env2_data, n_voxels, nres)
global cells_in_voxels

env1_voxels_m = zeros(n_voxels, nres);
env2_voxels_m = zeros(n_voxels, nres);

num_grids = length(env1_data(:,1));
grids_in_voxel = num_grids / n_voxels;

for n1 = 1:n_voxels % only random split

    if n1 == n_voxels
        selected_grids_env1 = env1_data(index_random_v(1+(n1-1)*cells_in_voxels:end), :);
        selected_grids_env2 = env2_data(index_random_v(1+(n1-1)*cells_in_voxels:end), :);
    else
        selected_grids_env1 = env1_data(index_random_v(1+(n1-1)*cells_in_voxels:cells_in_voxels*n1), :);
        selected_grids_env2 = env2_data(index_random_v(1+(n1-1)*cells_in_voxels:cells_in_voxels*n1), :);
    end
    if length(selected_grids_env1(:, 1)) < grids_in_voxel - 500 | length(selected_grids_env1(:, 1))>grids_in_voxel + 500
        disp(["env1: number of cells within voxel should be approximatly: ", num2str(grids_in_voxel)])
        warning(append("env1: number of cells: ", num2str(length(selected_grids_env1(:, 1)))))
    end
    grid2average1 = mean(selected_grids_env1);
    env1_voxels_m(n1, :) = grid2average1 + a_noise*randn(size(grid2average1));

    
    if length(selected_grids_env2(:, 1)) < grids_in_voxel - 500 | length(selected_grids_env1(:, 1))>grids_in_voxel + 500
        warning(append("env 2: number of cells: ", num2str(length(selected_grids_env2(:, 1)))))
    end

    grid2average2 = mean(selected_grids_env2);
    env2_voxels_m(n1, :) = grid2average2 + a_noise*randn(size(grid2average2));

    snr = std([grid2average1(:); grid2average2(:)])/a_noise;

end

end

