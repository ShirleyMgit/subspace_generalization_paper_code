clear all
close all
clc


% Calculates subspace generalization between simulated voxels.
% Each voxel is an average over simulated grid cells according to their
% phase. grid cells with closer phases are averaged within the same voxels.

global dv_voxels r_voxel a_noise

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
more_code = fullfile(exp_root, 'cells_data_muller');

addpath(genpath(more_code));

% load the simulated data:

gridscale = 2; %modules: gridscale = 2, 4
file_name1 = ['simulated_grid_cells_env_1_module_s', num2str(gridscale),'.mat'];
file_name2 = ['simulated_grid_cells_env_2_module_s', num2str(gridscale),'.mat'];
load(file_name1)
nres = length(grid_cells_env2save(1, :));
num_grids = length(grid_cells_env2save( :, 1));
env1_data_module1 = grid_cells_env2save;
load(file_name2)
env2_data_module1 = grid_cells_env2save;

gridscale = 4; %modules: gridscale = 2, 4
file_name1 = ['simulated_grid_cells_env_1_module_s', num2str(gridscale),'.mat'];
file_name2 = ['simulated_grid_cells_env_2_module_s', num2str(gridscale),'.mat'];
load(file_name1)
env1_data_module2 = grid_cells_env2save;
load(file_name2)
env2_data_module2 = grid_cells_env2save;

grid_dim = sqrt(num_grids);

%%% calculate pseudo voxels %%%
a_noise = 0.05;

v_max = sqrt(num_grids);
dv = 1;
shift_vec_1 = linspace(0, v_max-dv , v_max) ;
[shift_vecX,shift_vecY] = meshgrid(shift_vec_1 ,shift_vec_1 );

div_voxels = 2;
cell_in_voxels1D = grid_dim / div_voxels;

n_voxels1 = ceil(div_voxels * div_voxels);
dv_voxels = v_max / div_voxels;

cells_in_voxels = num_grids/ n_voxels1;

n_voxels = 2* n_voxels1;

v_r_random = 0:0.25:1;
mean_auc = zeros(length(v_r_random),1);
p_val_effect = zeros(length(v_r_random),1);

for n=1:length(v_r_random)

    shift_vecXv = shift_vecX(:);
    shift_vecYv = shift_vecY(:);

    r_random = v_r_random(n); % 1 if only random segregated voxels
    r_voxel = ceil(num_grids*r_random / n_voxels1);
    max_ind_random = r_voxel * n_voxels1;
    index_all = 1:num_grids;

    index_random = randperm(num_grids);
    index_random_v = index_random(1:max_ind_random );
    index_phase = index_all;
    index_phase(index_random_v) = [];

    shift_vecXv = shift_vecXv(index_phase);
    shift_vecYv = shift_vecYv(index_phase);

    [env1_voxels_m1, env2_voxels_m1] = cal_voxels(index_random_v, index_phase, env1_data_module1, env2_data_module1, shift_vecXv, shift_vecYv, index_random, n_voxels1, nres);
    [env1_voxels_m2, env2_voxels_m2] = cal_voxels(index_random_v, index_phase, env1_data_module2, env2_data_module2, shift_vecXv, shift_vecYv, index_random, n_voxels1, nres);

    env1_voxels_m = [env1_voxels_m1; env1_voxels_m2];
    env2_voxels_m = [env2_voxels_m1; env2_voxels_m2];

    % Calculate subspace generalization:
    cumsum_var_grid11 = cal_projection_plot(transpose(env1_voxels_m), transpose(env1_voxels_m));
    cumsum_var_grid12 = cal_projection_plot(transpose(env1_voxels_m), transpose(env2_voxels_m));

    cumsum_var_grid22 = cal_projection_plot(transpose(env2_voxels_m), transpose(env2_voxels_m));
    cumsum_var_grid21 = cal_projection_plot(transpose(env2_voxels_m), transpose(env1_voxels_m));

    cumsum_var_grid_same = 0.5 * (cumsum_var_grid11 + cumsum_var_grid22);
    cumsum_var_grid_dif = 0.5 * (cumsum_var_grid12 + cumsum_var_grid21);

    num_grid_cells = length(cumsum_var_grid_same);
    pc_fraction_grid = linspace(0,1,length(cumsum_var_grid_same)+1);

    mean_auc(n) = 0.5 * ((sum(cumsum_var_grid_same) + sum(cumsum_var_grid_dif))/n_voxels);

    num_permutation = 8000;
    proj_auc_grid = cal_proj_auc_random_dist_2sides(transpose(env1_voxels_m), transpose(env2_voxels_m), num_permutation);

    % calculate stats:
    grid_effect = (sum(cumsum_var_grid_same) - sum(cumsum_var_grid_dif))/n_voxels;
    null_dist = sum(cumsum_var_grid_same)/n_voxels - proj_auc_grid;
    [count_null_grid, auc_dif_nul] = hist(null_dist, 100);
    p_null = count_null_grid / num_permutation;
    cdf_null = cumsum(p_null);

    n_effect = sum(auc_dif_nul < grid_effect);
    if n_effect > 0
        p_val_grid = cdf_null(n_effect);
    else
        p_val_grid = 0;
    end

    p_val_effect(n) = p_val_grid;
end

figure(11)
subplot(1,2,1)
plot(v_r_random, mean_auc, 'k')
ylabel('mean(AUC)')
subplot(1,2,2)
plot(v_r_random, p_val_effect, 'k')
ylabel('p_{value}-subspace generalization')
xlabel('franction random')

figure(1)
subplot(2, 2, 1)
imagesc(reshape(squeeze(env1_voxels_m1(1, :)), sqrt(nres), sqrt(nres)))
colorbar
subplot(2, 2, 2)
imagesc(reshape(squeeze(env1_voxels_m1(2, :)), sqrt(nres), sqrt(nres)))
colorbar
subplot(2, 2, 3)
imagesc(reshape(squeeze(env1_voxels_m1(3, :)), sqrt(nres), sqrt(nres)))
colorbar
subplot(2, 2, 4)
imagesc(reshape(squeeze(env1_voxels_m1(4, :)), sqrt(nres), sqrt(nres)))
colorbar

figure(2)
subplot(2, 2, 1)
imagesc(reshape(squeeze(env1_voxels_m2(1, :)), sqrt(nres), sqrt(nres)))
colorbar
subplot(2, 2, 2)
imagesc(reshape(squeeze(env1_voxels_m2(2, :)), sqrt(nres), sqrt(nres)))
colorbar
subplot(2, 2, 3)
imagesc(reshape(squeeze(env1_voxels_m2(3, :)), sqrt(nres), sqrt(nres)))
colorbar
subplot(2, 2, 4)
imagesc(reshape(squeeze(env1_voxels_m2(4, :)), sqrt(nres), sqrt(nres)))
colorbar

figure(3)
subplot(2,1,1)
imagesc(reshape(mean(env1_voxels_m1), sqrt(nres), sqrt(nres)))
colorbar
subplot(2,1,2)
imagesc(reshape(mean(env1_voxels_m2), sqrt(nres), sqrt(nres)))
colorbar

figure(4)
plot(pc_fraction_grid, [0, cumsum_var_grid_same], 'k', 'LineWidth', 2)
hold on
plot(pc_fraction_grid, [0, cumsum_var_grid_dif], 'g:', 'LineWidth', 2)
hold on
plot(0:0.1:1, 0:0.1:1, ':k')

figure(10)
plot(auc_dif_nul, cdf_null)
title(['grid permutation cdf, p = ', num2str(p_val_grid)])

function [env1_voxels_m, env2_voxels_m] = cal_voxels(index_random_v, index_phase, env1_data, env2_data, shift_vecXv, shift_vecYv, index_random, n_voxels, nres)
global dv_voxels r_voxel a_noise

env1_voxels_m = zeros(n_voxels, nres);
env2_voxels_m = zeros(n_voxels, nres);

env1_data_phase = env1_data(index_phase, :);
env2_data_phase = env2_data(index_phase, :);

nv = 1;
for n1 = 1:sqrt(n_voxels)
    for n2 = 1:sqrt(n_voxels)

        ind2averageL = (shift_vecXv < n1*dv_voxels) .* (shift_vecYv(:) < n2*dv_voxels);
        ind2averageS = (shift_vecXv(:) >= (n1-1)*dv_voxels) .* (shift_vecYv(:) >= (n2-1)*dv_voxels);
        ind2average = (ind2averageS .* ind2averageL);

        grid2average = mean([env1_data_phase(ind2average>0,  :); env1_data(index_random_v(1+(nv-1)*r_voxel:r_voxel*nv), :)]);
        env1_voxels_m(n1 + (n2-1)*sqrt(n_voxels), :) = grid2average + a_noise*randn(size(grid2average));
        grid2average = mean([env2_data_phase(ind2average>0,  :); env2_data(index_random_v(1+(nv-1)*r_voxel:r_voxel*nv), :)]);
        env2_voxels_m(n1 + (n2-1)*sqrt(n_voxels), :) = grid2average + a_noise*randn(size(grid2average));
        nv = nv + 1;
    end
end

end

