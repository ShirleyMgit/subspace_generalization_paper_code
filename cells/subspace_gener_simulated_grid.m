clear all
close all
clc


% Calculates subspace generalization between simulated voxels.
% Each voxel is an average over simulated grid cells according to their
% phase. grid cells with closer phases are averaged within the same voxels.


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
env1_data = grid_cells_env2save;
load(file_name2)
env2_data = grid_cells_env2save;

%%% calculate pseudo voxels %%%

dx = 0.01;
max_x = 2/sqrt(3);%because of the gridscale/ actual grid scale ratio (reciprocal lattice)
max_y = max_x*sqrt(3) / 2;
dy = dx* sqrt(3) / 2;
shift_vec_x1 = linspace(0, max_x -dx , ceil(max_x/dx)) ;
shift_vec_y1 = linspace(0, max_y -dy, ceil(max_y/dy)) ;
[shift_vecX,shift_vecY] = meshgrid(shift_vec_x1,shift_vec_y1); 

div_voxels = 4;
cell_in_voxels1D = sqrt(num_grids) / div_voxels;

n_voxels = ceil(div_voxels * div_voxels);
dx_voxels = max_x / div_voxels;
dy_voxels = max_y / div_voxels;

cells_in_voxels = num_grids/ n_voxels;

env1_voxels_m = zeros(n_voxels, nres);
env2_voxels_m = zeros(n_voxels, nres);

for n1 = 1:sqrt(n_voxels)
    for n2 = 1:sqrt(n_voxels)
        ind2averageL = (shift_vecX(:) <= n1*dx_voxels) .* (shift_vecY(:) <= n2*dy_voxels);
        ind2averageS = (shift_vecX(:) > (n1-1)*dx_voxels) .* (shift_vecY(:) > (n2-1)*dy_voxels);
        ind2average = (ind2averageS .* ind2averageL);
        grid2average = mean(env1_data(ind2average>0,  :));
        env1_voxels_m(n1 + (n2-1)*sqrt(n_voxels), :) = grid2average;
        grid2average = mean(env2_data(ind2average>0,  :));
        env2_voxels_m(n1 + (n2-1)*sqrt(n_voxels), :) = grid2average;
    end
end

figure(1)
subplot(2, 2, 1)
imagesc(reshape(squeeze(env1_voxels_m(1, :)), sqrt(nres), sqrt(nres)))
subplot(2, 2, 2)
imagesc(reshape(squeeze(env1_voxels_m(2, :)), sqrt(nres), sqrt(nres)))
subplot(2, 2, 3)
imagesc(reshape(squeeze(env1_voxels_m(3, :)), sqrt(nres), sqrt(nres)))
subplot(2, 2, 4)
imagesc(reshape(squeeze(env1_voxels_m(4, :)), sqrt(nres), sqrt(nres)))

% Calculate subspace generalization:
cumsum_var_grid11 = cal_projection_plot(transpose(env1_voxels_m), transpose(env1_voxels_m));
cumsum_var_grid12 = cal_projection_plot(transpose(env1_voxels_m), transpose(env2_voxels_m));

cumsum_var_grid22 = cal_projection_plot(transpose(env2_voxels_m), transpose(env2_voxels_m));
cumsum_var_grid21 = cal_projection_plot(transpose(env2_voxels_m), transpose(env1_voxels_m));

cumsum_var_grid_same = 0.5 * (cumsum_var_grid11 + cumsum_var_grid22);
cumsum_var_grid_dif = 0.5 * (cumsum_var_grid12 + cumsum_var_grid21);

num_grid_cells = length(cumsum_var_grid_same);
pc_fraction_grid = linspace(0,1,length(cumsum_var_grid_same)+1);
figure(2)
plot(pc_fraction_grid, [0, cumsum_var_grid_same], 'k', 'LineWidth', 2)
hold on
plot(pc_fraction_grid, [0, cumsum_var_grid_dif], 'g:', 'LineWidth', 2)

num_permutation = 500;
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

figure(10)
plot(auc_dif_nul, cdf_null)
title(['grid permutation cdf, p = ', num2str(p_val_grid)])

