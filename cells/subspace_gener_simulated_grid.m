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
v_max = sqrt(num_grids);
dv = 1;
shift_vec_1 = linspace(0, v_max-dv , v_max) ;
[shift_vecX,shift_vecY] = meshgrid(shift_vec_1 ,shift_vec_1 ); 

div_voxels = 4;
cell_in_voxels1D = sqrt(num_grids) / div_voxels;

n_voxels = ceil(div_voxels * div_voxels);
cells_in_voxels = num_grids/ n_voxels;
dv_voxels = v_max / div_voxels;


env1_voxels_m = zeros(n_voxels, nres);
env2_voxels_m = zeros(n_voxels, nres);
nall = zeros(n_voxels, 1);
nv = 1;
for n1 = 1:div_voxels
    for n2 = 1:div_voxels
        ind2averageL = (shift_vecX(:) < n1*dv_voxels) .* (shift_vecY(:) < n2*dv_voxels);
        ind2averageS = (shift_vecX(:) >= (n1-1)*dv_voxels) .* (shift_vecY(:) >= (n2-1)*dv_voxels);
        ind2average = (ind2averageS .* ind2averageL);
        
        nall(nv) = sum(sum( ind2average));
        grid2average = mean(env1_data(ind2average>0,  :));
        env1_voxels_m(n1 + (n2-1)*sqrt(n_voxels), :) = grid2average;

        grid2average = mean(env2_data(ind2average>0,  :));
        env2_voxels_m(n1 + (n2-1)*sqrt(n_voxels), :) = grid2average;
        nv = nv+1;
    end
end

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

figure(2)
subplot(3,1,1)
imagesc(reshape(mean(env1_voxels_m), sqrt(nres), sqrt(nres)))
title('pseudo voxels 1')
colorbar
subplot(3,1,2)
imagesc(reshape(mean(env2_voxels_m), sqrt(nres), sqrt(nres)))
title('pseudo voxels 2')
colorbar
subplot(3,1,3)
imagesc(reshape(mean(env1_data), sqrt(nres), sqrt(nres)))
title('simulated grids')
colorbar
% Calculate subspace generalization:
cumsum_var_grid11 = cal_projection_plot(transpose(env1_voxels_m), transpose(env1_voxels_m));
cumsum_var_grid12 = cal_projection_plot(transpose(env1_voxels_m), transpose(env2_voxels_m));

cumsum_var_grid22 = cal_projection_plot(transpose(env2_voxels_m), transpose(env2_voxels_m));
cumsum_var_grid21 = cal_projection_plot(transpose(env2_voxels_m), transpose(env1_voxels_m));

cumsum_var_grid_same = 0.5 * (cumsum_var_grid11 + cumsum_var_grid22);
cumsum_var_grid_dif = 0.5 * (cumsum_var_grid12 + cumsum_var_grid21);

num_grid_cells = length(cumsum_var_grid_same);
pc_fraction_grid = linspace(0,1,length(cumsum_var_grid_same)+1);
figure(3)
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

