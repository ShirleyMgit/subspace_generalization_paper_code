clear all
close all
clc

% plot simulated grid cells

% load the simulated data:
% modules: gridscale = 2, 4
gridscale = 2;
file_name1 = ['simulated_grid_cells_env_1_module_s', num2str(gridscale),'.mat'];
file_name2 = ['simulated_grid_cells_env_2_module_s', num2str(gridscale),'.mat'];
load(file_name1)
nres = length(grid_cells_env2save(1, :));
num_grids = length(grid_cells_env2save( :, 1));
env1_data_module1 = grid_cells_env2save;
load(file_name2)
env2_data_module1 = grid_cells_env2save;

gridscale = 4;
file_name1 = ['simulated_grid_cells_env_1_module_s', num2str(gridscale),'.mat'];
file_name2 = ['simulated_grid_cells_env_2_module_s', num2str(gridscale),'.mat'];
load(file_name1)
env1_data_module2 = grid_cells_env2save;
load(file_name2)
env2_data_module2 = grid_cells_env2save;

figure(1)
subplot(2,4,1)
imagesc(reshape(env1_data_module1(1,:,:), 50, 50))
colorbar
subplot(2,4,2)
imagesc(reshape(env1_data_module1(100,:,:), 50, 50))
colorbar
subplot(2,4,3)
imagesc(reshape(env2_data_module1(1,:,:), 50, 50))
colorbar
subplot(2,4,4)
imagesc(reshape(env2_data_module1(100,:,:), 50, 50))
colorbar
title("module 1")
subplot(2,4,5)
imagesc(reshape(env1_data_module2(1,:,:), 50, 50))
colorbar
subplot(2,4,6)
imagesc(reshape(env1_data_module2(100,:,:), 50, 50))
colorbar
subplot(2,4,7)
imagesc(reshape(env2_data_module2(1,:,:), 50, 50))
colorbar
subplot(2,4,8)
imagesc(reshape(env2_data_module2(100,:,:), 50, 50))
colorbar
title("module 2")