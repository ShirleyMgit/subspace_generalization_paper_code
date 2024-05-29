clear all
close all
clc

% Calculates subspace generalization between simulated voxels.
% Each voxel is an average over simulated grid cells. Fraction of cells are selected
% according to their phase such that grid cells with closer phases are averaged within the same voxels
% and fraction that is selected randomly.
% noise is added to the voxels' activity map
% 2 voxels per module, 4 modules
global dv_voxels r_voxel flag_2D_split median_m_shift

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
more_code = fullfile(exp_root, 'cells_data_muller');

saving_folder = 'C:\Users\User\Documents\fMRI_EXP\cells_data_muller\wetransfer_neil-s-data_2023-09-25_1057';
gridscale_all = [4, 3, 2, 5];
phase_env12 = [pi/2, pi/3, 0, pi/4; pi/4, pi/5, pi/3, pi/6];
shift_x_env12 = [0, 0.1, 0, 0.2; 0.1, 0, 0.1, 0];

addpath(genpath(more_code));

flag_2D_split = false; % split grids into voxels both by X and Y (if false split only by X)
n_modules = 4;
n_voxels1 = 2;

nres = 2500;
num_grids = 13456; % in module

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
v_noise = [0, 0.1];

v_max = sqrt(num_grids);
dv = 1;
shift_vec_1 = linspace(0, v_max-dv , v_max) ;
[shift_vecX,shift_vecY] = meshgrid(shift_vec_1 ,shift_vec_1 );
m_shift = shift_vecY.*shift_vecX;
median_m_shift = median(m_shift(:));
m_shift_cut = m_shift < median_m_shift;

if flag_2D_split % split a module according to phase on the XY plane (otherwise using x-axis alone)
    div_voxels = sqrt(n_voxels1);
    cell_in_voxels1D = grid_dim / div_voxels;
else
    div_voxels = n_voxels1;
end

dv_voxels = v_max / div_voxels;

cells_in_voxels = num_grids/ n_voxels1;

n_voxels = length(gridscale_all)* n_voxels1;

v_r_random = [0, 1];
index_all = 1:num_grids;

is_plot_auc = true;
is_plot_voxels = false;
for isnoise = 1:length(v_noise)
    a_noise = v_noise(isnoise);

    for n=1:length(v_r_random)

        r_random = v_r_random(n); % 1 if only random segregated voxels
        r_voxel = ceil(num_grids*r_random / n_voxels1);
        max_ind_random = r_voxel * n_voxels1;

        env1_voxels_m = [];
        env2_voxels_m = [];

        for module = 1:n_modules

            shift_vecXv = shift_vecX(:);
            shift_vecYv = shift_vecY(:);

            index_random = randperm(num_grids);
            index_random_v = index_random(1:max_ind_random);
            index_phase = index_all;
            index_phase(index_random_v) = [];

            shift_vecXv = shift_vecXv(index_phase);
            shift_vecYv = shift_vecYv(index_phase);
            [env1_voxels_m1, env2_voxels_m1] = cal_voxels(a_noise, index_random_v, index_phase, squeeze(module_env(module, 1, :, :)), squeeze(module_env(module, 2, :, :)), shift_vecXv, shift_vecYv, n_voxels1, nres);

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

        grid_effect1 = (sum(cumsum_var_grid_same) - sum(cumsum_var_grid_dif))/n_voxels;

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

        if is_plot_auc

            figure(4)
            subplot(1,2,n)

            if a_noise==0
                plot(pc_fraction_grid, [0, cumsum_var_grid_same], 'k', 'LineWidth', 2)
                hold on
                plot(pc_fraction_grid, [0, cumsum_var_grid_dif], 'g', 'LineWidth', 2)
                hold on
                plot(0:0.1:1, 0:0.1:1, ':k')
                title(['fraction random = ', num2str(r_random)])
            else
                plot(pc_fraction_grid, [0, cumsum_var_grid_same], '-.k', 'LineWidth', 2)
                hold on
                plot(pc_fraction_grid, [0, cumsum_var_grid_dif], '-.g', 'LineWidth', 2)
                hold on
                plot(0:0.1:1, 0:0.1:1, ':k')
            end

        end
        disp(["done random f", num2str(n)])
        if is_plot_voxels
            figure(n)
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
        end
    end
    disp("done isnoise")
end


for n = 1:n_modules
    figure(3)
    subplot(2,2,n)
    imagesc(reshape(mean(env1_voxels_m(1+(n-1)*n_voxels1:n*n_voxels1, :)), sqrt(nres), sqrt(nres)))
    colorbar
end

function [env1_voxels_m, env2_voxels_m] = cal_voxels(a_noise, index_random_v, index_phase, env1_data, env2_data, shift_vecXv, shift_vecYv, n_voxels, nres)
global dv_voxels r_voxel flag_2D_split median_m_shift

env1_voxels_m = zeros(n_voxels, nres);
env2_voxels_m = zeros(n_voxels, nres);

env1_data_phase = env1_data(index_phase, :);
env2_data_phase = env2_data(index_phase, :);

num_grids = length(env1_data(:,1));
grids_in_voxel = num_grids / n_voxels;

shift_vec = shift_vecXv.*shift_vecYv;
if flag_2D_split  % seperate into voxels according to both X and Y
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
else
    for n1 = 1:n_voxels

        if n1==1
            ind2averageL = shift_vec < median_m_shift;
            ind2averageS = shift_vec < median_m_shift;
        else
            ind2averageL = shift_vec >= median_m_shift;
            ind2averageS = shift_vec >= median_m_shift;
        end
        ind2average = (ind2averageS .* ind2averageL);

        selected_grids_env1 = [env1_data_phase(ind2average>0,  :); env1_data(index_random_v(1+(n1-1)*r_voxel:r_voxel*n1), :)];
        if length(selected_grids_env1(:, 1)) < 6000 | length(selected_grids_env1(:, 1))>7000
            disp(["env1: number of cells within voxel should be approximatly: ", num2str(grids_in_voxel)])
            disp(["env1: voxels according to phase: ", num2str(sum(sum(ind2average)))])
            warning(append("env1: number of cells: ", num2str(length(selected_grids_env1(:, 1)))))
        end
        grid2average = mean(selected_grids_env1);
        env1_voxels_m(n1, :) = grid2average + a_noise*randn(size(grid2average));

        selected_grids_env2 = [env2_data_phase(ind2average>0,  :); env2_data(index_random_v(1+(n1-1)*r_voxel:r_voxel*n1), :)];
        if length(selected_grids_env2(:, 1)) < 6000 | length(selected_grids_env1(:, 1))>7000
            disp(["env2: voxels according to phase: ", num2str(sum(sum(ind2average)))])
            warning(append("env 2: number of cells: ", num2str(length(selected_grids_env2(:, 1)))))
        end

        grid2average = mean(selected_grids_env2);
        env2_voxels_m(n1, :) = grid2average + a_noise*randn(size(grid2average));

    end
end
end

