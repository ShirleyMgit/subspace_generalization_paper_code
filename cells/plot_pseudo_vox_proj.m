clear all
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
more_code = fullfile(exp_root, 'cells_data_muller');
data_folder = fullfile(more_code, 'wetransfer_neil-s-data_2023-09-25_1057');
addpath(genpath(more_code));
addpath(genpath(data_folder));

grid_rats = [1, 3, 4];
num_grid_rats = length(grid_rats);

place_cells_rats = [4, 6, 8];
num_place_rats = length(place_cells_rats);

v_num = 7;
psedo_voxels_number = 2;
num_psudo_subjects = 10;
npsuedo = psedo_voxels_number*num_grid_rats;
psudo_voxels_grid = zeros(num_psudo_subjects, 2, npsuedo, 4096);

% create voxel like data:
for n=1:num_grid_rats
    load(['grid_animal',num2str(grid_rats(n)),'.mat']);
    mean_cell = mean(cells_array,2);
    std_cell_array = std(cells_array,[],2);
    norm_cells_arr = (cells_array - mean_cell) ./ std_cell_array;
    psudo_voxels1 = cal_psudo_voxels(norm_cells_arr, num_psudo_subjects, psedo_voxels_number, v_num);
    psudo_voxels_grid(:, :, (n - 1)*psedo_voxels_number + 1:n*psedo_voxels_number, :) = psudo_voxels1;
end

psudo_voxels_place = zeros(num_psudo_subjects, 2, psedo_voxels_number*num_place_rats, 4096);
for n=1:num_place_rats
    load(['place_cells_animal',num2str(place_cells_rats(n)),'.mat']);
    mean_cell = mean(cells_array,2);
    std_cell_array = std(cells_array,[],2);
    norm_cells_arr = (cells_array - mean_cell) ./ std_cell_array;
    psudo_voxels1 = cal_psudo_voxels(norm_cells_arr, num_psudo_subjects, psedo_voxels_number, v_num);
    psudo_voxels_place(:, :, (n - 1)*psedo_voxels_number + 1:n*psedo_voxels_number, :) = psudo_voxels1;
end
disp('done creating psudo voxels')

cumsum_var_11_place = zeros(num_psudo_subjects,psedo_voxels_number*num_place_rats);
cumsum_var_12_place = zeros(num_psudo_subjects,psedo_voxels_number*num_place_rats); 
cumsum_var_22_place = zeros(num_psudo_subjects,psedo_voxels_number*num_place_rats);
cumsum_var_21_place = zeros(num_psudo_subjects,psedo_voxels_number*num_place_rats);

cumsum_var_11_grid = zeros(num_psudo_subjects,psedo_voxels_number*num_grid_rats);
cumsum_var_12_grid = zeros(num_psudo_subjects,psedo_voxels_number*num_grid_rats); 
cumsum_var_22_grid = zeros(num_psudo_subjects,psedo_voxels_number*num_grid_rats);
cumsum_var_21_grid = zeros(num_psudo_subjects,psedo_voxels_number*num_grid_rats);

for s = 1:num_psudo_subjects

    psequo_place = squeeze(psudo_voxels_place(s, :, :, :));
    [cumsum_var_11_place(s, :), cumsum_var_12_place(s, :), cumsum_var_22_place(s, :),cumsum_var_21_place(s, :)] = plot_proj(psequo_place, 'g');
    psequo_grid = squeeze(psudo_voxels_grid(s, :, :, :));
    [cumsum_var_11_grid(s, :), cumsum_var_12_grid(s, :), cumsum_var_22_grid(s, :),cumsum_var_21_grid(s, :)] = plot_proj(psequo_grid, 'k');

end

cumsum_var_same_grid = 0.5 * (cumsum_var_11_grid + cumsum_var_22_grid);
cumsum_var_dif_grid = 0.5 * (cumsum_var_12_grid + cumsum_var_21_grid);
cumsum_var_same_place = 0.5 * (cumsum_var_11_place+ cumsum_var_22_place);
cumsum_var_dif_place = 0.5 * (cumsum_var_12_place + cumsum_var_21_place);

mean_grid_same = mean(cumsum_var_same_grid);
mean_grid_dif = mean(cumsum_var_dif_grid);
mean_place_same = mean(cumsum_var_same_place);
mean_place_dif = mean(cumsum_var_dif_place);

std_grid_same = std(cumsum_var_same_grid)/sqrt(num_psudo_subjects);
std_grid_dif = std(cumsum_var_dif_grid)/sqrt(num_psudo_subjects);
std_place_same = std(cumsum_var_same_place)/sqrt(num_psudo_subjects);
std_place_dif = std(cumsum_var_dif_place)/sqrt(num_psudo_subjects);

ci_grid_same_p = mean_grid_same + std_grid_same;
ci_grid_same_m = mean_grid_same- std_grid_same;
ci_grid_dif_p = mean_grid_dif + std_grid_dif;
ci_grid_dif_m = mean_grid_dif - std_grid_dif;

ci_place_same_p = mean_place_same + std_place_same;
ci_place_same_m = mean_place_same- std_place_same;
ci_place_dif_p = mean_place_dif + std_place_dif;
ci_place_dif_m = mean_place_dif - std_place_dif;

mean_grid_11 = mean(cumsum_var_11_grid);
mean_grid_12 = mean(cumsum_var_12_grid);
mean_grid_21 = mean(cumsum_var_21_grid);
mean_grid_22 = mean(cumsum_var_22_grid);

std_grid_11 = std(cumsum_var_11_grid)/sqrt(num_psudo_subjects);
std_grid_12 = std(cumsum_var_12_grid)/sqrt(num_psudo_subjects);
std_grid_21 = std(cumsum_var_21_grid)/sqrt(num_psudo_subjects);
std_grid_22 = std(cumsum_var_22_grid)/sqrt(num_psudo_subjects);

mean_place_11 = mean(cumsum_var_11_place);
mean_place_12 = mean(cumsum_var_12_place);
mean_place_21 = mean(cumsum_var_21_place);
mean_place_22 = mean(cumsum_var_22_place);

std_place_11 = std(cumsum_var_11_place)/sqrt(num_psudo_subjects);
std_place_12 = std(cumsum_var_12_place)/sqrt(num_psudo_subjects);
std_place_21 = std(cumsum_var_21_place)/sqrt(num_psudo_subjects);
std_place_22 = std(cumsum_var_22_place)/sqrt(num_psudo_subjects);

ci_grid_11p = mean_grid_11 + std_grid_11;
ci_grid_11m = mean_grid_11 - std_grid_11;
ci_grid_12p = mean_grid_12 + std_grid_12;
ci_grid_12m = mean_grid_12 - std_grid_12;
ci_grid_21p = mean_grid_21 + std_grid_21;
ci_grid_21m = mean_grid_21 - std_grid_21;
ci_grid_22p = mean_grid_22 + std_grid_22;
ci_grid_22m = mean_grid_22 - std_grid_22;

ci_place_11p = mean_place_11 + std_place_11;
ci_place_11m = mean_place_11 - std_place_11;
ci_place_12p = mean_place_12 + std_place_12;
ci_place_12m = mean_place_12 - std_place_12;
ci_place_21p = mean_place_21 + std_place_21;
ci_place_21m = mean_place_21 - std_place_21;
ci_place_22p = mean_place_22 + std_place_22;
ci_place_22m = mean_place_22 - std_place_22;

x = linspace(0, 1 ,npsuedo);
figure(10)
plot(x, mean_grid_same, 'k', 'LineWidth', 2)
hold on
plot(x, mean_grid_dif, ':k', 'LineWidth', 2)
hold on
plot(x, mean_place_same, 'g', 'LineWidth', 2)
hold on
plot(x, mean_place_dif, ':g', 'LineWidth', 2)
hold on
patch([x'; flipud(x')], [ci_grid_same_m'; flipud(ci_grid_same_p')], [1 1 1]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)
hold on
patch([x'; flipud(x')], [ci_grid_dif_m'; flipud(ci_grid_dif_p')], [1 1 1]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)
hold on
patch([x'; flipud(x')], [ci_place_same_m'; flipud(ci_place_same_p')], [0 1 0]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)
hold on
patch([x'; flipud(x')], [ci_place_dif_m'; flipud(ci_place_dif_p')], [0 1 0]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)
ylim([0, 1])
legend('grid within-env', 'grid across-env','place within-env', 'place across-env')

figure(1)
plot(x, mean_grid_11, 'g')
hold on
patch([x'; flipud(x')], [ci_grid_11m'; flipud(ci_grid_11p')], [0 1 0]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)
hold on
plot(x, mean_grid_12, 'b')
hold on
patch([x'; flipud(x')], [ci_grid_12m'; flipud(ci_grid_12p')], [0 0 1]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)
hold on
plot(x, mean_place_11, 'r')
hold on
patch([x'; flipud(x')], [ci_place_11m'; flipud(ci_place_11p')], [1 0 0]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)
hold on
plot(x, mean_place_12, 'k')
hold on
patch([x'; flipud(x')], [ci_place_12m'; flipud(ci_place_12p')], [1 1 1]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)


figure(2)
plot(x, mean_grid_22, 'g')
hold on
patch([x'; flipud(x')], [ci_grid_22m'; flipud(ci_grid_22p')], [0 1 0]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)
hold on
plot(x, mean_grid_21, 'b')
hold on
patch([x'; flipud(x')], [ci_grid_21m'; flipud(ci_grid_21p')], [0 0 1]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)
hold on
plot(x, mean_place_22, 'r')
hold on
patch([x'; flipud(x')], [ci_place_22m'; flipud(ci_place_22p')], [1 0 0]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)
hold on
plot(x, mean_place_21, 'k')
hold on
patch([x'; flipud(x')], [ci_place_21m'; flipud(ci_place_21p')], [1 1 1]*0.5, 'EdgeColor','none', 'FaceAlpha',0.2)

figure(10)
for f=1:6
    subplot(6,2,f)
    imagesc(reshape(squeeze(psudo_voxels_place(1,1,f,:)), 64,64))
    colorbar
    title('R')
    subplot(6,2,f+6)
    imagesc(reshape(squeeze(psudo_voxels_place(1,2,f,:)), 64,64))
    colorbar
    title('VR')
end

figure(11)
for f=1:6
    subplot(6,2,f)
    imagesc(reshape(squeeze(psudo_voxels_grid(1,1,f,:)), 64,64))
    colorbar
    title('R')
    subplot(6,2,f+6)
    imagesc(reshape(squeeze(psudo_voxels_grid(1,2,f,:)), 64,64))
    colorbar
    title('VR')
end

function psudo_voxels1 = cal_psudo_voxels(cells_array, num_psudo_subjects, psedo_voxels_number, v_num)
    env1_data = squeeze(cells_array(1, :, :));
    env1_data(isnan(env1_data)) = 0;
    env2_data = squeeze(cells_array(2, :, :));
    env2_data(isnan(env2_data)) = 0;
    
    psudo_voxels1 = zeros(num_psudo_subjects, 2, psedo_voxels_number, 4096);
    num_cells = length(env1_data(:, 1));
    for s = 1:num_psudo_subjects
        p = randperm(num_cells);
       % v_num = floor(num_cells / psedo_voxels_number);
        for n = 1:psedo_voxels_number
            env1_data_p = env1_data(p, :);
            env2_data_p = env2_data(p, :);
            psudo_voxels1(s, 1, n, :) = mean(env1_data_p((1+(n-1)*v_num):(n*v_num) , :));
            psudo_voxels1(s, 2, n, :) = mean(env2_data_p((1+(n-1)*v_num):(n*v_num) , :));
        end
    end
end


function [cumsum_var_11, cumsum_var_12, cumsum_var_22,cumsum_var_21] = plot_proj(psudo_voxels, c)

    env1_data = squeeze(psudo_voxels(1, :, :));
    env2_data = squeeze(psudo_voxels(2, :, :));

    num_cells1 = length(env1_data(:, 1));
    cumsum_var_11 = cal_projection_plot(transpose(env1_data), transpose(env1_data));
    cumsum_var_12 = cal_projection_plot(transpose(env1_data), transpose(env2_data));

    num_cells2 = length(env2_data(:, 1));
    cumsum_var_22 = cal_projection_plot(transpose(env2_data), transpose(env2_data));
    cumsum_var_21 = cal_projection_plot(transpose(env2_data), transpose(env1_data));

    % figure(100)
    % hold on
    % plot(cumsum_var_11, c)
    % hold on
    % plot(cumsum_var_12, ['-.', c])
    % 
    % figure(101)
    % hold on
    % plot(cumsum_var_22, c)
    % hold on
    % plot(cumsum_var_21, ['-.', c])
end