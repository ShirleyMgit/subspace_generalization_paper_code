% checking subjects responses during the fMRI session
% both catch questions and last question in each block
clear all
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
data_path = fullfile(exp_root, 'alon', 'beh');

n_catch = 5;
num_subjects = 28;
num_runs = 4;
num_maps= 4; 
num_blocks = 5;
all_responses = zeros(num_subjects, num_runs, num_maps);
catch_iscor = zeros(num_subjects, num_runs, num_maps);

for s = 1:num_subjects
    if s < 10
        sub_dir = ['sub-0', num2str(s)];
    else
        sub_dir = ['sub-', num2str(s)];
    end
    for r = 1:num_runs
        load(fullfile(data_path, sub_dir,['run-0', num2str(r)]))
        maps = RSA_maps.maps;
        for b = 1:num_blocks
            if maps(b) ~= 5
                graph_num = RSA_maps(b).block_map;
                catch_iscor(s, r, graph_num) = sum(RSA_maps(b).res.choice)/n_catch; % 1 correcr, 0 incorect
                if s > 1 % I asked the first participants when this set appears during the frist day as she did not describe any clusters in the briefing
                    all_responses(s, r, graph_num) = RSA_maps(b).pic.whichDayS;
                else
                    if RSA_maps(b).pic.whichDayS > -1
                        all_responses(s, r, graph_num) = 1 - RSA_maps(b).pic.whichDayS;
                    else
                        all_responses(s, r, graph_num) = RSA_maps(b).pic.whichDayS;
                    end
                end
            end
        end
    end
end

all_hex = reshape(all_responses(:, :, 1:2), num_subjects, 8);
all_cluster = reshape(all_responses(:, :, 3:4), num_subjects, 8);

hex_catch = reshape(catch_iscor(:, :, 1:2), num_subjects, 8);
cl_catch = reshape(catch_iscor(:, :, 3:4), num_subjects, 8);

% figure(1)
% subplot(2,1,1)
% imagesc(all_hex)
% colorbar
% subplot(2,1,2)
% imagesc(all_cluster)
% colorbar

% calcualte fraction correct:
hex_correct = sum(all_hex == 0, 2)/8;
cluster_correct = sum(all_cluster == 1, 2)/8;
hex_correct_catch = mean(cl_catch, 2);
cluster_correct_catch = mean(hex_catch, 2);

sem_cluster = std(cluster_correct)/sqrt(num_subjects);
sem_hex = std(cluster_correct)/sqrt(num_subjects);

sem_cluster_catch = std(cluster_correct_catch)/sqrt(num_subjects);
sem_hex_catch = std(hex_correct_catch)/sqrt(num_subjects);

% a = get(gca, 'XTickLabel');
% set(gca, 'XTickLabel', a, 'fontsize', 12);
figure(1)
subplot(1,2,1)
bar([1,2], [mean(hex_correct_catch), mean(cluster_correct_catch)])
hold on
errorbar([1,2], [mean(hex_correct_catch), mean(cluster_correct_catch)], [sem_hex_catch/2, sem_cluster_catch/2],'k.')
xticklabels({'hex', 'cluster'})
ylabel('p_{correct}','FontSize', 12)
title('within block catch trials')
subplot(1,2,2)
bar([1,2], [mean(hex_correct), mean(cluster_correct)])
hold on
errorbar([1,2], [mean(hex_correct), mean(cluster_correct)], [sem_hex/2, sem_cluster/2],'k.')
xticklabels({'hex', 'cluster'})
title('graph type question','FontSize', 12)

[h_hex_catch, p_hex_catch, ~, stats_hex_catch] = ttest(hex_correct_catch, 0.5, 'tail', 'right')
[h_cl_catch, p_cl_catch, ~, stats_cl_catch] = ttest(cluster_correct_catch, 0.5, 'tail', 'right')

[h_hex, p_hex, ~, stats_hex] = ttest(hex_correct, 0.5, 'tail', 'right')
[h_cl, p_cl, ~, stats_cl] = ttest(cluster_correct, 0.5, 'tail', 'right')