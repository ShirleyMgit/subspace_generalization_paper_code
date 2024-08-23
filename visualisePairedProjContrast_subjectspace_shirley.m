clear
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
root = 'C:\Users\User\Documents\fMRI_EXP\Alon\';

% where to get the data from - only one should be true
peakFlag = false;
avgEffectT2_EC_Mask = true;
OFC_T = false;
visual_effect = false;
%avgEffect0p01Mask = false;
%avgEffect0p05Mask = true;
alon_mask = false;
alon_peak = false;
EHR_julich = false;

if avgEffectT2_EC_Mask
    masks_path = fullfile(root,'subspaceGener','fsl_normalization','groupStats','masks');
else
    masks_path = fullfile(root,'subspaceGener','fsl_normalization','groupStats','masks','masks_in_sub_numbers');
end
spm_path =  "C:\Users\User\Documents\spm12";
addpath(spm_path);

subjects = cell(1,28);
for iSub = 1:length(subjects)
    subjects{iSub}= ['sub-' num2str(iSub,'%02.f')];   
end

plotSwarmChartFlag = false; 

if peakFlag
    if alon_peak
        titleStr = 'alon_peak';
        peakVoxIndsFSL = [32, 60,22];
        peakVoxInds = peakVoxIndsFSL + 1; % Matlab indexing
    else
        titleStr = 'peak';
        peakVoxIndsFSL = [31,58,16];
        peakVoxInds = peakVoxIndsFSL + 1; % Matlab indexing
    end

else
    if avgEffectT2_EC_Mask
        titleStr = 'EC T2';
        mask_name = 'tStat_hexOnHexMinusHexOnClusterT2_TemporalFusiform_mask_R_warp.nii';
    elseif OFC_T
        titleStr = 'semi OFC T2';
        mask_name = 'tStat_projSameStr_allOthersT2_Subcallosal_mask.nii'; 
    elseif visual_effect 
        titleStr='LOC'
        mask_name = 'tStat_visual_sameStructSameStimMinusSameStructDiffStimT2_LOC_mask_L.nii';
    elseif EHR_julich
        titleStr = 'right EC masks';
        mask = fullfile(exp_root,'ROImasks','Entorhinal_R_juelich.nii');
    elseif alon_mask
        titleStr = 'Alon effect';
        mask =  fullfile(root,'masks','mni', 'relationalStructure_smth5_rh_nPerm10000_clstrTh3p1_clusterm_mask.nii');
    end
end
allElem = 1:16;
proHexElem = [1, 2, 5, 6];
proHexClElem = [9, 10, 13, 14];
proClClElem = [11, 12, 15, 16];
proClHexElem = [3, 4, 7, 8];
proAllSameElm = [1, 6, 11, 16];
proSameStrElm = [2, 5, 12, 15];

allData = nan(length(allElem),length(subjects));
mask_size = zeros(length(subjects), 1);
for iSub=1:length(subjects)
    sub = subjects{iSub};
    subspaceGenerDir = fullfile(root,'subspaceGener', sub);
    if ~peakFlag
        sub_mask_dir = fullfile(masks_path, sub);
        mask = fullfile(sub_mask_dir, mask_name);
        maskData = niftiread(mask);
    end
    % Full matrix
    for iElem = 1:16
        fname{iElem} = fullfile(subspaceGenerDir,['L100_projMat' num2str(allElem(iElem)) '.nii']);
        V_all = spm_vol(fname{iElem});
        if peakFlag
            allData(iElem,iSub) = spm_sample_vol(V_all,peakVoxInds(1),peakVoxInds(2),peakVoxInds(3),0);
        else
            tmp = niftiread(fname{iElem});%spm_read_vols(V_all);
            % mask by effect thresholded mask
            tmp = tmp(:);
            maskData = logical(maskData(:)>0.99);    
            mask_size(iSub) = sum(maskData(:),"omitnan");
            blob = tmp(maskData);
            % take the mean in th the mask
            allData(iElem,iSub) = mean(blob,"omitnan");             
        end
    end
end

proHex_allData = allData(proHexElem,:);
proHexCl_allData = allData(proHexClElem,:);
proClCl_allData = allData(proClClElem, :);
proClHex_allData = allData(proClHexElem, :);
proAllSame_allData = allData(proAllSameElm, :);
proSameStr_allData = allData(proSameStrElm, :);

mean_visualeffect = mean(proAllSame_allData - proSameStr_allData);
figure(200)
plot_warm(1, mean_visualeffect, '');

mean_matrixHex_suqres = mean(proHex_allData);
mean_matrixHexCl_suqres = mean(proHexCl_allData);
mean_matrixClCl_suqres = mean(proClCl_allData);
mean_matrixClHex_suqres = mean(proClHex_allData);

mean_main_effect = mean_matrixHex_suqres - mean_matrixHexCl_suqres;

% figure(100)
% plot_warm(1, mean_main_effect, '');

[h_hex_cl, p_hex_cl, ~, stats] = ttest(mean_matrixHex_suqres, mean_matrixHexCl_suqres, 'Tail', 'right');
[h_cl_hex, p_cl_hex, ~, stats_clhex] = ttest(mean_matrixHex_suqres, mean_matrixClHex_suqres, 'Tail', 'right');
[h_clcl_hex, p_clcl_hex] = ttest(mean_matrixClCl_suqres, mean_matrixHexCl_suqres, 'Tail', 'right');
[h_clcl_clhex, p_clcl_clhex] = ttest(mean_matrixClCl_suqres, mean_matrixClHex_suqres, 'Tail', 'right');
[h_clcl_clhex_sym, p_clcl_clhex_sym] = ttest(mean_matrixClCl_suqres, mean_matrixClHex_suqres);



cl_contrast = mean_matrixClCl_suqres - mean_matrixClHex_suqres;
% load(fullfile(dataDir,'pConCsCcor2.mat'))
% load(fullfile(dataDir,'pConCsCncor2.mat'))
% load(fullfile(dataDir,'pConCsCncor1.mat'))
% 
% cor_minus_notcor = pConCsCcor2 - pConCsCncor2;
% 
% [c, p] = corrcoef(cl_contrast, cor_minus_notcor);
% [c1, p1] = corrcoef(cl_contrast, pConCsCncor1);
% [c2, p2] = corrcoef(cl_contrast, pConCsCncor2);

if plotSwarmChartFlag
    nSub=28;
    
    proHex = mean(proHex_allData, 1);
    proHexCl = mean(proHexCl_allData, 1);

    proClCl = mean(proClCl_allData, 1);
    proClHex= mean(proClHex_allData, 1);

    % Sample data
    figure(1010)
    plot_warm(1, proHex - proHexCl, 'HexOnHex - HexOnComm');
    title(['p_{val}=', num2str(p_hex_cl)])
    plot_warm(2, proHex - proClHex, 'HexOnHex - CommOnHex');
    title(['p_{val}=', num2str(p_cl_hex)])
    plot_warm(3, proClCl - proClHex, 'CommOnComm - CommOnHex');
    title(['p_{val}=', num2str(p_clcl_clhex)])
    plot_warm(4, proClCl - proHexCl, 'CommOnComm - HexOnComm');
    title(['p_{val}=',num2str(p_clcl_hex)])

end

%%

projMatAllSubj = reshape(allData,[4,4,length(subjects)]);

%[c, p ] = corrcoef(hex11_data, cluster11_data);


projMatAllSubj2 = projMatAllSubj;
projMatAllSubj2(: ,4 ,:) = projMatAllSubj(:, 3, :);
projMatAllSubj2(: ,3 ,:) = projMatAllSubj(:, 4, :);
v = projMatAllSubj2(4, :, :);
projMatAllSubj2(4 ,: ,:) = projMatAllSubj2(3, :, :);
projMatAllSubj2(3 ,: ,:) = v;

projMatAll26Subj = projMatAllSubj2(:,:,[1, 3:end]);
figure(100)
imagesc(mean(projMatAll26Subj, 3))
colormap gray
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij

projMatAllSubjC = projMatAllSubj2;
projMatMean = squeeze(mean(projMatAllSubjC,3));

figure(1000)
hold on
title([titleStr ', mean'])
imagesc(projMatMean);
colormap gray
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
%colorbar
axis ij


hex11_data = squeeze(projMatAllSubjC(1, 1, :));
hex22_data = squeeze(projMatAllSubjC(2, 2, :));
cluster11_data = squeeze(projMatAllSubjC(3, 3, :));
cluster22_data = squeeze(projMatAllSubjC(4, 4, :));

hex1cl1_data = squeeze(projMatAllSubjC(1, 3, :));%project Hex1(big) on Cl
hex2cl2_data = squeeze(projMatAllSubjC(2, 4, :));
hex1cl2_data = squeeze(projMatAllSubjC(1, 4, :));%project on Hex1 (big)
hex2cl1_data = squeeze(projMatAllSubjC(2, 3, :));
cluster1hex1_data = squeeze(projMatAllSubjC(3, 1, :));
cluster2hex2_data = squeeze(projMatAllSubjC(4, 2, :));
cluster1hex2_data = squeeze(projMatAllSubjC(3, 2, :));
cluster2hex1_data = squeeze(projMatAllSubjC(4, 1, :));

hex12_data = squeeze(projMatAllSubjC(1, 2, :));
hex21_data = squeeze(projMatAllSubjC(2, 1, :));
cluster12_data = squeeze(projMatAllSubjC(3, 4, :));
cluster21_data = squeeze(projMatAllSubjC(4, 3, :));


hex_same = hex11_data + hex22_data;
hex_dif = hex12_data + hex21_data;
cluster_same = cluster11_data + cluster22_data;
cluster_dif = cluster12_data + cluster21_data;

hex_cluster_same =  hex1cl1_data + hex2cl2_data; % hex is projected
hex_cluster_dif = hex1cl2_data + hex2cl1_data; % hex is projected
hex1onCluster = hex1cl1_data + hex1cl2_data;
hex2onCluster = hex2cl1_data + hex2cl2_data;
cluster_hex_same = cluster1hex1_data + cluster2hex2_data;% cluster projected 
cluster_hex_dif = cluster1hex2_data + cluster2hex1_data; 
cluster1onHex = cluster1hex1_data + cluster1hex2_data;
cluster2onHex = cluster2hex1_data + cluster2hex2_data;
HexOnCluster = hex1onCluster + hex2onCluster;
clusterOnHex = cluster1onHex + cluster2onHex;

figure(5)
plot(HexOnCluster - clusterOnHex,'.')
hold on
plot(zeros(size(cluster_same)),'k:')
ylabel('HexOnCluster - clusterOnHex')
xlabel('subjects')


figure(2)
subplot(2,2,1)
plot(hex_cluster_same, cluster_hex_same,'.')
hold on
plot(100:102,100:102,'k:')
ylabel('cluster is projected on same stim')
xlabel('hex is projected on same stim')
subplot(2,2,2)
plot(hex_cluster_dif, cluster_hex_dif,'.')
hold on
plot(100:102,100:102,'k:')
ylabel('cluster is projected on dif stim')
xlabel('hex is projected on dif stim')
subplot(2,2,3)
plot(hex1onCluster, hex2onCluster,'.')
hold on
plot(100:102,100:102,'k:')
ylabel('hex2 on cluster both')
xlabel('hex1 on cluster both')
subplot(2,2,4)
plot(cluster1onHex, cluster2onHex,'.')
hold on
plot(100:102,100:102,'k:')
ylabel('cluster2 on hex both')
xlabel('cluster1 on hex both')


figure(1)
subplot(2,2,1)
plot(cluster_same, cluster_dif,'.')
hold on
plot(100:102,100:102,'k:')
ylabel('cluster dif')
xlabel('cluster same')
subplot(2,2,2)
plot(hex_same, hex_dif,'.')
hold on
plot(96:102,96:102,'k:')
ylabel('hex dif')
xlabel('hex same')
subplot(2,2,3)
plot(cluster_same, hex_same,'.')
hold on
plot(100:102,100:102,'k:')
ylabel('hex same')
xlabel('cluster same')
subplot(2,2,4)
plot(cluster_dif, hex_dif,'.')
hold on
plot(100:102,100:102,'k:')
ylabel('hex dif')
xlabel('cluster dif')

figure(11)
subplot(2,2,1)
plot(cluster_same - cluster_dif,'.')
hold on
plot(zeros(size(cluster_same)),'k:')
ylabel('cluster same - cluster dif')
xlabel('subjects')
subplot(2,2,2)
plot(hex_same - hex_dif,'.')
hold on
plot(zeros(size(cluster_same)),'k:')
ylabel('hex same - hex dif')
xlabel('subject')
subplot(2,2,3)
plot(cluster_same - hex_same,'.')
hold on
plot(zeros(size(cluster_same)),'k:')
ylabel('cluster same - hex same')
xlabel('subjects')
subplot(2,2,4)
plot(cluster_dif - hex_dif,'.')
hold on
plot(zeros(size(cluster_same)),'k:')
ylabel('cluster dif - hex dif')
xlabel('subjects')

figure(110)
subplot(2,1,1)
plot(cluster_same - cluster_dif,'.')
hold on
plot(zeros(size(cluster_same)),'k:')
ylabel('cluster same - cluster dif')
xlabel('subjects')
subplot(2,1,2)
plot(cluster_same, cluster_dif,'.')
hold on
plot(100:102,100:102,'k:')
ylabel('cluster dif')
xlabel('cluster same')

figure(20)
subplot(2,2,1)
plot(hex_cluster_same - cluster_hex_same,'.')
hold on
plot(zeros(size(cluster_same)),'k:')
ylabel('hex-cluster_{both on same stim}')
xlabel('subject')
subplot(2,2,2)
plot(hex_cluster_dif- cluster_hex_dif,'.')
hold on
plot(zeros(size(cluster_same)),'k:')
ylabel('hex-cluster_{both on dif stim}')
xlabel('subject')
subplot(2,2,3)
plot(hex1onCluster - hex2onCluster,'.')
hold on
plot(zeros(size(cluster_same)),'k:')
ylabel('hex1 on cluster-hex2 on cluster')
xlabel('subject')
subplot(2,2,4)
plot(cluster1onHex - cluster2onHex,'.')
hold on
plot(zeros(size(cluster_same)),'k:')
ylabel('cluster1 on hex- cluster2 on hex')
xlabel('subject')

allOnCluster = hex1onCluster +hex_cluster_dif + cluster_same + cluster_dif;
allOnHex  = cluster_hex_same + cluster_hex_dif + hex_dif + hex_same;

allOnClusterNoStrDif = hex1onCluster +hex_cluster_dif + cluster_same;
allOnHexNoStrDif  = cluster_hex_same + cluster_hex_dif + hex_same;

figure(40)
subplot(2,1,1)
plot(allOnHex - allOnCluster ,'.')
hold on
plot(zeros(size(cluster_same)),'k:')
ylabel('allOnHex - allOnCluster')
subplot(2,1,2)
plot(allOnHexNoStrDif - allOnClusterNoStrDif ,'.')
hold on
plot(zeros(size(cluster_same)),'k:')
ylabel('allOnHex - allOnCluster_{NoStrDif}')

figure(30)
subplot(2,2,1)
plot(hex11_data, '.')
hold on
plot(projMatMean(1,1)*ones(size(hex11_data)))
title(['hex 11, mean: ', num2str(projMatMean(1, 1))])
subplot(2,2,2)
plot(hex22_data, '.')
hold on
plot(projMatMean(2,2)*ones(size(hex22_data)))
title(['hex 22, mean: ', num2str(projMatMean(2,2))])
subplot(2,2,3)
plot(cluster11_data, '.')
hold on
plot(projMatMean(3,3)*ones(size(cluster11_data)))
title(['cluster 11, mean: ', num2str(projMatMean(3,3))])
subplot(2,2,4)
plot(cluster22_data, '.')
hold on
plot(projMatMean(4,4)*ones(size(cluster22_data)))
title(['cluster 22, mean: ', num2str(projMatMean(4,4))])

figure(21)
subplot(2,2,1)
plot(hex12_data, '.')
hold on
plot(projMatMean(1,2)*ones(size(hex12_data)))
title(['hex 12, mean: ', num2str(projMatMean(1, 2))])
subplot(2,2,2)
plot(hex21_data, '.')
hold on
plot(projMatMean(2,1)*ones(size(hex12_data)))
title(['hex 21, mean: ', num2str(projMatMean(2,1))])
subplot(2,2,3)
plot(cluster12_data, '.')
hold on
plot(projMatMean(3,4)*ones(size(cluster12_data)))
title(['cluster 12, mean: ', num2str(projMatMean(3,4))])
subplot(2,2,4)
plot(cluster22_data, '.')
hold on
plot(projMatMean(4,3)*ones(size(cluster21_data)))
title(['cluster 21, mean: ', num2str(projMatMean(4,3))])

%projMatAll_p_Subj = projMatAllSubj2(:,:,[2:12,14:22,24:27]);
projMatAll_p_Subj = projMatAllSubj2(:,:,[1:12,14:23,25:27]);
figure(101)
imagesc(mean(projMatAll_p_Subj, 3))
colormap gray
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij


figure
hold on
title([titleStr ', t-value'])
[~,~,~,S] = ttest(projMatAllSubj2,50,'dim',3);
projMatT = S.tstat;
imagesc(projMatT);
colormap default
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij

figure(10)
subplot(1,2,1)
projMatMean = squeeze(mean(projMatAllSubj2,3));
imagesc(projMatMean);
title([titleStr ', mean'])
colormap default
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij

subplot(1,2,2)
[~,~,~,S] = ttest(projMatAllSubj2,50,'dim',3);
projMatT = S.tstat;
imagesc(projMatT);
title([titleStr ', t-value'])
colormap default
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij

function plot_warm(num_figure, data, ylabel_plot)
    % Calculate mean and standard deviation
    meanValue = mean(data);
    stdValue = std(data);

    subplot(2,2,num_figure)
    % Create a swarm plot
    swarmchart(ones(size(data)),data, 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none','SizeData', 5);
    
    hold on;
    
    % Plot mean value
    plot(1, meanValue, 'r*', 'MarkerSize', 10);
    
    % Plot standard deviation range
    xRange = [0.8, 1.2]; % x-axis range for std indicator
    yRange = [meanValue - stdValue, meanValue + stdValue]; % y-axis range for std indicator
    errorbar(1,meanValue,stdValue/sqrt(length(data)));
    xlim([-0.5,2.5])
    plot([-0.5,2.5],[0,0],'k--')
    hold off;
    
    % Adjust plot labels and legend
    xticks([]);
    yticks(0);
    ylabel(ylabel_plot)
end