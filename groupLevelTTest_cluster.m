function groupLevelTTest_cluster(output_dir)

nFiles = 16;
rsa_tool_path = 'rsatoolbox';

spm_path = '/data/holly-host/smark/spm';
data_path = fullfile('/data/holly-host/smark/fmri_sub_preproc_dir',output_dir);
addpath(genpath(fullfile(data_path,rsa_tool_path)))
addpath(spm_path)
ana_name = 'standard_L100_FSLcleaned280920b';

group_level_output_path = fullfile(data_path, 'groupStats');
if ~exist(group_level_output_path)
    mkdir(group_level_output_path);
end
num_subjects = 28;
subjects = cell(num_subjects);
for s=1:28
    if s<10
        subjects{s} = ['sub-0',num2str(s)];
    else
        subjects{s} = ['sub-',num2str(s)];
    end
end

niiDims = [91,109,91]; % dims of standard brain
projMat = zeros(niiDims(1),niiDims(2),niiDims(3),nFiles,length(subjects)); 
for iSub=1:length(subjects)
    sub = subjects{iSub};
    nameDir = fullfile(data_path,sub);
    projMatFiles = cell(nFiles,1);
    for f=1:nFiles
        projMatFiles{f} = fullfile(nameDir,['smth6',ana_name,num2str(f),'__brain.nii']);
        % V           = spm_vol(projMatFiles{f}); % heatder info      
        % [projMat_currElement,~] = spm_read_vols(V,0); 
        %projMat1 = read_avw(projMatFiles{f});
        %projMat(:,:,:,f,iSub) = projMat1;
        %projMat(:,:,:,f,iSub) = projMat_currElement;
        projMat(:,:,:,f,iSub) = niftiread( projMatFiles{f});
    end
end
V           = spm_vol(projMatFiles{f});

nVox = prod(niiDims);

% projection of Hex on hex: 1: HlHl 2: HsHl 5: HlHs 6: HsHs
projHex = squeeze(projMat(:,:,:,1,:) + projMat(:,:,:,2,:) + projMat(:,:,:,5,:) + projMat(:,:,:,6,:));
% projections of Hex on Cluster: 9: HlCl 10: HsCl 13: HlCs 14 HsCs
projHexCl = squeeze(projMat(:,:,:,9,:) + projMat(:,:,:,10,:) + projMat(:,:,:,13,:) + projMat(:,:,:,14,:) );
% projections of Cluster on Hex:
projClHex = squeeze(projMat(:,:,:,3,:) + projMat(:,:,:,4,:) + projMat(:,:,:,7,:) + projMat(:,:,:,8,:));
% projections of Cluster on Cluster:
projClCl = squeeze(projMat(:,:,:,11,:) + projMat(:,:,:,12,:) + projMat(:,:,:,15,:) + projMat(:,:,:,16,:));

Hex_same_minus_same_stim = squeeze(projMat(:,:,:,1,:) + projMat(:,:,:,6,:) - projMat(:,:,:,10,:) - projMat(:,:,:,13,:));
Hex_str_minus_nothing = squeeze(projMat(:,:,:,2,:) + projMat(:,:,:,5,:) - projMat(:,:,:,9,:) - projMat(:,:,:,14,:) );

Cluster_same_minus_same_stim = squeeze(projMat(:,:,:,11,:) + projMat(:,:,:,16,:) - projMat(:,:,:,4,:) - projMat(:,:,:,7,:));
Cluster_str_minus_nothing = squeeze( projMat(:,:,:,12,:) + projMat(:,:,:,15,:) - projMat(:,:,:,3,:) - projMat(:,:,:,8,:));

Cluster_same_minus_same_stim_linVox = reshape(Cluster_same_minus_same_stim,[nVox,length(subjects)]);
Cluster_str_minus_nothing_linVox = reshape(Cluster_str_minus_nothing,[nVox,length(subjects)]);
% tStats_linVox = zeros(nVox,1);
% tStats_same_linVox = zeros(nVox,1);
% for iVox = 1:nVox
%     [~,~,~,stats] = ttest(Cluster_same_minus_same_stim_linVox(iVox,:),0,'tail','right');
%     tStats_same_linVox(iVox) = stats.tstat;
%     [~,~,~,stats] = ttest(Cluster_str_minus_nothing_linVox(iVox,:),0,'tail','right');
%     tStats_linVox(iVox) = stats.tstat;
% end
% tStats_same = reshape(tStats_same_linVox,niiDims);
% save4Dnii(tStats_same,V,fullfile(data_path,'groupStats','tStat_ClSameMinusSameStim.nii'))
% tStats = reshape(tStats_linVox,niiDims);
% save4Dnii(tStats,V,fullfile(group_level_output_path,'tStat_ClStrMinusNothing.nii'))

% % contrasts:
projHex_contrast = projHex - projHexCl; 
projHex_ClHex_contrast = projHex - projClHex; 
projClCl_ClHex_contrast = projClCl - projClHex; 
projSameStr_allOthers = projClCl_ClHex_contrast + projHex_contrast;
clear proHex proHexCl projClHex projClCl

%% Calculate stim set ID contrast
% HlHl, HsHs, ClCl, CsCs
projSameStructSameStim = squeeze(projMat(:,:,:,1,:) + projMat(:,:,:,6,:) + projMat(:,:,:,11,:) + projMat(:,:,:,16,:));
% HsHl, HlHs, CsCl, ClCs
projSameStructDiffStim = squeeze(projMat(:,:,:,2,:) + projMat(:,:,:,5,:) + projMat(:,:,:,12,:) + projMat(:,:,:,15,:));
projStimSet_contrast = projSameStructSameStim - projSameStructDiffStim;
clear projMat projSameStructSameStim projSameStructDiffStim

%% stats
nVox = prod(niiDims);
projHex_contrast_linVox = reshape(projHex_contrast,[nVox,length(subjects)]);
tStats_linVox = zeros(nVox,1);
for iVox = 1:nVox
    [~,~,~,stats] = ttest(projHex_contrast_linVox(iVox,:),0,'tail','right');
    tStats_linVox(iVox) = stats.tstat;
end
tStats = reshape(tStats_linVox,niiDims);
save4Dnii(tStats,V,fullfile(group_level_output_path,'tStat_hexOnHexMinusHexOnCluster.nii'))


%nVox = prod(niiDims);
projClCl_ClHex_contrast_linVox = reshape(projClCl_ClHex_contrast,[nVox,length(subjects)]);
tStats_linVox = zeros(nVox,1);
for iVox = 1:nVox
    [~,~,~,stats] = ttest(projClCl_ClHex_contrast_linVox(iVox,:),0,'tail','right');
    tStats_linVox(iVox) = stats.tstat;
end
tStats = reshape(tStats_linVox,niiDims);
save4Dnii(tStats,V,fullfile(group_level_output_path,'tStat_ClOnClMinusClusterOnHex.nii'))

nVox = prod(niiDims);
projSameStr_allOthers_linVox = reshape(projSameStr_allOthers,[nVox,length(subjects)]);
tStats_linVox = zeros(nVox,1);
for iVox = 1:nVox
    [~,~,~,stats] = ttest(projSameStr_allOthers_linVox(iVox,:),0,'tail','right');
    tStats_linVox(iVox) = stats.tstat;
end
tStats = reshape(tStats_linVox,niiDims);
save4Dnii(tStats,V,fullfile(group_level_output_path,'tStat_projSameStr_allOthers.nii'))

Hex_same_minus_same_stim_linVox = reshape(Hex_same_minus_same_stim,[nVox,length(subjects)]);
Hex_str_minus_nothing_linVox = reshape(Hex_str_minus_nothing,[nVox,length(subjects)]);
tStats_linVox = zeros(nVox,1);
tStats_same_linVox = zeros(nVox,1);
for iVox = 1:nVox
    [~,~,~,stats] = ttest(Hex_same_minus_same_stim_linVox(iVox,:),0,'tail','right');
    tStats_same_linVox(iVox) = stats.tstat;
    [~,~,~,stats] = ttest(Hex_str_minus_nothing_linVox(iVox,:),0,'tail','right');
    tStats_linVox(iVox) = stats.tstat;
end
tStats_same = reshape(tStats_same_linVox,niiDims);
save4Dnii(tStats_same,V,fullfile(data_path,'groupStats','tStat_hexSameMinusSameStim.nii'))
tStats = reshape(tStats_linVox,niiDims);
save4Dnii(tStats,V,fullfile(group_level_output_path,'tStat_hexStrMinusNothing.nii'))

% nVox = prod(niiDims);
projHex_ClHex_contrast_linVox = reshape(projHex_ClHex_contrast,[nVox,length(subjects)]);
tStats_linVox = zeros(nVox,1);
for iVox = 1:nVox
    [~,~,~,stats] = ttest(projHex_ClHex_contrast_linVox(iVox,:),0,'tail','right');
    tStats_linVox(iVox) = stats.tstat;
end
tStats = reshape(tStats_linVox,niiDims);
save4Dnii(tStats,V,fullfile(group_level_output_path,'tStat_hexOnHexMinusClusterOnHex.nii'))

projStimSet_contrast_linVox = reshape(projStimSet_contrast,[nVox,length(subjects)]);
tStats_linVox = zeros(nVox,1);
for iVox = 1:nVox
    [~,~,~,stats] = ttest(projStimSet_contrast_linVox(iVox,:),0,'tail','right');
    tStats_linVox(iVox) = stats.tstat;
end
tStats = reshape(tStats_linVox,niiDims);
save4Dnii(tStats,V,fullfile(group_level_output_path,'tStat_visual_sameStructSameStimMinusSameStructDiffStim.nii'))

