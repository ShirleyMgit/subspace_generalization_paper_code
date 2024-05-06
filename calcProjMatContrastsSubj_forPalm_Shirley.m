function calcProjMatContrastsSubj_forPalm_Shirley(matFname, analysis_name , n_sub)
% save single subject contrusts of Tim's PCA projections analysis
spm_path = '/data/holly-host/smark/spm';
data_path = '/data/holly-host/smark/fmri_sub_preproc_dir';
addpath(genpath(spm_path));
util_path = fullfile(data_path,'alon_util');
addpath(genpath(util_path));

nFiles = 16;
smoothKernel = 6;

projMat = zeros(91,109,91,nFiles); % dims of standard brain

if n_sub < 10
    sub = ['sub-0',num2str(n_sub)];
else
    sub = ['sub-',num2str(n_sub)];
end

results_dir = fullfile(data_path, analysis_name,'fsl_norm_new_names' ,sub);
projMatFiles = cell(nFiles,1);
for f=1:nFiles
    %projMatFiles{f} = [results_dir,'/smth',num2str(smoothKernel) 'standard_' matFname num2str(f),'__brain.nii'];
    projMatFiles{f} = [results_dir,'/smth',num2str(smoothKernel) 'standard_' matFname num2str(f),'.nii'];
    V           = spm_vol(projMatFiles{f}); % header info
    [projMat_currElement,~] = spm_read_vols(V,0);
    projMat(:,:,:,f) = projMat_currElement;
end

subjContrastsDir = fullfile(results_dir,'contrasts');
if ~exist(subjContrastsDir,'dir')
    mkdir(subjContrastsDir);
end

%% Calculate main Hex contrast - same labels as in the paper (e.g.
% HsCs denotes the element of the matrix corresponding to activity from the
% small hexagonal graph projected on eigenvectors calculated from the small
% (same stimuli-set) community-structure graph.
% 1:4 : projections on eigenvectors from (3 averaged runs of) Hl.
% 5:8 : projections on eigenvectors from (3 averaged runs of) Hs.
% 9:12 : projections on eigenvectors from (3 averaged runs of) Cl.
% 13:16 : projections on eigenvectors from (3 averaged runs of) Cs.

% projection of Hex on hex: 1: HlHl 2: HsHl 5: HlHs 6: HsHs
proHex = squeeze(projMat(:,:,:,1) + projMat(:,:,:,2) + projMat(:,:,:,5) + projMat(:,:,:,6));
% projections of Hex on Cluster: 9: HlCl 10: HsCl 13: HlCs 14 HsCs
proHexCl = squeeze(projMat(:,:,:,9) + projMat(:,:,:,10) + projMat(:,:,:,13) + projMat(:,:,:,14) );
projHex_contrast = proHex - proHexCl;

% projections of Cluster on Hex:
projClHex = squeeze(projMat(:,:,:,3) + projMat(:,:,:,4) + projMat(:,:,:,7) + projMat(:,:,:,8));
% projections of Cluster on Cluster:
projClCl = squeeze(projMat(:,:,:,11) + projMat(:,:,:,12) + projMat(:,:,:,15) + projMat(:,:,:,16));
projClCl_ClHex_contrast = projClCl - projClHex; 

projSameStr_allOthers = projClCl_ClHex_contrast + projHex_contrast;

pathContrast_All = fullfile(subjContrastsDir,'con_SameStr_allOthers.nii');
save4Dnii(projSameStr_allOthers,V,pathContrast_All)

pathContrast_Cl = fullfile(subjContrastsDir,'con_ClOnClMinusClustOnHex.nii');
save4Dnii(projClCl_ClHex_contrast,V,pathContrast_Cl)

pathContrast_Hex = fullfile(subjContrastsDir,'con_hexOnHexMinusHexOnCl.nii');
save4Dnii(projHex_contrast,V,pathContrast_Hex)


%% Calculate stim set ID contrast
% HlHl, HsHs, ClCl, CsCs
projSameStructSameStim = squeeze(projMat(:,:,:,1) + projMat(:,:,:,6) + projMat(:,:,:,11) + projMat(:,:,:,16));
% HsHl, HlHs, CsCl, ClCs
projSameStructDiffStim = squeeze(projMat(:,:,:,2) + projMat(:,:,:,5) + projMat(:,:,:,12) + projMat(:,:,:,15));

projStimSet_contrast = projSameStructSameStim - projSameStructDiffStim;
pathContrast = fullfile(subjContrastsDir,'con_visual_sameStructSameStimMinusSameStructDiffStim.nii');
save4Dnii(projStimSet_contrast,V,pathContrast)






