function createSl(root,sub,glm,nVoxelSl)

spm_path = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(spm_path)
util_path = fullfile(root,'code','util');
addpath(genpath(util_path))

cd(fullfile(root,'masks',sub))

%% read ROI/brain mask:
wholeBrainMask = fullfile(root,'masks',sub,'brain.nii');
ROI         = spm_vol(wholeBrainMask);
ROI.data    = spm_read_vols(ROI);

nanMask = fullfile(root,'masks',sub,[glm '_noNanNoInf.nii']);

% nameNanMask = [based,'\fmri_sub_preproc_dir\masks\nativeNotNanMaskSub',num2str(vnpar(sb)),'300419betas.nii'];
nanMask     = spm_vol(nanMask);
nanMask.data= spm_read_vols(nanMask);

L = rsa.defineSearchlight_volume(ROI,nanMask,'sphere',[70 nVoxelSl]); % maximum radius 70 to find nVoxelL voxels

slDir = fullfile(root,'sl',sub);
mkdir(slDir);
save([slDir,'/L',num2str(nVoxelSl),'.mat'],'L');
disp(['subject ',sub,' done searchlight'])
