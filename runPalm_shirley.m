function runPalm_shirley(conName,nPerm, maskName, output_dir)
% conName: name of contrast to run

root = '/data/holly-host/smark';
spm_path = '/data/holly-host/smark/spm';
addpath(spm_path)
data_path0 = '/data/holly-host/smark/fmri_sub_preproc_dir';
data_path = fullfile(data_path0, output_dir, 'fsl_norm_new_names');

% make sure palm is installed an in path.
palmDir = fullfile(root, 'palm-alpha116');
addpath(genpath(palmDir));

inputFile = fullfile(data_path,'groupStats','palm','stackedInputs',['con_' conName '.nii']);
disp(inputFile);
outputDir = fullfile(data_path,'groupStats','palm');
outputName = [conName '_nPerm_' nPerm '_' maskName];
maskFile = fullfile(data_path0,'ROI_masks',[maskName '.nii']);

palmStr = ['palm -i ' inputFile ' -T -n ' nPerm ' -o ' fullfile(outputDir,outputName) ' -ise -save1-p -m ' maskFile];
eval(palmStr)