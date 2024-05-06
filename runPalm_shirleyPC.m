function runPalm_shirleyPC(conName,nPerm, maskName, output_dir)
% conName: name of contrast to run

root = 'C:\Users\User\Documents\fMRI_EXP\';
spm_path =  "C:\Users\User\Documents\spm12";
addpath(spm_path)
data_path = fullfile(root, 'Alon');
data_path = fullfile(data_path, output_dir, 'fsl_normalization');

% make sure palm is installed an in path.
palmDir = fullfile(root, 'palm-alpha119');
addpath(genpath(palmDir));

inputFile = fullfile(data_path,'groupStats','palm','stackedInputs',['con_' conName '.nii']);
disp(inputFile);
outputDir = fullfile(data_path,'groupStats','palm');
outputName = [conName '_nPerm_' nPerm '_' maskName];
maskFile = fullfile(root,'roi_masks_fsl',[maskName '.nii']);

palmStr = ['palm -i ' inputFile ' -T -n ' nPerm ' -o ' fullfile(outputDir,outputName) ' -ise -save1-p -m ' maskFile];
eval(palmStr)