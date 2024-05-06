% group stats pipline

% normalized using fsl and gunzip
% smooth using spm
% calculate group stats

% Inputs:
% slName - the search light name
% output_dir - where the results that need to be norm and smoothed are located
%              it should be after the data_path
close all
clear all
clc

slName = 'L100_projMat';
output_dir = 'subspaceGener'; 
doNorm = false; 
doSmooth = false;
make_group_stat= false;
save_contrast_palm = false;%true;
stack4palm = false;%true;
run_palm = true;

status = 0;
nSm = 6;
nFiles = 16;

home_dir = 'C:\Users\User\Documents\fMRI_EXP\';
spm_path =  "C:\Users\User\Documents\spm12";
my_code_path = fullfile(home_dir, 'afterSearchLightCodes');
data_path = fullfile(home_dir, 'Alon');
data_path = fullfile(data_path, output_dir, 'fsl_normalization');
util_path = fullfile(data_path, 'from_git', 'util');
addpath(spm_path)
addpath(my_code_path)
addpath(fullfile(home_dir, 'Alon'))
subjects = cell(1,28);
for iSub = 1:length(subjects)
    subjects{iSub}= ['sub-' num2str(iSub,'%02.f')];
end

if doNorm
    %%% fsl normalization and gunzip:
    disp('doing normalization')
    scriptPath = 'normalize_results_fsl_zappa_matlab.sh';
    % Check if the script exists
    if exist(scriptPath, 'file') ~= 2
        error('Noramlization script file does not exist.');
    end

    % Execute the shell script
    [status, result] = system(['sh ' scriptPath, ' ' output_dir ' ' slName]);
end
% Check if execution was successful
if status ~= 0
    error('Error executing the shell script:\n%s', result);
else
    disp('Shell script executed successfully.');

    %%% smooth: %%%%%
    if doSmooth
        disp('doing smoothing')
        for s= 1:28
            %s = 1;
            if s<10
                sub = ['sub-0',num2str(s)];
            else
                sub = ['sub-',num2str(s)];
            end
            searchlight_dir = fullfile(data_path, 'fsl_norm_new_names', sub);
            %slName = 'L100_projMat';

            searchLightResFiles = cell(nFiles,1);
            for f=1:nFiles
                disp(fullfile(searchlight_dir,['standard_',slName,num2str(f),'__brain.nii']))
                searchLightResFiles{f} = fullfile(searchlight_dir,['standard_',slName,num2str(f),'__brain.nii']);
            end
            disp('loaded searchlight files')
            clear matlabbatch

            matlabbatch{1}.spm.spatial.smooth.data =  searchLightResFiles;%cellstr(spm_select('FPList', searchlight_dir, '^standard_L100_projMat*.nii'));
            matlabbatch{1}.spm.spatial.smooth.fwhm   = [nSm nSm nSm];
            matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
            matlabbatch{1}.spm.spatial.smooth.im     = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = ['smth',num2str(nSm)];

            spm_jobman('run', matlabbatch)
        end
    end
    if make_group_stat
        % The directory fullfile(output_dir, 'fsl_norm_new_names') should
        % be created at the above stages
        disp('doing group stats')
        groupLevelTTest_cluster(fullfile(output_dir, 'fsl_norm_new_names'));
    end
    if save_contrast_palm
        disp('saving contrast for palm')
        for n_sub = 1:28
            calcProjMatContrastsSubj_forPalm_ShirleyPC(slName, n_sub, spm_path, data_path, util_path);
            disp(['done subject ', num2str(n_sub)])
        end
    end
    if stack4palm
        disp('stack contrast for palm')
        stackSubjectsContrasts_shirley(output_dir,subjects)
    end
    if run_palm
       % cd(home_dir)
        disp('do palm')
        nPerm = '10000';
        %conName = 'SameStr_allOthers';
        conName = 'hexOnHexMinusHexOnCl';%'visual_sameStructSameStimMinusSameStructDiffStim';%'hexOnHexMinusHexOnCl';
         maskName =  'standard_mask';
       % maskName = 'juelich_EC10';%'LOC'; %'harvardOxford_iLOC50';%'LOC'; %'juelich_EC10';
        %maskName = 'harvardoxfrod_frontal_medial_subcallosal_cortex_maskT25';
        runPalm_shirleyPC(conName,nPerm,maskName, output_dir)
    end
end


