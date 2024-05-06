function [stepPerTrialWithSkips,stepsPerTrial] = calc_nSteps_navig(block,choice,initDist)
% block, choice and initDist are all nChoices x 1 vectors, where each row
% is a single choice made  the subject in any trial (a path from start to
% target image).  
% blocks:  the block number of the trial (1-32, depends on the training
% day)
% choice: 1 and 2 for image options, 0 for pressing enter to sample more
% options 
% initDist: the initial distance between start and target images in the
% current trial. 

blockInds = unique(block);

stepPerTrialWithSkips = zeros(length(blockInds),3);
stepsPerTrial = zeros(length(blockInds),3);

for t=1:length(blockInds)
    for d = 1:3
        currTrialLogical = (block==(blockInds(1)-1+t)).*(initDist==d+1); % 1s for all steps of current trial
        stepPerTrialWithSkips(t,d) = sum(currTrialLogical); % number of steps for this trial 
        currTrialLogical = currTrialLogical.*(choice>0); % 
        stepsPerTrial(t,d) = sum(currTrialLogical);
    end
end
