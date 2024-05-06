close all
clear all
clc

%% compute choose connecting node for wrong answer
%% show that this number is significantly larger than 0.5

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
root = 'C:\Users\User\Documents\fMRI_EXP\Alon\';
behvioural_analysis_codes = 'C:\Users\User\Documents\psychExpLaptopOld';

addpath(genpath(behvioural_analysis_codes));

dataDir = fullfile(root,"beh","training");

A1 = createAll2Allcluster(35,5);%
Distm1 = calDistM(A1);

A2 = createAll2Allcluster(42,6);%
Distm2 = calDistM(A2);

end0 = 16;%1;%end trial of previous structure, or in hex exclude the first trial as the number of step to the target is given
end1 = 24;%8;%end trial day1
end2 = 32;%16;%end trial day2

ConNode = [1,7,8,14,15,21,22,28,29,35,42];

load(fullfile(dataDir,'allTrainingData.mat'),'D');

nSub = length(D.middle.hex_day1_ratioCorrect);

lenn= nSub;

mT = readtable(fullfile(dataDir, 'navigTable.xlsx'));

%trials to check:
Tcheck1 = (end0+1):end1;
Tcheck2 = (end1+1):end2;

sConCsCcor = zeros(lenn,8);
nCCsCcor = zeros(lenn,8);
sConCsCncor = zeros(lenn,8);
nCCsCncor = zeros(lenn,8);

sConCsCcor2 = zeros(lenn,8);
nCCsCcor2 = zeros(lenn,8);
sConCsCncor2 = zeros(lenn,8);
nCCsCncor2 = zeros(lenn,8);

for n=1:lenn
    if n<10
        subName = ['sub-0', num2str(n)];
    else
        subName = ['sub-', num2str(n)];
    end
    namesT = mT.Var1;
    namesT = lower(namesT);
    subVT = strcmp(subName,namesT);
    trialT = mT(subVT,:).Var2;

    targ =  mT(subVT,:).Var5;
    isCorT = mT(subVT,:).Var11;
    im2 = mT(subVT,:).Var9;
    im1 = mT(subVT,:).Var10;
%    [im1,im2] = makeMat(imL,imR);
    choice = mT(subVT,:).Var7;

   [sConCcor,sConCncor,nCcor,nCncor] = isChcConTaskCorNoCorNew_fMRI(Distm1, im1, im2, targ, choice, trialT, ConNode, end0, 20, isCorT);
    sConCsCcor(n,1:4) = sConCcor;
    nCCsCcor(n,1:4) = nCcor;
    sConCsCncor(n,1:4) = sConCncor;
    nCCsCncor(n,1:4) = nCncor;
    [sConCcor,sConCncor,nCcor,nCncor] = isChcConTaskCorNoCorNew_fMRI(Distm2, im1, im2, targ, choice, trialT, ConNode, 20, end1, isCorT);
    sConCsCcor(n,5:8) = sConCcor;
    nCCsCcor(n,5:8) = nCcor;
    sConCsCncor(n,5:8) = sConCncor;
    nCCsCncor(n,5:8) = nCncor;
    [sConCcor,sConCncor,nCcor,nCncor]  = isChcConTaskCorNoCorNew_fMRI(Distm1, im1, im2, targ, choice, trialT, ConNode, end1, 28, isCorT);
    sConCsCcor2(n,1:4) = sConCcor;
    nCCsCcor2(n,1:4) = nCcor;
    sConCsCncor2(n,1:4) = sConCncor;
    nCCsCncor2(n,1:4) = nCncor;
    [sConCcor,sConCncor,nCcor,nCncor]  = isChcConTaskCorNoCorNew_fMRI(Distm2, im1, im2, targ, choice, trialT, ConNode, 28, end2, isCorT);
    sConCsCcor2(n,5:8) = sConCcor;
    nCCsCcor2(n,5:8) = nCcor;
    sConCsCncor2(n,5:8) = sConCncor;
    nCCsCncor2(n,5:8) = nCncor;
end

nCCsCncor2a = nCCsCncor2;
nCCsCcor2a = nCCsCcor2;
for n=1:nSub
    nCCsCncor2a(n,nCCsCncor2(n,:)==0) = 1;
    nCCsCcor2a(n,nCCsCcor2(n,:)==0) = 1;
end
pConCsCncor1 =  sum(sConCsCncor,2)./sum(nCCsCncor,2);
pConCsCncor2 = sum(sConCsCncor2,2)./sum(nCCsCncor2a ,2);
pConCsCcor1 =  sum(sConCsCcor,2)./sum(nCCsCcor,2);
pConCsCcor2 = sum(sConCsCcor2,2)./sum(nCCsCcor2a, 2);

save(fullfile(dataDir,'pConCsCncor1.mat'),'pConCsCncor1')
save(fullfile(dataDir,'pConCsCncor2.mat'),'pConCsCncor2')
save(fullfile(dataDir,'pConCsCcor1.mat'),'pConCsCcor1')
save(fullfile(dataDir,'pConCsCcor2.mat'),'pConCsCcor2')

[h1, p1, ~,stats] = ttest(pConCsCncor1, 0.5, 'Tail', 'right')
[h2, p2] = ttest(pConCsCncor2, 0.5, 'Tail', 'right')
[h12, p12, ~, stats12] = ttest2(pConCsCncor1, pConCsCncor2)

figure(14)
errorbar(1:2,[mean(pConCsCncor1), mean(pConCsCncor2)], [std(pConCsCncor1) / sqrt(length(pConCsCncor1)), std(pConCsCncor2)/ sqrt(length(pConCsCncor2))],'.k')
hold on
scatter(ones(size(pConCsCncor1)),pConCsCncor1,3,'fill','k')
hold on
scatter(2*ones(size(pConCsCncor2)),pConCsCncor2,3,'fill','k')
xlim([0.5,2.5])
ylim([0.2,0.95])


