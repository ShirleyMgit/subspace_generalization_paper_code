close all
clear
clc

% dataDir = '/Volumes/Scratch_abaram/shirley/alon/beh/training';
dataDir = '/vols/Scratch/abaram/shirley/alon/beh/training';

load(fullfile(dataDir,'allTrainingData.mat'),'D');

nSub = length(D.middle.hex_day1_ratioCorrect);

figure

%% is in middle task

% plot all datapoints
subplot(2,3,1)
hold on
title('Is in middle?')
% Plot individual data points
scatter(ones(nSub,1),100*D.middle.hex_day1_ratioCorrect,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(2*ones(nSub,1),100*D.middle.hex_day2_ratioCorrect,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(3*ones(nSub,1),100*D.middle.clus_day3_ratioCorrect,3,'filled','b','MarkerFaceAlpha',0.5)
scatter(4*ones(nSub,1),100*D.middle.clus_day4_ratioCorrect,3,'filled','b','MarkerFaceAlpha',0.5)

% plot means and error bars
wBar = 0.15;
plot([1+0.2-wBar, 1+0.2+wBar], 100*[mean(D.middle.hex_day1_ratioCorrect),mean(D.middle.hex_day1_ratioCorrect)],'k','lineWidth',2)
plot([2-0.2-wBar, 2-0.2+wBar], 100*[mean(D.middle.hex_day2_ratioCorrect),mean(D.middle.hex_day2_ratioCorrect)],'k','lineWidth',2)
plot([3+0.2-wBar, 3+0.2+wBar], 100*[mean(D.middle.clus_day3_ratioCorrect),mean(D.middle.clus_day3_ratioCorrect)],'k','lineWidth',2)
plot([4-0.2-wBar, 4-0.2+wBar], 100*[mean(D.middle.clus_day4_ratioCorrect),mean(D.middle.clus_day4_ratioCorrect)],'k','lineWidth',2)

errorbar(1+0.2,100*mean(D.middle.hex_day1_ratioCorrect),100*std(D.middle.hex_day1_ratioCorrect)/sqrt(nSub),'k')
errorbar(2-0.2,100*mean(D.middle.hex_day2_ratioCorrect),100*std(D.middle.hex_day2_ratioCorrect)/sqrt(nSub),'k')
errorbar(3+0.2,100*mean(D.middle.clus_day3_ratioCorrect),100*std(D.middle.clus_day3_ratioCorrect)/sqrt(nSub),'k')
errorbar(4-0.2,100*mean(D.middle.clus_day4_ratioCorrect),100*std(D.middle.clus_day4_ratioCorrect)/sqrt(nSub),'k')

% Plot lines connecting the data points
for iSub=1:nSub
    plot([1,2],100*[D.middle.hex_day1_ratioCorrect(iSub),D.middle.hex_day2_ratioCorrect(iSub)],'color',[1 0 0 0.2])
    plot([3,4],100*[D.middle.clus_day3_ratioCorrect(iSub),D.middle.clus_day4_ratioCorrect(iSub)],'color',[0 0 1 0.2])
end
xlim([0,5])
labels = {'hex day1','hex day2','comm day3', 'comm day4'};
set(gca,'FontSize',14,'XTick',1:4,'XTickLabels',labels,'XTickLabelRotation', 45)
ylabel('% correct')


% stats for whether answers oneach day were different from chanve
[~, pMiddleHexDay1,~,s] = ttest(D.middle.hex_day1_ratioCorrect,0.5,'Tail','right');
fprintf('Middle, day 1 p = %d, t = %d \n',pMiddleHexDay1,s.tstat)
[~, pMiddleHexDay2,~,s] = ttest(D.middle.hex_day2_ratioCorrect,0.5,'Tail','right');
fprintf('Middle, day 2 p = %d, t = %d \n',pMiddleHexDay2,s.tstat)
[~, pMiddleClusDay3,~,s] = ttest(D.middle.clus_day3_ratioCorrect,0.5,'Tail','right');
fprintf('Middle, day 3 p = %d, t = %d \n',pMiddleClusDay3,s.tstat)
[~, pMiddleClusDay4,~,s] = ttest(D.middle.clus_day4_ratioCorrect,0.5,'Tail','right');
fprintf('Middle, day 4 p = %d, t = %d \n',pMiddleClusDay4,s.tstat)

% stats for improvement beween days
[~, pMiddleHexDays,~,s] = ttest2(D.middle.hex_day2_ratioCorrect,D.middle.hex_day1_ratioCorrect,'Tail','right');
fprintf('Middle, days 1-2 p = %d, t = %d \n',pMiddleHexDays,s.tstat)
[~, pMiddleClusDays,~,s] = ttest2(D.middle.clus_day4_ratioCorrect,D.middle.clus_day3_ratioCorrect,'Tail','right');
fprintf('Middle, days 3-4 p = %d, t = %d \n',pMiddleClusDays,s.tstat)

%% piles task
subplot(2,3,2)
hold on
title('Is this image next?')
% Plot individual data points
scatter(ones(nSub,1),100*D.piles.hex_day1_ratioCorrect,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(2*ones(nSub,1),100*D.piles.hex_day2_ratioCorrect,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(3*ones(nSub,1),100*D.piles.clus_day3_ratioCorrect,3,'filled','b','MarkerFaceAlpha',0.5)
scatter(4*ones(nSub,1),100*D.piles.clus_day4_ratioCorrect,3,'filled','b','MarkerFaceAlpha',0.5)

% plot means and SEM error bars
wBar = 0.15;
plot([1+0.2-wBar, 1+0.2+wBar], 100*[mean(D.piles.hex_day1_ratioCorrect),mean(D.piles.hex_day1_ratioCorrect)],'k','lineWidth',2)
plot([2-0.2-wBar, 2-0.2+wBar], 100*[mean(D.piles.hex_day2_ratioCorrect),mean(D.piles.hex_day2_ratioCorrect)],'k','lineWidth',2)
plot([3+0.2-wBar, 3+0.2+wBar], 100*[mean(D.piles.clus_day3_ratioCorrect),mean(D.piles.clus_day3_ratioCorrect)],'k','lineWidth',2)
plot([4-0.2-wBar, 4-0.2+wBar], 100*[mean(D.piles.clus_day4_ratioCorrect),mean(D.piles.clus_day4_ratioCorrect)],'k','lineWidth',2)

errorbar(1+0.2,100*mean(D.piles.hex_day1_ratioCorrect),100*std(D.piles.hex_day1_ratioCorrect)/sqrt(nSub),'k')
errorbar(2-0.2,100*mean(D.piles.hex_day2_ratioCorrect),100*std(D.piles.hex_day2_ratioCorrect)/sqrt(nSub),'k')
errorbar(3+0.2,100*mean(D.piles.clus_day3_ratioCorrect),100*std(D.piles.clus_day3_ratioCorrect)/sqrt(nSub),'k')
errorbar(4-0.2,100*mean(D.piles.clus_day4_ratioCorrect),100*std(D.piles.clus_day4_ratioCorrect)/sqrt(nSub),'k')

% Plot lines connecting the data points
for iSub=1:nSub
    plot([1,2],100*[D.piles.hex_day1_ratioCorrect(iSub),D.piles.hex_day2_ratioCorrect(iSub)],'color',[1 0 0 0.2])
    plot([3,4],100*[D.piles.clus_day3_ratioCorrect(iSub),D.piles.clus_day4_ratioCorrect(iSub)],'color',[0 0 1 0.2])
end

xlim([0,5]);
ylabel('% correct')
labels = {'hex day1','hex day2','comm day3', 'comm day4'};
set(gca,'FontSize',14,'XTick',1:4,'XTickLabels',labels,'XTickLabelRotation', 45)



% stats for whether answers oneach day were different from chanve
[~, pPilesHexDay1,~,s] = ttest(D.piles.hex_day1_ratioCorrect,0.333,'Tail','right');
fprintf('Piles, day 1 p = %d, t = %d \n',pPilesHexDay1,s.tstat)
[~, pPilesHexDay2,~,s] = ttest(D.piles.hex_day2_ratioCorrect,0.333,'Tail','right');
fprintf('Piles, day 2 p = %d, t = %d \n',pPilesHexDay2,s.tstat)
[~, pPilesClusDay3,~,s] = ttest(D.piles.clus_day3_ratioCorrect,0.333,'Tail','right');
fprintf('Piles, day 3 p = %d, t = %d \n',pPilesClusDay3,s.tstat)
[~, pPilesClusDay4,~,s] = ttest(D.piles.clus_day4_ratioCorrect,0.333,'Tail','right');
fprintf('Piles, day 4 p = %d, t = %d \n',pPilesClusDay4,s.tstat)

% stats for improvement beween days
[~, pPilesHexDays,~,s] = ttest2(D.piles.hex_day2_ratioCorrect,D.piles.hex_day1_ratioCorrect,'Tail','right');
fprintf('Piles, days 1-2 p = %d, t = %d \n',pPilesHexDays,s.tstat)
[~, pPilesClusDays,~,s] = ttest2(D.piles.clus_day4_ratioCorrect,D.piles.clus_day3_ratioCorrect,'Tail','right');
fprintf('Piles, days 3-4 p = %d, t = %d \n',pPilesClusDays,s.tstat)


%% which is closer task

subplot(2,3,3)
hold on
title('Which is closer?')
% Plot individual data points
scatter(ones(nSub,1),100*D.closer.hex_day1_ratioCorrect,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(2*ones(nSub,1),100*D.closer.hex_day2_ratioCorrect,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(3*ones(nSub,1),100*D.closer.clus_day3_ratioCorrect,3,'filled','b','MarkerFaceAlpha',0.5)
scatter(4*ones(nSub,1),100*D.closer.clus_day4_ratioCorrect,3,'filled','b','MarkerFaceAlpha',0.5)

% plot means and SEM error bars
wBar = 0.15;
plot([1+0.2-wBar, 1+0.2+wBar], 100*[mean(D.closer.hex_day1_ratioCorrect),mean(D.closer.hex_day1_ratioCorrect)],'k','lineWidth',2)
plot([2-0.2-wBar, 2-0.2+wBar], 100*[mean(D.closer.hex_day2_ratioCorrect),mean(D.closer.hex_day2_ratioCorrect)],'k','lineWidth',2)
plot([3+0.2-wBar, 3+0.2+wBar], 100*[mean(D.closer.clus_day3_ratioCorrect),mean(D.closer.clus_day3_ratioCorrect)],'k','lineWidth',2)
plot([4-0.2-wBar, 4-0.2+wBar], 100*[mean(D.closer.clus_day4_ratioCorrect),mean(D.closer.clus_day4_ratioCorrect)],'k','lineWidth',2)

errorbar(1+0.2,100*mean(D.closer.hex_day1_ratioCorrect),100*std(D.closer.hex_day1_ratioCorrect)/sqrt(nSub),'k')
errorbar(2-0.2,100*mean(D.closer.hex_day2_ratioCorrect),100*std(D.closer.hex_day2_ratioCorrect)/sqrt(nSub),'k')
errorbar(3+0.2,100*mean(D.closer.clus_day3_ratioCorrect),100*std(D.closer.clus_day3_ratioCorrect)/sqrt(nSub),'k')
errorbar(4-0.2,100*mean(D.closer.clus_day4_ratioCorrect),100*std(D.closer.clus_day4_ratioCorrect)/sqrt(nSub),'k')

% Plot lines connecting the data points
for iSub=1:nSub
    plot([1,2],100*[D.closer.hex_day1_ratioCorrect(iSub),D.closer.hex_day2_ratioCorrect(iSub)],'color',[1 0 0 0.2])
    plot([3,4],100*[D.closer.clus_day3_ratioCorrect(iSub),D.closer.clus_day4_ratioCorrect(iSub)],'color',[0 0 1 0.2])
end

xlim([0,5]);
ylabel('% correct')
labels = {'hex day1','hex day2','comm day3', 'comm day4'};
set(gca,'FontSize',14,'XTick',1:4,'XTickLabels',labels,'XTickLabelRotation', 45)



% stats for whether answers oneach day were different from chanve
[~, pCloserHexDay1,~,s] = ttest(D.closer.hex_day1_ratioCorrect,0.5,'Tail','right');
fprintf('Closer, day 1 p = %d, t = %d \n',pCloserHexDay1,s.tstat)
[~, pCloserHexDay2,~,s] = ttest(D.closer.hex_day2_ratioCorrect,0.5,'Tail','right');
fprintf('Closer, day 2 p = %d, t = %d \n',pCloserHexDay2,s.tstat)
[~, pCloserClusDay3,~,s] = ttest(D.closer.clus_day3_ratioCorrect,0.5,'Tail','right');
fprintf('Closer, day 3 p = %d, t = %d \n',pCloserClusDay3,s.tstat)
[~, pCloserClusDay4,~,s] = ttest(D.closer.clus_day4_ratioCorrect,0.5,'Tail','right');
fprintf('Closer, day 4 p = %d, t = %d \n',pCloserClusDay4,s.tstat)

% stats for improvement beween days
[~, pCloserHexDays,~,s] = ttest2(D.closer.hex_day2_ratioCorrect,D.closer.hex_day1_ratioCorrect,'Tail','right');
fprintf('Closer, days 1-2 p = %d, t = %d \n',pCloserHexDays,s.tstat)
[~, pCloserClusDays,~,s] = ttest2(D.closer.clus_day4_ratioCorrect,D.closer.clus_day3_ratioCorrect,'Tail','right');
fprintf('Closer, days 3-4 p = %d, t = %d \n',pCloserClusDays,s.tstat)

%% navigation task


for iDist = 1:3
    subplot(2,3,iDist+3)
    hold on
    
    title(['Navig dist ' num2str(iDist + 1)]) % between 2 and 4
    
    % individual datapoints for current initial distance
    day1 = mean(D.navig.hex_day1_nSteps(:,:,iDist),2);
    day2 = mean(D.navig.hex_day2_nSteps(:,:,iDist),2);
    day3 = mean(D.navig.clus_day3_nSteps(:,:,iDist),2);
    day4 = mean(D.navig.clus_day4_nSteps(:,:,iDist),2);
    
    % Plot individual data points
    scatter(ones(nSub,1),day1,3,'filled','r','MarkerFaceAlpha',0.5)
    scatter(2*ones(nSub,1),day2,3,'filled','r','MarkerFaceAlpha',0.5)
    scatter(3*ones(nSub,1),day3,3,'filled','b','MarkerFaceAlpha',0.5)
    scatter(4*ones(nSub,1),day4,3,'filled','b','MarkerFaceAlpha',0.5)
    
    % plot means and SEM error bars
    wBar = 0.15;
    plot([1+0.2-wBar, 1+0.2+wBar], [mean(day1),mean(day1)],'k','lineWidth',2)
    plot([2-0.2-wBar, 2-0.2+wBar], [mean(day2),mean(day2)],'k','lineWidth',2)
    plot([3+0.2-wBar, 3+0.2+wBar], [mean(day3),mean(day3)],'k','lineWidth',2)
    plot([4-0.2-wBar, 4-0.2+wBar], [mean(day4),mean(day4)],'k','lineWidth',2)
    
    errorbar(1+0.2,mean(day1),std(day1)/sqrt(nSub),'k')
    errorbar(2-0.2,mean(day2),std(day2)/sqrt(nSub),'k')
    errorbar(3+0.2,mean(day3),std(day3)/sqrt(nSub),'k')
    errorbar(4-0.2,mean(day4),std(day4)/sqrt(nSub),'k')
    
    % Plot lines connecting the data points
    for iSub=1:nSub
        plot([1,2],[day1(iSub),day2(iSub)],'color',[1 0 0 0.2])
        plot([3,4],[day3(iSub),day4(iSub)],'color',[0 0 1 0.2])
    end
    
    xlim([0,5]);
    ylabel('nSteps')
    labels = {'hex day1','hex day2','comm day3', 'comm day4'};
    set(gca,'FontSize',14,'XTick',1:4,'XTickLabels',labels,'XTickLabelRotation', 45)
    
    % stats for improvement beween days
    [~, pNavigHexDays,~,s] = ttest2(day1,day2,'Tail','right');
    fprintf('initDist %d, diff days Hex p = %d, t = %d \n',iDist+1,pNavigHexDays,s.tstat)
    [~, pNavigClusDays,~,s] = ttest2(day3,day4,'Tail','right');
    fprintf('initDist %d, diff days Comm p = %d, t = %d \n',iDist+1,pNavigClusDays,s.tstat)
end






