close all
clear 
clc

origSubName = {'Martina_Kavanova','Yuechen_Ren','Antilia_Virginie','Ujjayanta_Bhaumik','Delia_ciobotaru','Wu_Trista','Zixi_Liu','Josh_Spowage','Han Jien_Chang','Anusha_Leonard','Rachel_Loh','YIFEI_LANG','Neza_Vehar','zalike_keskin','Bernice_Tan','Jallene_Chua','Aslam_Norhan','Jingyi_Ye','Sofya_Vakhonina','Lucia_Canto','Saul_Westbrook','Nidhi_Singh','Eoin_Bentick','zahid_matzin','Zhen Jie_Low','Ruslan_Seletskiy','Shirley_Feng','Candice_Hughes'};
origSubNum = [51,50,49,48,46,45,44,43,42,40,38:-1:34,32:-1:28,26:-1:24,22:-1:18];


% first swap the order of the subject names and indeces
nSub = numel(origSubName);%nsb
subName = cell(nSub, 1); % vsubName
for i = 1:nSub
    subName{i} = origSubName{nSub - i + 1};
end
subNum = origSubNum(end:-1:1);

subNameNew = cell(1,28);
for iSub = 1:length(subNameNew)
    subNameNew{iSub}= ['sub-' num2str(iSub,'%02.f')];    
end

% read tables
origDir = '/vols/Scratch/abaram/shirley/codes/beh';
middleTable = readtable(fullfile(origDir,'isMiddleTable260821.csv')); %mM
pilesTable  = readtable(fullfile(origDir,'isPileTable260821.csv')); %mP
closerTable = readtable(fullfile(origDir,'qTable101220.csv')); %mQ
navigTable = readtable(fullfile(origDir,'TaskTable260821.csv')); %mT

middleTable.Var1 = lower(middleTable.Var1); %namesM
pilesTable.Var1 = lower(pilesTable.Var1); %namesP
closerTable.Var1 = lower(closerTable.Var1); %namesQ
navigTable.Var1 = lower(navigTable.Var1); %namesT

newDir = '/vols/Scratch/abaram/shirley/alon/beh/training';

indsSubMiddleAll=[];
indsSubPilesAll=[];
indsSubCloserAll=[];
indsSubNavigAll=[];

for iSub = 1:length(subName)
    indsSubMiddle = find(strcmpi(middleTable.Var1,subName{iSub}));
    indsSubMiddleAll = vertcat(indsSubMiddleAll,indsSubMiddle);
    middleTable.Var1(indsSubMiddle) = subNameNew(iSub);
    
    indsSubPiles = find(strcmpi(pilesTable.Var1,subName{iSub}));
    indsSubPilesAll = vertcat(indsSubPilesAll,indsSubPiles);    
    pilesTable.Var1(indsSubPiles) = subNameNew(iSub);
    
    indsSubCloser = find(strcmpi(closerTable.Var1,subName{iSub}));
    indsSubCloserAll = vertcat(indsSubCloserAll,indsSubCloser);
    closerTable.Var1(indsSubCloser) = subNameNew(iSub);    
        
    indsSubNavig = find(strcmpi(navigTable.Var1,subName{iSub}));
    indsSubNavigAll = vertcat(indsSubNavigAll,indsSubNavig);    
    navigTable.Var1(indsSubNavig) = subNameNew(iSub);
end
middleTableNew = middleTable(indsSubMiddleAll,:);
pilesTableNew = pilesTable(indsSubPilesAll,:);
closerTableNew = closerTable(indsSubCloserAll,:);
navigTableNew = navigTable(indsSubNavigAll,:);

writetable(middleTableNew, fullfile(newDir,'middleTable.csv'));
writetable(pilesTableNew, fullfile(newDir,'pilesTable.csv'));
writetable(closerTableNew, fullfile(newDir,'closerTable.csv'));
writetable(navigTableNew, fullfile(newDir,'navigTable.csv'));
