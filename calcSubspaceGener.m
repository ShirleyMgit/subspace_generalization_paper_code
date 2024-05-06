function H = calcSubspaceGener(rawData,SPM)
% rawData: raw data from searchlight
% calculate the subspace generalisation (area under the curve of explained
% variance as a function of number of EVs)

flagCal = 1; % flag whether to caclulate subspace generalisatoin, will change to 0 only if there is a problem with NaNs or Infs. 

% if there are any voxels in the SL that have data with NaN of Inf, get rid
% of them 
isNanOrInf = isnan(mean(rawData,1)) + isinf(mean(rawData,1));
rawData = rawData(:,isNanOrInf==0);

% pre - whitening and calculate betas:
% default: normalize using all runs (overall):
[betaAll00,~,~,~,~,~]=rsa.spm.noiseNormalizeBeta(rawData,SPM);

% which betas: (for later analysis)
nPiles= 10;
nMaps = 4;

pilesBetasInds = ones(1,nPiles*nMaps);
zeroBetas = 10 + 5 +2 +7; % 10: piles in map 5; 5: catch trials regs, one per block; 2: catch trial question and start experiment reg; 7: nuisance regressor (6 motion + CSF)
betasToUse_singleRun =  [zeros(1,5) pilesBetasInds zeros(1,zeroBetas)]; % the first 5 0s: The regs modelling the random walks at the start of each map. 
betasToUse = [betasToUse_singleRun betasToUse_singleRun betasToUse_singleRun betasToUse_singleRun];
% get indeces of betas to use
b_idx0 = 1:length(betasToUse);
b_idx = b_idx0(betasToUse==1);

% check there are voxels with no betas that are NaNs or Infs
betasToUse = betaAll00(b_idx,:);
isNanOrInf = isnan(mean(betasToUse,1))  + isinf(mean(betasToUse,1));

if sum(isNanOrInf>0)
    disp(' \n \n \n!!!!! There was a NaN or Inf, did not calculate !!!!! \n \n \n ')
    warning(' \n \n \n!!!!! There was a NaN or Inf, did not calculate !!!!! \n \n \n ')
    flagCal = 0;
end


if flagCal==1 % if betas are OK calculate, otherwise put -10000 at that voxel
    nVoxels = length(betasToUse(1,:));
    betasToUse = reshape(betasToUse,[40,4,nVoxels]); % nPiles*nMaps, nRuns, nVoxels
    betaMeanLOO = zeros(40,4,nVoxels);
    
    nr = 1:4;
    
    %% LOO averages
    for r = 1:4
        betaMeanLOO(:,r,:) = mean(betasToUse(:,nr~=r,:),2); % mean of the betas in 3 of the 4 runs
    end
    
    s1L = zeros(4,1);
    s2L = zeros(4,1);
    s3L = zeros(4,1);
    s4L = zeros(4,1);
    
    for r = 1:4
        
        % Do SVD on the data averaged from the 3 runs except current run. 
        % U is the spatial components, S the scalings (eigenvalues), V the
        % temporal components. 
        
        %hex 1
        hex1avLoo = squeeze(betaMeanLOO(1:10,r,:));
        hex1avLoo = hex1avLoo - repmat(mean(hex1avLoo),nPiles,1); % substract the mean over condition
        [U1,~,~] = svd(hex1avLoo'*hex1avLoo,'econ');%
        
        % hex 2
        hex2avLoo = squeeze(betaMeanLOO(11:20,r,:));
        hex2avLoo = hex2avLoo - repmat(mean(hex2avLoo),nPiles,1);
        [U2,~,~] = svd(hex2avLoo'*hex2avLoo,'econ');%
        
        % cluster 1:
        clus1avLoo = squeeze(betaMeanLOO(21:30,r,:));
        clus1avLoo = clus1avLoo - repmat(mean(clus1avLoo),nPiles,1);
        [U3,~,~] = svd(clus1avLoo'*clus1avLoo,'econ');%
        
        % cluster 2:
        clus2avLoo = squeeze(betaMeanLOO(31:40,r,:));
        clus2avLoo = clus2avLoo - repmat(mean(clus2avLoo),nPiles,1);
        [U4,~,~] = svd(clus2avLoo'*clus2avLoo,'econ');%
        
        % Calculate SVD on the data from the run itself. 
        % Sr are the scalings (eigenvalues)
        
        %hex 1
        hex1run = squeeze(betasToUse(1:10,r,:));
        hex1run =  hex1run - repmat(mean(hex1run,1),nPiles,1);
        [~,S1r,~] = svd(hex1run'*hex1run,'econ');%

        %hex 2
        hex2run = squeeze(betasToUse(11:20,r,:));
        hex2run = hex2run - repmat(mean(hex2run,1),nPiles,1);
        [~,S2r,~] = svd(hex2run'*hex2run,'econ');

        %cluster 1:
        clus1run = squeeze(betasToUse(21:30,r,:));
        clus1run = clus1run - repmat(mean(clus1run),nPiles,1);
        [~,S3r,~] = svd(clus1run'*clus1run,'econ');

        %cluster 2:
        clus2run = squeeze(betasToUse(31:40,r,:));
        clus2run = clus2run - repmat(mean(clus2run,1),nPiles,1);
        [~,S4r,~] = svd(clus2run'*clus2run,'econ');
        
        %projections:
        % To hex 1:
        pr1 = zeros(4,nVoxels);
        % project to eigenvectors
        P11 = (hex1run*U1).^2;
        P21 = (hex2run*U1).^2;
        P31 = (clus1run*U1).^2;
        P41 = (clus2run*U1).^2;
        % sum variance explained
        pr1(1,:) = sum(P11)/trace(S1r);
        pr1(2,:) = sum(P21)/trace(S2r);
        pr1(3,:) = sum(P31)/trace(S3r);
        pr1(4,:) = sum(P41)/trace(S4r);
        % calculate area under the curve
        v1L = cumsum(pr1,2);
        s1L(:,r)  = sum(v1L,2);
        
        % To hex 2:
         pr2 = zeros(4,nVoxels);
        % project to eigenvectors
        P12 = (hex1run*U2).^2;
        P22 = (hex2run*U2).^2;
        P32 = (clus1run*U2).^2;
        P42 = (clus2run*U2).^2;
      % sum variance explained
        pr2(1,:) = sum(P12)/trace(S1r);
        pr2(2,:) = sum(P22)/trace(S2r);
        pr2(3,:) = sum(P32)/trace(S3r);
        pr2(4,:) = sum(P42)/trace(S4r);
        % calculate area under the curve
        v2L = cumsum(pr2,2);
        s2L(:,r)  = sum(v2L,2);
        
        % To cluster 1:
        pr3 = zeros(4,nVoxels);
        % project to eigenvectors
        P13 = (hex1run*U3).^2;
        P23 = (hex2run*U3).^2;
        P33 = (clus1run*U3).^2;
        P43 = (clus2run*U3).^2;
      % sum variance explained
        pr3(1,:) = sum(P13)/trace(S1r);
        pr3(2,:) = sum(P23)/trace(S2r);
        pr3(3,:) = sum(P33)/trace(S3r);
        pr3(4,:) = sum(P43)/trace(S4r);
        % calculate area under the curve
        v3L = cumsum(pr3,2);
        s3L(:,r)  = sum(v3L,2);
        
        % To cluster 2:
        pr4 = zeros(4,nVoxels);
        % project to eigenvectors
        P14 = (hex1run*U4).^2;
        P24 = (hex2run*U4).^2;
        P34 = (clus1run*U4).^2;
        P44 = (clus2run*U4).^2;
      % sum variance explained, and normalise by the total variance in the
      % projected data 
        pr4(1,:) = sum(P14)/trace(S1r);
        pr4(2,:) = sum(P24)/trace(S2r);
        pr4(3,:) = sum(P34)/trace(S3r);
        pr4(4,:) = sum(P44)/trace(S4r);
        % calculate area under the curve
        v4L = cumsum(pr4,2);
        s4L(:,r)  = sum(v4L,2);
        
    end
    
    % get the average of the area under the curve across all runs
    ms1 = mean(s1L,2); 
    ms2 = mean(s2L,2);
    ms3 = mean(s3L,2);
    ms4 = mean(s4L,2);
    
    H = zeros(1,16);
    H(1:4) = ms1; % projections on eigenvectors from Hex1 (actually from the mean betas of 3 runs of Hex 1)
    H(5:8) = ms2; % projections on eigenvectors from Hex2
    H(9:12) = ms3; % projections on eigenvectors from cluster1 
    H(13:16) = ms4; % projections on eigenvectors from cluster2
    
else
    H = -10000*ones(1,16);
end
