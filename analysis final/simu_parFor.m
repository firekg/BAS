%%
parpool(10);

%% add pixel to SIM (data saved analysisSIMwPix-2015-06-16)
for simu = 1:10
for parti = 1:3
    for cond = 2:5 %1 is empty
        SIMwPix(parti,cond,simu) = addPixToData(SIM_batch(parti,cond,simu));
    end
end
fprintf('simu %d done.\n', simu);
end

% fill SIM(:,1) with ffAL
load('SCdata.mat', 'SCAL');
load('EPdata.mat', 'EPAL');
load('BZdata.mat', 'BZAL');
ffAL = feedForwardRevPos(SCAL);
SIMwPix(1,1,1) = addPixToData(ffAL);
ffAL = feedForwardRevPos(EPAL);
SIMwPix(2,1,1) = addPixToData(ffAL);
ffAL = feedForwardRevPos(BZAL);
SIMwPix(3,1,1) = addPixToData(ffAL);
for simu = 2:10
    SIMwPix(:,1,simu) = SIMwPix(:,1,1);
end

%% simulation prop for the above
PROP = [];
ntask = 10;
taskChoice = 5; 
parfor task = 1:ntask
    rng('shuffle');
    rng;
    PROP(task).prop = prop_sim(SIMwPix(:,:,taskChoice)); %with SIM
    %PROP(task).prop = prop_sim([]); %without SIM
end

%% add pixel to DActive (data saved analysisDActivewPix-2015-06-25)
for session = 1:11
    for parti = 1:3
        DALwPix(parti,session) = addPixToData(DActive(parti,session));
    end
end

%% information efficiency of above
INFOA = [];
parS = [];
ntask = 10;
parfor task = 1:ntask
    rng('shuffle');
    rng;    
    [INFOA(task).info, parS(task).par] = info_progress_v2(DALwPix);
end

%% 2015-07-02 ideal BAS with (cond=1) and without (cond=2) Prior
for parti = 1:3
    for cond = 1:2
        SIM_idealBAS_PBwPix(parti,cond) = addPixToData(SIM_idealBAS_PB(parti,cond));
    end
end

%% continue
PROP = [];
ntask = 10;
parfor task = 1:ntask
    rng('shuffle');
    rng;
    PROP(task).prop = prop_sim(SIM_idealBAS_PBwPix);
end

%% put old prop_sim into this new format
% nsimu = 10;
% for isimu = 1:nsimu;
%     load(['expPROP_data/expPROPRefitNoise',num2str(isimu),'.mat']);
%     PROPre(isimu).prop = PROP;
% end

%% 2015-09-04 to 24
% simulating MaxEnt in SIM_MaxEnt(:,2,1)
%[SIM_MaxEnt] = simu_BASs();
% % saved as SIM_MaxEnt-2015-09-04.mat

[REV, ~, ~] = revmap_sim(SIM_MaxEnt);
ntask = 10;
parfor task = 1:ntask
    [SMboot] = SM_bootstrap(REV);
    CombSM(task).SMboot = SMboot;
end
% saved as CombSM-MaxEnt-2015-09-07.mat

[CORR] = corr_comb_SMboot([], CombSM);
plot_corr(CORR);
% export_fig('corrMaxEnt.pdf');

%% 2015-09-07 to 24: dealing with order shuffling
clear;
load('analysisSIMwPix-2015-06-16.mat')
[REV, ~, ~] = revmap_sim(SIMwPix);
ntask = 10;
parfor task = 1:ntask
    % 2015-09-24 check: I think I have gotten the wrong weights 
    %[SMboot] = SM_bootstrap(REV);
    [SMboot] = SM_bootstrap_shuffle2(REV);
    CombSM(task).SMboot = SMboot;
end
% saved CombSM as CombSM-Shuffle-2015-09-07.mat
[CORR] = corr_comb_SMboot([], CombSM);
plot_corr(CORR);

%% 2015-09-10: dealing with correct vs incorrect density maps
load('analysisSIMwPix-2015-06-16.mat')
% % correct-incorrect trials - changes inside revmap_sim
[REV, RrevMap, RdrevMap] = revmap_sim(SIMwPix);
plot_rev_maps(RrevMap, RdrevMap);

% ntask = 10;
% parfor task = 1:ntask
%     [SMboot] = SM_bootstrap(REV);
%     CombSM(task).SMboot = SMboot;
% end
% % saved CombSM as CombSM-correct-2015-09-10.mat
% % saved CombSM as CombSM-incorrect-2015-09-10.mat
% [CORR] = corr_comb_SMboot(RdrevMap, CombSM);
% plot_corr(CORR);

%% 2015-09-11 BAS and MaxEnt percentile
% B = BASmap();
% Bcomb = [];
% ntask = 10;
% parfor task = 1:ntask
%     rng('shuffle');
%     rng;
%     Bcomb(task).B = BASptile();
% end
% % saved Bcomb as BASptile-2015-09-11.mat

%% 2015-09-22: Shuffling revealing order in info curves analysis
% PROP = [];
% ntask = 10;
% taskChoice = 5; 
% parfor task = 1:ntask
%     rng('shuffle');
%     rng;
%     PROP(task).prop = prop_sim();
% end

% [PERF, Pmix, Bias, BiasMix, mPERF, mPmix, lpxPERF, lpxPmix, INFO, INFOC] = perf_participants(PROP);
% plot_perf_info(PERF, INFO, mPERF);

%% 2015-09-24: redo density maps and correlation because of wrong weights

% load('analysisSIMwPix-2015-06-16.mat')
% [REV, RrevMap, RdrevMap] = revmap_sim(SIMwPix);
% plot_rev_maps(RrevMap, RdrevMap);

% load('analysisSIMwPix-2015-06-16.mat')
% [REV, ~, ~] = revmap_sim(SIMwPix(:,:,5));
% ntask = 10;
% parfor task = 1:ntask
%     [SMboot] = SM_bootstrap(REV);
%     CombSM(task).SMboot = SMboot;
% end
% [CORR] = corr_comb_SMboot([], CombSM);
% plot_corr(CORR);

