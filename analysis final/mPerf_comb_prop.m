function [mPERF, mPERFC] = mPerf_comb_prop()
% Depreciated -- because it's a SUBSET of perf_participants

% PERF has 5 layers: 
% revealing number x25 (1:25); 
% parti x4; 
% condiitons x 3 (AL,PLR,BAS);
% moments x2 (mean,SD)
% number of prop

load('SCdata.mat', 'SCAL', 'SCPLRbt', 'SCPLBbt', 'SCPLaBbt'); 
load('EPdata.mat', 'EPAL', 'EPPLRbt', 'EPPLBbt', 'SCPLaBbt'); 
load('BZdata.mat', 'BZAL', 'BZPLRbt', 'BZPLBbt', 'SCPLaBbt'); 
D(1,1).ansReal = SCAL.AnswerReal;
D(2,1).ansReal = EPAL.AnswerReal;
D(3,1).ansReal = BZAL.AnswerReal;
D(1,2).ansReal = SCPLRbt.AnswerReal;
D(2,2).ansReal = EPPLRbt.AnswerReal;
D(3,2).ansReal = BZPLRbt.AnswerReal;
D(1,3).ansReal = SCPLBbt.AnswerReal;
D(2,3).ansReal = EPPLBbt.AnswerReal;
D(3,3).ansReal = BZPLBbt.AnswerReal;

nprop = 20;
nsimu = 10;
ncond = 3;

mPERFC = nan(25, 4, ncond, 2, nprop*nsimu);
for isimu = 1:nsimu
    %load(['PROP_data/PROP',num2str(isimu),'.mat']);
    %load(['newNoisePROP_data/newPNoisePROP',num2str(isimu),'.mat']);
    load(['expPROP_data/expPROP',num2str(isimu),'.mat']);
    for iprop = 1:nprop
        mPERF = mPerf_prop(PROP(:,:,iprop), D);        
        mPERFC(:,:,:,:, (isimu-1)*nprop+iprop) = mPERF;
    end
end

mPERF = nanmean(mPERFC,5);

end

%% Compute model performance
function [mPERF] = mPerf_prop(PROP, D)
% PROP has 2 layers: participants x3 (SC,EP,BZ); conditions x4 (ALprop, PLRprop, BASprop, mBASprop)
% mPERF has 4 layers: revealing number x25 (1:25); parti x4; condiitons x3 (AL,PLR,BAS); stat x2 (mean, sd)

ncond = 3;
mPERF = nan(25,4,ncond,2);
for parti = 1:3
    for cond = 1:ncond
        m = PROP(parti, cond).mla; % marginal likelihood for patchy        
        ans = D(parti,cond).ansReal;
        correct = calCorret(m, ans);
        mPERF(:, parti, cond, 1) = nanmean(correct);  
        mPERF(:, parti, cond, 2) = sqrt(nanvar(correct)./sum(~isnan(correct)));        
    end
   %disp(sprintf('parti=%d;',parti));
end
for parti=4
    for cond = 1:ncond
        m=[PROP(1, cond).mla; PROP(2, cond).mla; PROP(3, cond).mla]; 
        ans = [D(1, cond).ansReal; D(2, cond).ansReal; D(3, cond).ansReal];
        correct = calCorret(m, ans);
        mPERF(:, parti, cond, 1) = nanmean(correct);  
        mPERF(:, parti, cond, 2) = sqrt(nanvar(correct)./sum(~isnan(correct)));
    end
    %disp(sprintf('parti=%d;',parti));
end

end
%%
function correct = calCorret(m, ans)
    nrev = 25;
    [nrow, ncol] = size(m);
    choice = nan(nrow, ncol);
    choice(m>0.5) = 1;
    choice(m<0.5) = 2;
    choice(:, 1) = round(rand(nrow,1)) + 1;
    ansM = repmat(ans, 1, nrev);
    correct = (choice == ansM) + 0; %convert to double
    correct(isnan(m)) = nan;
end


