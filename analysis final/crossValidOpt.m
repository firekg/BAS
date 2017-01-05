%% 2015-09-08 Cross-validation of models-response to eLife review
% based on saved data:
% probParalShifts-Lpr-Sigo0d1-1d0.mat (shiftVec in units of 770/20 pixels...)
% probScaleDown-Lpr-Sigo0d1-0d5.mat & probScaleDown-Lpr-Sigo0d6-1d0.mat
% code mostly copied from optimize_para

%% split trials into 10 parts
clear;
% for one model, one parti (~10 min)
% see if enough memory to use parfor; yes!
load('probParalShifts-Lpr-Sigo0d1-1d0.mat', 'LprEP');
% L1 = load('probScaleDown-Lpr-Sigo0d1-0d5.mat', 'LprBZ');
% L2 = load('probScaleDown-Lpr-Sigo0d6-1d0.mat', 'LprBZ');
% LprBZ = cat(2, L1.LprBZ, L2.LprBZ);
% These are the lpr for the choice (a_Get_Model_Perf)
% so after softmaxing, these are the prob for choices (optfunc_lpr)
% so 1 - this prob = prob of making the other choice, ie. prediction error
% so should be right

%%
% Lpr 4D: trials, sigo, shift, samples
% Lpr = LprEP(:,:,1,:); %for model 1-2
Lpr = LprEP; %for model 3-6
[ntrials, nsigo, nshift, nsamp] = size(Lpr)
 
% setup permutation
ntr_sub = floor(ntrials/10);
perm = randperm(ntrials);
permLpr = Lpr(perm,:,:,:); %is this legit?

% setup optimization
seedn = 1; % CHECK: should match that inside opt_lpr
paran = 2; % lbd, lapse
opt_para = nan(seedn,paran,nsigo,nshift);
opt_val = nan(seedn,nsigo,nshift);
x0_para = nan(seedn,paran,nsigo,nshift);

crossErr = zeros(10,1);
for crossV = 1:10

    ind1 = (crossV-1)*ntr_sub + 1;
    ind2 = (crossV)*ntr_sub;
    % test data
    testLpr = permLpr(ind1:ind2,:,:,:);    
    % training data
    trainLpr = permLpr;
    trainLpr(ind1:ind2,:,:,:) = [];
    
    % optimize lbd and lapse
    % CHECK a_Model has the right parameter shift/scale method
    % CHECK opt_lpr has the right bounds on lbd and lapse
    for isigo = 1:nsigo
        parfor ishift = 1:nshift
            [opt_temp, x0_para(:,:,isigo,ishift)] = opt_lpr( trainLpr(:,isigo,ishift,:) );
            opt_para(:,:,isigo,ishift) = opt_temp(:,1:2);
            opt_val(:,isigo,ishift) = opt_temp(:,3);
        end
        fprintf('isigo %d\n', isigo); 
    end
    
    % find optimal val and parameter
    val = mean(opt_val, 1); % unncessary with seedn=1
    val = reshape(val, nsigo, nshift);
    [minval, minInd] = min(val(:));
    [indrow, indcol] = ind2sub([nsigo, nshift], minInd);
    parabd = reshape(opt_para(1,1,:,:),nsigo, nshift);
    paralapse = reshape(opt_para(1,2,:,:),nsigo, nshift);
    % fprintf('opt lbd: %f; \n', parabd(indrow, indcol));
    % fprintf('opt lapse: %f; \n', paralapse(indrow, indcol));
    % fprintf('min val = %f; \n', minval);

    % prediction error
    optBd = exp(parabd(indrow, indcol));
    optLapse = paralapse(indrow, indcol);
    softlik = 1./(1+exp(-optBd*testLpr(:, indrow, indcol, :)));
    lapse_soft_lik = (1-optLapse)*softlik + 0.5*optLapse;
    trialAcc = mean(lapse_soft_lik, 4); %marginalize over perception noise
    predAcc = geomean(trialAcc); 
    crossErr(crossV) = 1-predAcc;
    
    fprintf('cross val num %d\n', crossV); 

end
fprintf('average CV err = %f\n', mean(crossErr)); 

%% record data by hand
% model number as indexed in the paper

% Arithmetic mean; geometric mean
% parti SC: trial size = 1502
% parti SC, model 1 = 0.4265; 0.435209
% parti SC, model 2 = 0.4115; 0.435718
% parti SC, model 3 = 0.3812; 0.426977
% parti SC, model 4 = 0.3820; 0.429026
% parti SC, model 5 = 0.3946; 0.437074
% parti SC, model 6 = 0.3931; 0.436872
% SC = [0.4265, 0.4115, 0.3812, 0.3820, 0.3946, 0.3931];
SC = [0.4352, 0.4357, 0.4269, 0.4290, 0.4370, 0.4368];

% Arithmetic mean; geometric mean
% parti EP: trial size = 1530
% parti EP, model 1 = 0.4482; 0.439324
% parti EP, model 2 = 0.4206; 0.438568
% parti EP, model 3 = 0.3859; 0.433775
% parti EP, model 4 = 0.3853; 0.432116
% parti EP, model 5 = 0.3955; 0.439230
% parti EP, model 6 = 0.3963; 0.439986
% EP = [0.4482, 0.4206, 0.3859, 0.3853, 0.3955, 0.3963];
EP = [0.4393, 0.4385, 0.4337, 0.4321, 0.4392, 0.4399];

% Arithmetic mean; geometric mean
% parti BZ: trial size = 1376
% parti BZ, model 1 = 0.4102; 0.401808
% parti BZ, model 2 = 0.3512; 0.397120
% parti BZ, model 3 = 0.3308; 0.395252
% parti BZ, model 4 = 0.3264; 0.389997
% parti BZ, model 5 = 0.3249; 0.392066
% parti BZ, model 6 = 0.3250; 0.388483
% BZ = [0.4102, 0.3512, 0.3308, 0.3264, 0.3249, 0.3250]; 
BZ = [0.4018, 0.3971, 0.3952, 0.3899, 0.3920, 0.3884]; 

trs = [1502, 1530, 1376];
w = trs/sum(trs);
CV = w(1)*SC + w(2)*EP + w(3)*BZ;

