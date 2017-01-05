%% gradient-based local optimization -- trial-by-trial performance
D=a_Combine_Data(SCPLRbt,SCPLBbt,SCPLaBbt);
% DIMSCAL=a_DIM(SCAL);  
% DIMSCPLR=a_DIM(SCPLRbt);  
% DIMSCPLB=a_DIM(SCPLBbt);  
% DIMSCPLaB=a_DIM(SCPLaBbt);
DIM=[DIMSCPLR,DIMSCPLB,DIMSCPLaB];
n=1; w=1;

% initialize AUX
% Moments should be defined:
% defined in Test_revealing_vs_cateogry_stat
AUX.Moments=Moments;
AUX.pvp=[]; AUX.pvd=[];

opt_data_all=zeros(n,10);
x0_data_all=zeros(n,9);
for i=1:n
    %lb = [lsigo,    phit,     lphis,  phis0,        lbx,        lbd, lapse, prpa];
    lb= [log(0.1),   0.5,  log(0.1),    0.5,  log(1e-3),  log(1e-2),  0.01,  0.5];
    ub= [log(1.0),     5,   log(20),     10,     log(1),    log(10),   0.5,  0.5];
    %x0 = [log(0.2),   100,    log(1),    100,  log(1e-5),  log(1),    0,   0.5];
    x0 = rand(1,8).*(ub-lb)+lb;
    f = @(x)a_AllTrial_SoftLikLapse(x,w,D,DIM,AUX);
    options = optimset('Algorithm','interior-point');
    options = optimset(options,'Display','iter');
    options = optimset(options,'GradObj','on');
    options = optimset(options,'TolX',1e-6,'TolFun',1e-6);
    [x,fval,exitflag,output] = fmincon(f,x0,[],[],[],[],lb,ub,[],options);
    opt_data_all(i,:)=[x,w,fval];
    x0_data_all(i,:)=[x0,w];
    disp(sprintf('i=%d;',i));
end
%% only want to change lsigo and fix all other parameters
% after making sure the above code runs, can then do this

AUX.Moments=Moments;
AUX.pvp=[]; AUX.pvd=[];

% % load all 3 of SCdata, EPdata, and BZdata, then load the following:
% AL = a_Combine_Data(SCAL, EPAL, BZAL);
% PLRbt = a_Combine_Data(SCPLRbt, EPPLRbt, BZPLRbt);
% PLBbt = a_Combine_Data(SCPLBbt, EPPLBbt, BZPLBbt);
% PLaBbt = a_Combine_Data(SCPLaBbt, EPPLaBbt, BZPLaBbt);
% DIMAL=a_DIM(AL);  
% DIMPLR=a_DIM(PLRbt);  
% DIMPLB=a_DIM(PLBbt);  
% DIMPLaB=a_DIM(PLaBbt);
% 
% disp(sprintf('Done loading.')); 

opt = [log(0.1),   100,    log(1),    100,  log(1e-5),  log(100),    0,   0.5, 1];

testmode = 2;  %1 is least-square; 2 is likelihood
if(testmode==1)
    [PerfAL, sdAL, ~, x] = a_Get_Human_Perf(AL);
    [PerfPLR, sdPLR, ~, x] = a_Get_Human_Perf(PLRbt);
    [PerfPLB, sdPLB, ~, x] = a_Get_Human_Perf(PLBbt);
    [PerfPLaB, sdPLaB, ~, x] = a_Get_Human_Perf(PLaBbt);
    test = zeros(10,10);
end
nsample = 50;
nincr = 10;
ProbAL = nan(AL.Trials, nincr, nsample);
ProbPLR = nan(PLRbt.Trials, nincr, nsample);
ProbPLB = nan(PLBbt.Trials, nincr, nsample);
ProbPLaB = nan(PLaBbt.Trials, nincr, nsample);
j=0;
for i = linspace(0.6, 1.5, nincr)
    j=j+1;
    for isample = 1:nsample
        opt(1)=log(i);
        [mPerfAL, msdAL, ProbAL(:,j,isample)] = a_Get_Model_Perf(AL,DIMAL,AUX,opt);
        [mPerfPLR, msdPLR, ProbPLR(:,j,isample)] = a_Get_Model_Perf(PLRbt,DIMPLR,AUX,opt);
        [mPerfPLB, msdPLB, ProbPLB(:,j,isample)] = a_Get_Model_Perf(PLBbt,DIMPLB,AUX,opt);
        [mPerfPLaB, msdPLaB, ProbPLaB(:,j,isample)] = a_Get_Model_Perf(PLaBbt,DIMPLaB,AUX,opt);
        if(testmode==1)
            chi2 = sum( (PerfAL - mPerfAL).^2 ) +...
                    sum( (PerfPLR - mPerfPLR).^2 ) +...
                    sum( (PerfPLB - mPerfPLB).^2 ) +...
                    sum( (PerfPLaB - mPerfPLaB).^2 );
            test(j,isample)=chi2;
            disp(sprintf('sigo=%f; chi2=%f;',i, chi2));    
        end
        disp(sprintf('sigo=%f; nsample=%d;',i, isample));        
    end
end

ProbSC = [ProbAL(1:600,:,:); ProbPLB(1:224,:,:); ProbPLaB(1:444,:,:); ProbPLR(1:234,:,:)];
ProbEP = [ProbAL(601:1200,:,:); ProbPLB(225:448,:,:); ProbPLaB(445:922,:,:); ProbPLR(235:462,:,:)];
ProbBZ = [ProbAL(1201:end,:,:); ProbPLB(449:end,:,:); ProbPLaB(923:end,:,:); ProbPLR(463:end,:,:)];
% llik = log(-nansum(max(-1000,log(nanmean(ProbBZ,3))),1));
% plot(llik,'-o');

%% Cleaning up the above two sections to fit lsigo just for AL (2015-04-13)
load('SCdata.mat', 'SCAL');
load('EPdata.mat', 'EPAL');
load('BZdata.mat', 'BZAL');
DIMSCAL=a_DIM(SCAL);  
DIMEPAL=a_DIM(EPAL);  
DIMBZAL=a_DIM(BZAL);  
fprintf('Finish loading data.\n');  

%% optimize for sigo 2015-05-06
AUX.Moments = [];
AUX.pvp = [];
AUX.pvd = [];

%opt= [  lsigo,    phit,   lphis, phis0,     lbx,       lbd, lapse, prpa, window, prxpa, prCommon];
opt =[log(0.1),   nan,     nan,   nan,    log(0),  log(100),     0,  1/2,      1,   1/3,      1/3];

ALTrials = 600;
nincr = 10;
% sigoVec = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2]; %minimum likely in [1e-1, 1e1]
% sigoVec = linspace(1e-1, 1e1, nincr); %minimum likely in [0.1, 2.3]
% sigoVec = linspace(0.1, 2.3, nincr); %minimum likely in [0.1, 2.3]
sigoVec = linspace(0.6, 1.5, nincr); %ans = 1.3, 0.8, 0.8;
nsample = 500;
nparti = 3;
ProbAL = nan(ALTrials, nincr, nsample, nparti);

load('expDataWPix.mat', 'SCALwPix', 'EPALwPix', 'BZALwPix');
%load('lpx.mat');
Gmap = zeros(ALTrials,25,3);

for parti = 1:3
    if(parti==1)      
        AL = SCALwPix;
    elseif(parti==2)  
        AL = EPALwPix;
    elseif(parti==3)  
        AL = BZALwPix;
    end
    for j = 1:nincr
        opt(1) = log(sigoVec(j));
        parfor isample = 1:nsample
            [~, ~, ProbAL(:,j,isample,parti)] = a_Get_Model_Perf(AL,Gmap,AUX,opt);
        end
        fprintf('parti = %d; sigo=%e\n',parti, exp(opt(1)));        
    end
end
%%
nllik = log(-nansum(max(-1000,log(nanmean(ProbAL,3))),1)); %want minimum
x = sigoVec;
plot(x, nllik(1,:,1,1),'r-o'); hold on;
plot(x, nllik(1,:,1,2),'g-o');
plot(x, nllik(1,:,1,3),'b-o');
[~, ind] = min(nllik(1,:,1,1));
fprintf('parti 1: sigo = %f;\n', x(ind));  
[~, ind] = min(nllik(1,:,1,2));
fprintf('parti 2: sigo = %f;\n', x(ind));  
[~, ind] = min(nllik(1,:,1,3));
fprintf('parti 3: sigo = %f;\n', x(ind));  

%% optimize for lbx 2015-05-06

AUX.Moments = [];
AUX.pvp = [];
AUX.pvd = [];
load('expDataWPix.mat');
load('lpx.mat');

%opt= [  lsigo,    phit,   lphis, phis0,     lbx,       lbd, lapse, prpa, window, prxpa, prCommon];
opt =[log(0.1),   nan,     nan,   nan,    log(0),  log(100),     0,  1/2,      1,   1/3,      1/3];

nincr = 10;
lbxVec = log([1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]);
% lml + lpx: nllik ~ 6
% SC: not as good >1e-5, almost equally <1e-5 --> favors 1e-8 (basically 0)
% EP: still decreasing <1e-9 --> favors 0
% BZ: not as good >1e-4, almost equally <1e-4 --> favors 1e-8 (basically 0)
% causal infer: nllik ~ 6
% SC: not as good >1e-5, almost equally <1e-5 --> favours 1e-9
% EP: a somewhat clear minimum around 1e-6 or 1e-7 --> favours 1e-8
% BZ: not as good >1e-4, almost equally <1e-4 --> favours 1e-10
% nincr = 10;
% lbxVec = linspace(-5, 0, nincr); % as a first try
nsample = 500;

for parti = 1:3
    if(parti==1)
        PLR = SCPLRwPix;  
        PLB = SCPLBwPix;  
        PLaB = SCPLaBwPix;
        opt(1) = log(1.3);  % from max likelihood fit to AL only
    elseif(parti==2)  
        PLR = EPPLRwPix;  
        PLB = EPPLBwPix;  
        PLaB = EPPLaBwPix;  
        opt(1) = log(0.8);  % from max likelihood fit to AL only
    elseif(parti==3)  
        PLR = BZPLRwPix;  
        PLB = BZPLBwPix;  
        PLaB = BZPLaBwPix;  
        opt(1) = log(0.8);  % from max likelihood fit to AL only
    end
    % for validation        
    PLR.mAnswerChoice = PROP(parti, 2).mAnswerChoice;
    PLB.mAnswerChoice = PROP(parti, 3).mAnswerChoice;
    PLaB.mAnswerChoice = PROP(parti, 4).mAnswerChoice;
    ProbPLR = nan(PLR.Trials, nincr, nsample);
    ProbPLB = nan(PLB.Trials, nincr, nsample);
    ProbPLaB = nan(PLaB.Trials, nincr, nsample);
    for j = 1:nincr
        opt(5) = lbxVec(j);
        parfor isample = 1:nsample
            Gmap = lpx(parti,1).lik;
            [~, ~, ProbPLR(:,j,isample)] = a_Get_Model_Perf(PLR,Gmap,AUX,opt);
            Gmap = lpx(parti,2).lik;
            [~, ~, ProbPLB(:,j,isample)] = a_Get_Model_Perf(PLB,Gmap,AUX,opt);
            Gmap = lpx(parti,3).lik;
            [~, ~, ProbPLaB(:,j,isample)] = a_Get_Model_Perf(PLaB,Gmap,AUX,opt);
        end
        fprintf('parti=%d; bx=%e\n',parti, exp(opt(5)));        
    end
    ProbPL = [ProbPLR(1:end,:,:); ProbPLB(1:end,:,:); ProbPLaB(1:end,:,:)];   
    if(parti==1)
        ProbSC = ProbPL;
    elseif(parti==2)  
        ProbEP = ProbPL;
    elseif(parti==3)      
        ProbBZ = ProbPL;
    end
end
%%
for parti=1:3
    if (parti==1)     sty = 'r-o'; Prob = ProbSC(:,:,:);
    elseif (parti==2) sty = 'g-o'; Prob = ProbEP(:,:,:);
    elseif (parti==3) sty = 'b-o'; Prob = ProbBZ(:,:,:);
    end
    nllik = log(-nansum(max(-1000,log(nanmean(Prob,3))),1)); %want minimum
    x = lbxVec/log(10);
    plot(x, nllik(1,:,1), sty); hold on;
    [~, ind] = min(nllik(1,:,1));
    fprintf('parti %d: bx=%e;\n', parti, x(ind));
    ylabel('-log(likelihood)');
    xlabel('log(\gamma)');
end


%% eye movement maximum likelihood
D=SCAL;
DIM=DIMSCAL;
%p = [  lsigo, phit,     lphis,  phis0,        lbx,        lbd, lapse, prpa, w];
opt=[log(0.1),    4,    log(3),      3,   log(0.1),     log(1),  0.1,   0.5, 1];
n=1; w=1;
opt_data=zeros(n,4);
x0_data=zeros(n,3);
yrM=rand(D.Trials,25*w*w);
for i=1:n
    %lb = [lbm,   km,  nm];
    lb =  [  1,  -50,   1];
    ub =  [  1,  -50, 100];
    x0 = [  1,   -50,  20];
    %x0 = rand(1,3).*(ub-lb)+lb;
    f = @(x)a_Get_BALDscoreProp(D,DIM,opt,x,yrM,1);
    options = optimset('Algorithm','interior-point');
    options = optimset(options,'Display','iter');
    options = optimset(options,'GradObj','off');
    options = optimset(options,'TolX',1e-6,'TolFun',1e-3);
    [x,fval,exitflag,output] = fmincon(f,x0,[],[],[],[],lb,ub,[],options);
    opt_data(i,:)=[x,fval];
    x0_data(i,:)=x0;
    disp(sprintf('i=%d;',i));
end
%     [              bm,              nm,              obj]
% SC: [1.13707823108535,42.0808656889251,6659755.93904558;]

%% global optimization (to compare results)
% but don't have any built-in solver now
% so check with different starting points

% Fit to only SCPLBbt; only phis0; obj=115.4093
% lb = [lsigo,     phit,  phis,  phis0, lphismax,   bd, lapse, prpa];
% lb = [log(0.1),     0,    10,      0,       10,    1,     0,  0.5];
% ub = [log(0.1),     0,    10,     20,       10,    1,  1e-3,  0.5];
% x0 = [log(0.1),     0,    10, 5.4816,       10,    1,     0,  0.5];

% % Fit to only SCPLBbt; temporal and spatial modulation; obj=74.6022
% % lb = [lsigo,     phit,  phis,  phis0, lphismax,   bd, lapse, prpa];
% % lb = [log(0.1),     0,     0,      0,        0,    1,     0,  0.5];
% % ub = [log(0.1),     1,    30,     20,       30,    1,  1e-3,  0.5];
% % x0 = [log(0.1),0.5097,10.812,1.43768,4.7290081,    1,     0,  0.5];
% % the above looks nice on the performance curve
% 
% % Compare to above, this shows spatial modulation is needed
% % lb = [lsigo,     phit,  phis,  phis0, lphismax,   bd, lapse, prpa];
% % x0 = [log(0.1),0.5097,     1,    -30,       20,    1,     0,  0.5];
% 
% % Fit to SCPLRbt+SCPLBbt+SCPLaBbt; temporal and spatial modulation; obj=519.9672
% % lb = [lsigo,     phit,  phis,  phis0, lphismax,   bd, lapse, prpa];
% % lb = [log(0.1),     0,     0,      0,        0,    1,     0,  0.5];
% % ub = [log(0.1),     1,    30,     20,       30,    1,  1e-3,  0.5];
% % x0 = [log(0.1),0.5223,0.6109,19.9970,  14.8949,    1,     0,  0.5];

% Fit to SCPLRbt+SCPLBbt+SCPLaBbt; only temporal modulation; obj=499.5042
% lb = [lsigo,     phit,  phis,  phis0, lphismax,   bd, lapse, prpa];
% lb = [log(0.01),    0,     1,   -1e5,      1e5,     0,     0,  0.5];
% ub = [log(10),      1,     1,   -1e5,      1e5,     1,     1,  0.5];
% x0 = [-1.6781,0.84638,     1,   -1e5,      1e5,0.4378,0.0361,  0.5];
% includes spatial modulation: lphismax->big, phis0->0: spatial contribution disppear

% not a bad one
% lb = [lsigo,     phit,  phis,  phis0, lphismax,   bd, lapse, prpa, w, obj];
% opt_data=[-1.67551873372160,0.822224504873503,0.776312583906946,19.5831944973255,16.9095501416046,0.444586650077876,0.001,0.5,1,501.440016598183;]

% new model

%% (UNFINISHED and OUTDATED) gradient-based local optimization -- average performance
D=a_Combine_Data(SCPLRbt,SCPLBbt,SCPLaBbt);
DIM=[DIMSCPLR,DIMSCPLB,DIMSCPLaB];
[PerfPLR, sdPLR, ~, ~] = a_Get_Human_Perf(SCPLRbt);
[PerfPLB, sdPLB, ~, ~] = a_Get_Human_Perf(SCPLBbt);
[PerfPLaB, sdPLaB, ~, ~] = a_Get_Human_Perf(SCPLaBbt);
exPerfM=[PerfPLR';PerfPLB';PerfPLaB'];
exPerfMstd=[sdPLR';sdPLB';sdPLaB'];
nc1=SCPLRbt.Trials;
nc2=SCPLBbt.Trials;
nc3=SCPLaBbt.Trials;
cond=[1, nc1;  nc1+1, nc1+nc2;  nc1+nc2+1, nc1+nc2+nc3];

n=1;
w=1;
yrM=0*rand(D.Trials,25*w*w);
opt_data_all=zeros(n,10);
x0_data_all=zeros(n,9);
for i=1:n
    %lb = [lsigo,    phit,  phis,  phis0, lphismax,   bd, lapse, prpa];
    lb = [log(0.01),    0,     0,      0,        0,    0,     0,  0.5];
    ub = [log(10),      1,     1,      0,        0,    5,     1,  0.5];
    %x0 = [log(0.1),    0,     1,      2,        1,    1,  1e-4,  0.5];
    x0 = rand(1,8).*(ub-lb)+lb;
    f = @(x)a_AllTrial_SoftLikLapse(x,w,D,DIM,yrM,exPerfM,exPerfMstd);
    options = optimset('Algorithm','interior-point');
    options = optimset(options,'Display','iter');
    options = optimset(options,'GradObj','on');
    options = optimset(options,'TolX',1e-6,'TolFun',1e-6);
    [x,fval,exitflag,output] = fmincon(f,x0,[],[],[],[],lb,ub,[],options);
    opt_data_all(i,:)=[x,w,fval];
    x0_data_all(i,:)=[x0,w];
    disp(sprintf('i=%d;',i));
end

%% plot: before vs after softmax+lapse
x=-50:50;
b1=0.096809465;
lapse=0.075709668;

likvecAL1 = Get_AllTrial_Likelihood(AL1,7,-0.632869277);
lprAL1 = log(likvecAL1)-log(1-likvecAL1); 
likvecPL5 = Get_AllTrial_Likelihood(PL5,7,-0.632869277);
lprPL5 = log(likvecPL5)-log(1-likvecPL5);

f=1./(1+exp(-x));  f_softlapse=(1-lapse)./(1+exp(-b1*x))+0.5*lapse;
yAL=1./(1+exp(-lprAL1));  yAL_softlapse=(1-lapse)./(1+exp(-b1*lprAL1))+0.5*lapse;
yPL=1./(1+exp(-lprPL5));  yPL_softlapse=(1-lapse)./(1+exp(-b1*lprPL5))+0.5*lapse;

plot([lprAL1;lprPL5],[yAL;yPL],'g.'); hold on;
plot([lprAL1;lprPL5],[yAL_softlapse;yPL_softlapse],'b.');
plot(x,f,'k',x,f_softlapse,'r');
legend('original','softmax', 'Location','NorthWest');
axis([-50,50,-0.001,1.001]);
xlabel('Log posterior ratio of choice','Fontsize',12);
ylabel('Probability','Fontsize',12);
title('Softmax + Lapse','Fontsize',12)
set(gca,'Fontsize',10);
