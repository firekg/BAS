%% load data

[EPAL1,~,~] = DATAFILE_Read('data\EP20130903ALM1.DAT');
[EPAL2,~,~] = DATAFILE_Read('data\EP20130903ALM2.DAT');
[EPAL3,~,~] = DATAFILE_Read('data\EP20130903ALM3.DAT');
[EPAL4,~,~] = DATAFILE_Read('data\EP20130904ALM4.DAT');
% [EPAL5,~,~] = DATAFILE_Read('data\EP20130904ALM5.DAT');
[EPAL6,~,~] = DATAFILE_Read('data\EP20130904ALM6.DAT');
[EPAL7,~,~] = DATAFILE_Read('data\EP20130904ALM7.DAT');

[SCAL2,~,~] = DATAFILE_Read('data\SC20130903ALM2.DAT');
[SCAL3,~,~] = DATAFILE_Read('data\SC20130903ALM3.DAT');
[SCAL4,~,~] = DATAFILE_Read('data\SC20130903ALM4.DAT');
[SCAL5,~,~] = DATAFILE_Read('data\SC20130904ALM5.DAT');
[SCAL6,~,~] = DATAFILE_Read('data\SC20130904ALM6.DAT');
% [SCAL7,~,~] = DATAFILE_Read('data\SC20130904ALM7.DAT');
[SCAL8,~,~] = DATAFILE_Read('data\SC20130904ALM8.DAT');

[BZAL1,~,~] = DATAFILE_Read('data\BZ20130906ALM1.DAT');
[BZAL2,~,~] = DATAFILE_Read('data\BZ20130906ALM2.DAT');
[BZAL3,~,~] = DATAFILE_Read('data\BZ20130906ALM3.DAT');
[BZAL5,~,~] = DATAFILE_Read('data\BZ20130907ALM5.DAT'); BZAL5.FrameData=[];
[BZAL6,~,~] = DATAFILE_Read('data\BZ20130907ALM6.DAT'); BZAL6.FrameData=[];
[BZAL7,~,~] = DATAFILE_Read('data\BZ20130907ALM7.DAT'); BZAL7.FrameData=[];

[EPPL1,~,~] = DATAFILE_Read('data\EP20130905PLM1.DAT');
[EPPL2,~,~] = DATAFILE_Read('data\EP20130905PLM2.DAT');
[EPPL3,~,~] = DATAFILE_Read('data\EP20130905PLM3.DAT');
[EPPL4,~,~] = DATAFILE_Read('data\EP20130905PLM4.DAT');
[EPPL5,~,~] = DATAFILE_Read('data\EP20130905PLM5.DAT');
[EPPL6,~,~] = DATAFILE_Read('data\EP20130906PLM6.DAT');
[EPPL7,~,~] = DATAFILE_Read('data\EP20130906PLM7.DAT');
[EPPL8,~,~] = DATAFILE_Read('data\EP20130906PLM8.DAT');
[EPPL9,~,~] = DATAFILE_Read('data\EP20130906PLM9.DAT');
[EPPL10,~,~] = DATAFILE_Read('data\EP20130906PLM10.DAT');

[SCPL1,~,~] = DATAFILE_Read('data\SC20130905PLM1.DAT');
[SCPL2,~,~] = DATAFILE_Read('data\SC20130905PLM2.DAT');
[SCPL3,~,~] = DATAFILE_Read('data\SC20130905PLM3.DAT');
[SCPL4,~,~] = DATAFILE_Read('data\SC20130905PLM4.DAT');
[SCPL5,~,~] = DATAFILE_Read('data\SC20130906PLM5.DAT');
[SCPL6,~,~] = DATAFILE_Read('data\SC20130906PLM6.DAT'); %RevealID is wrong; prob didn't change it correctly during exp
[SCPL7,~,~] = DATAFILE_Read('data\SC20130906PLM7.DAT');
[SCPL8,~,~] = DATAFILE_Read('data\SC20130906PLM8.DAT');
[SCPL9,~,~] = DATAFILE_Read('data\SC20130906PLM9.DAT');
[SCPL10,~,~] = DATAFILE_Read('data\SC20130906PLM10.DAT');

[BZPL1,~,~] = DATAFILE_Read('data\BZ20130910PLM1.DAT'); %BZPL1.FrameData=[];
[BZPL2,~,~] = DATAFILE_Read('data\BZ20130910PLM2.DAT'); %BZPL2.FrameData=[];
[BZPL3,~,~] = DATAFILE_Read('data\BZ20130910PLM3.DAT'); %BZPL3.FrameData=[];
[BZPL4,~,~] = DATAFILE_Read('data\BZ20130910PLM4.DAT'); %BZPL4.FrameData=[];
[BZPL5,~,~] = DATAFILE_Read('data\BZ20130910PLM5.DAT'); %BZPL5.FrameData=[];
[BZPL6,~,~] = DATAFILE_Read('data\BZ20130910PLM6.DAT'); %BZPL6.FrameData=[];
[BZPL7,~,~] = DATAFILE_Read('data\BZ20130910PLM7.DAT'); %BZPL7.FrameData=[];
[BZPL8,~,~] = DATAFILE_Read('data\BZ20130910PLM8.DAT'); %BZPL8.FrameData=[];

% data from SY20130820 is big revealing window with w7, noise0.5
% data from SY20130822 is small revealing window with w1, noise0.1, nr 25
% data from RB20130823 is small revealing window with w1, noise0.1, nr 25; 
%           use the ones without a_, except for a_Combine_Data. Apply to
%           all previous data.
% data from SY20130828 is small revealing window with w1, noise0.1, nr 25; 
%           PL1 MaxGazeStableTrial --> MaxRevealingTrial;

%% EP data
EPAL = a_Combine_Data(EPAL1,EPAL2,EPAL3,EPAL4,EPAL6,EPAL7);
EPPL = a_Combine_Data(EPPL1,EPPL2,EPPL3,EPPL4,EPPL5,EPPL6,EPPL7,EPPL8,EPPL9,EPPL10);
EPPLR = a_Get_PL_of_RevealType(EPPL,0); EPPLRbt=a_PL_BalanceType(EPPLR);
[EPPLB,EPPLaB] = a_Get_PL_of_Bald_antiBald(EPPL); EPPLBbt=a_PL_BalanceType(EPPLB); EPPLaBbt=a_PL_BalanceType(EPPLaB);

%% SC data
SCAL = a_Combine_Data(SCAL2,SCAL3,SCAL4,SCAL5,SCAL6,SCAL8);
SCPL = a_Combine_Data(SCPL1,SCPL2,SCPL3,SCPL4,SCPL5,SCPL6,SCPL7,SCPL8,SCPL9,SCPL10);  SCPL.RevealType = EPPL.RevealType;
SCPLR = a_Get_PL_of_RevealType(SCPL,0); SCPLRbt=a_PL_BalanceType(SCPLR);
[SCPLB,SCPLaB] = a_Get_PL_of_Bald_antiBald(SCPL); SCPLBbt=a_PL_BalanceType(SCPLB); SCPLaBbt=a_PL_BalanceType(SCPLaB);

%% BZ data
BZAL = a_Combine_Data(BZAL1,BZAL2,BZAL3,BZAL5,BZAL6,BZAL7);
BZPL = a_Combine_Data(BZPL1,BZPL2,BZPL3,BZPL4,BZPL5,BZPL6,BZPL7,BZPL8);
BZPLR = a_Get_PL_of_RevealType(BZPL,0); BZPLRbt=a_PL_BalanceType(BZPLR);
[BZPLB,BZPLaB] = a_Get_PL_of_Bald_antiBald(BZPL); BZPLBbt=a_PL_BalanceType(BZPLB); BZPLaBbt=a_PL_BalanceType(BZPLaB);

%% combining active learning data
% AL = a_Combine_Data(EPAL,SCAL,BZAL);
%% extracting passive random
% PLR = a_Combine_Data(EPPLRbt,SCPLRbt,BZPLRbt);
%% extracting passive BALD and antiBALD
% PLB = a_Combine_Data(EPPLBbt,SCPLBbt,BZPLB);
% PLaB = a_Combine_Data(EPPLaBbt,SCPLaBbt,BZPLaB);

%% Performance Analysis part 1
DIMSCAL=a_DIM(SCAL);  
DIMSCPLR=a_DIM(SCPLRbt);  
DIMSCPLB=a_DIM(SCPLBbt);  
DIMSCPLaB=a_DIM(SCPLaBbt);  
%% Performance Analysis part 2
% codes involved for: opt < a_optimize_likelihood < a_AllTrial_SoftLikLapse < a_AllTrial_Likelihood < a_Model
% codes involbed for: a_Get_Model_Perf < a_Model
%opt=opt_data_all(1:9);
%opt=[log(0.1), 2, 2, 2, 0, 0.3, 0.1, 0.5, 1];
opt = [log(0.0001),   1e5,    log(1),    1e5,  log(1e-5),  log(1),    0,   0.5, 1];
%opt=[pinM(2,:),1];
 
[PerfAL, sdAL, ~, x] = a_Get_Human_Perf(SCAL);
[PerfPLR, sdPLR, ~, x] = a_Get_Human_Perf(SCPLRbt);
[PerfPLB, sdPLB, ~, x] = a_Get_Human_Perf(SCPLBbt);
[PerfPLaB, sdPLaB, ~, x] = a_Get_Human_Perf(SCPLaBbt);

[mPerfAL,msdAL,~] = a_Get_Model_Perf(SCAL,DIMSCAL,AUX,opt);
[mPerfPLR,msdPLR,~] = a_Get_Model_Perf(SCPLRbt,DIMSCPLR,AUX,opt);
[mPerfPLB,msdPLB,~] = a_Get_Model_Perf(SCPLBbt,DIMSCPLB,AUX,opt);
[mPerfPLaB,msdPLaB,~] = a_Get_Model_Perf(SCPLaBbt,DIMSCPLaB,AUX,opt);

chi2 = sum( (PerfAL - mPerfAL).^2 ) +...
        sum( (PerfPLR - mPerfPLR).^2 ) +...
        sum( (PerfPLB - mPerfPLB).^2 ) +...
        sum( (PerfPLaB - mPerfPLaB).^2 )
    
%%%  Performance Analysis part 3
% subplot('Position',[0.1,0.1,0.4,0.8]);  
subplot(1,2,1);  
errorbar(x,PerfPLB,sdPLB,'kd-','MarkerFaceColor','k');  hold on;
errorbar(x,PerfPLaB,sdPLaB,'^-','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5]);
errorbar(x,PerfPLR,sdPLR,'bs-','MarkerFaceColor','b');
errorbar(x,PerfAL,sdAL,'ro-','MarkerFaceColor','r'); hold off;
% plot(x,mPerfPLB,'k-');
% plot(x,mPerfPLaB,'-','Color',[0.5,0.5,0.5]);
% plot(x,mPerfPLR,'b-');
% plot(x,mPerfAL,'r-');
axis([0 30 0.4 1]);  xlabel('Number of revealing');  ylabel('Performance');
% legend('BALD','antiBALD','Random','AL','Location','EastOutside');
title('Subject')
% grid on

subplot(1,2,2);
errorbar(x,mPerfPLB,msdPLB,'kd-','MarkerFaceColor','k');  hold on;
errorbar(x,mPerfPLaB,msdPLaB,'^-','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5]);
errorbar(x,mPerfPLR,msdPLR,'bs-','MarkerFaceColor','b');
errorbar(x,mPerfAL,msdAL,'ro-','MarkerFaceColor','r'); hold off;
axis([0 30 0.4 1]);  xlabel('Number of revealing');  ylabel('Performance');
% leg=legend('BALD','anti-BALD','random','AL');
% set(leg, 'Position',[0.83,0.5,0.15,0.15]);
% set(leg, 'Location','NorthWest');
title('Model')
% grid on

%% BALD score Analysis part 1
pv=[-0.5427, 0.129, 0, 0, 0, 0.4447, 0.09, 0.5, 1];
pw=[0.959858337511111,-48.7463898814911,48.1723964639968];
yrM=0*rand(1000,25);
SCPLBsoftm=a_Get_BALDnoisy(SCAL,DIMSCAL,yrM,pv,pw,1);
SCPLBhardm=a_Get_BALDnoisy(SCAL,DIMSCAL,yrM,pv,pw,2);
pw(3)=1;
SCPLBsoft=a_Get_BALDnoisy(SCAL,DIMSCAL,yrM,pv,pw,1);
SCPLBhard=a_Get_BALDnoisy(SCAL,DIMSCAL,yrM,pv,pw,2);

%% Revealing position Analysis Part 2
%pv = [  lsigo, phit,     lphis,  phis0,        lbx,        lbd, lapse, prpa, w];
pv= [log(1e-5),    4,    log(3),      3,   log(0.1),     log(1),  0.1,   0.5, 1];
%pw = [    lbm,   km,  nm];
pw = [  log(1),  0.1,  20];
%ESCAL = a_Get_BALDscoreProp(SCAL,DIMSCAL,pv,pw,yrM,2);
ESCPLR = a_Get_BALDscoreProp(SCPLRbt,DIMSCPLR,pv,pw,yrM,2);
%ESCPLB = a_Get_BALDscoreProp(SCPLBbt,DIMSCPLB,pv,pw,yrM,2);

%% BALD score Analysis part 3: BALDscore percentile
M = ESCAL.PTILExs;
ALptile_m = nanmean(M);  
ALptile_v = nanvar(M);  
ALptile_n = sum(~isnan(M));

M = ESCPLR.PTILExs;
PLRptile_m = nanmean(M);  
PLRptile_v = nanvar(M);  
PLRptile_n = sum(~isnan(M));

M = ESCPLB.PTILExs;
PLBptile_m = nanmean(M);  
PLBptile_v = nanvar(M);  
PLBptile_n = sum(~isnan(M));

% M = SCPLBhard.PTILExs;
% PLBptile_m = nanmean(M);  
% PLBptile_v = nanvar(M);  
% PLBptile_n = sum(~isnan(M));
% 
% M = SCPLBhardm.PTILExs;
% PLBmptile_m = nanmean(M);  
% PLBmptile_v = nanvar(M);  
% PLBmptile_n = sum(~isnan(M));
% 
% M = SCPLBsoft.PTILExs;
% PLBsoftptile_m = nanmean(M);  
% PLBsoftptile_v = nanvar(M);  
% PLBsoftptile_n = sum(~isnan(M));
% 
% M = SCPLBsoftm.PTILExs;
% PLBsoftmptile_m = nanmean(M);  
% PLBsoftmptile_v = nanvar(M);  
% PLBsoftmptile_n = sum(~isnan(M));

Bx=1:25;
figure(1)
errorbar(Bx, PLBptile_m, sqrt(PLBptile_v./PLBptile_n), 'kd-','MarkerFaceColor','k'); hold on;
% errorbar(Bx, PLBmptile_m, sqrt(PLBmptile_v./PLBmptile_n), 'd-','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5]);
% errorbar(Bx, PLBsoftptile_m, sqrt(PLBsoftptile_v./PLBsoftptile_n), 'g-^','MarkerFaceColor','g');
% errorbar(Bx, PLBsoftmptile_m, sqrt(PLBsoftmptile_v./PLBsoftmptile_n), '-^','Color',[0,0.5,0],'MarkerFaceColor',[0,0.5,0]);
errorbar(Bx, ALptile_m, sqrt(ALptile_v./ALptile_n), 'r-o','MarkerFaceColor','r');
errorbar(Bx, PLRptile_m, sqrt(PLRptile_v./PLRptile_n), 'b-s','MarkerFaceColor','b'); hold off;

% a_ErrorShade(Bx, PLBptile_m, sqrt(PLBptile_v./PLBptile_n), [0,0,0],0.4); hold on;
% a_ErrorShade(Bx, ALptile_m, sqrt(ALptile_v./ALptile_n), [1,0,0],0.4);
% a_ErrorShade(Bx, PLBsoftptile_m, sqrt(PLBsoftptile_v./PLBsoftptile_n), [0,0.5,0],0.4);
% a_ErrorShade(Bx, PLRptile_m, sqrt(PLRptile_v./PLRptile_n), [0,0,1],0.4);
% h=zeros(1,4);
% h(1)=plot(Bx, PLBptile_m, 'k-', 'DisplayName','BALD'); hold on;
% h(2)=plot(Bx, ALptile_m,  'r-', 'DisplayName','BALD');
% h(3)=plot(Bx, PLBsoftptile_m, '-', 'Color', [0,0.5,0], 'DisplayName','BALD');
% h(4)=plot(Bx, PLRptile_m, 'b-', 'DisplayName','BALD');
% legend(h,'Hardmax', 'AL', 'Softmax', 'random', 'Location','SouthWest'); legend('boxoff');

axis([0 25 30 100]);  %axis square;  %grid on;  
set(gca,'Xtick',0:5:25);  xlabel('Revealing number', 'FontSize', 10);
set(gca,'Ytick',30:10:100);  ylabel('BALD percentile', 'FontSize', 10);
% legend('Hardmax', 'Hardmax *', 'Softmax', 'Softmax *', 'AL', 'random', 'Location','SouthWest'); legend('boxoff');

% title('Optimality measure: BALD percentile')

% disp(sprintf('chi2=%f',nansum((ALptile_m-PLBsoftptile_m).^2)));

%% BALD score Analysis part 4: cumulative BALD score
ESCAL.BALDxs(:,1)=0;  M=cumsum(ESCAL.BALDxs,2); M= -M/log(0.5);
ALcums_m = nanmean(M);
ALcums_v = nanvar(M);
ALcums_n = sum(~isnan(M));

ESCPLR.BALDxs(:,1)=0;  M=cumsum(ESCPLR.BALDxs,2); M= -M/log(0.5);
PLRcums_m = nanmean(M);
PLRcums_v = nanvar(M);
PLRcums_n = sum(~isnan(M));

% SCPLBhard.BALDxs(:,1)=0;  M=cumsum(SCPLBhard.BALDxs,2); M= -M/log(0.5);
% PLBcums_m = nanmean(M);
% PLBcums_v = nanvar(M);
% PLBcums_n = sum(~isnan(M));
% 
% SCPLBhardm.BALDxs(:,1)=0;  M=cumsum(SCPLBhardm.BALDxs,2); M= -M/log(0.5);
% PLBmcums_m = nanmean(M);
% PLBmcums_v = nanvar(M);
% PLBmcums_n = sum(~isnan(M));
% 
% SCPLBsoft.BALDxs(:,1)=0;  M=cumsum(SCPLBsoft.BALDxs,2); M= -M/log(0.5);
% PLBsoftcums_m = nanmean(M);
% PLBsoftcums_v = nanvar(M);
% PLBsoftcums_n = sum(~isnan(M));
% 
% SCPLBsoftm.BALDxs(:,1)=0;  M=cumsum(SCPLBsoftm.BALDxs,2); M= -M/log(0.5);
% PLBsoftmcums_m = nanmean(M);
% PLBsoftmcums_v = nanvar(M);
% PLBsoftmcums_n = sum(~isnan(M));

Bx=1:25;
% figure(2)
% errorbar(Bx, PLBcums_m, sqrt(PLBcums_v./PLBcums_n), 'k-d','MarkerFaceColor','k'); hold on;
% errorbar(Bx, PLBmcums_m, sqrt(PLBmcums_v./PLBmcums_n), 'd-','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5]);
% errorbar(Bx, PLBsoftcums_m, sqrt(PLBsoftcums_v./PLBsoftcums_n), 'g-^','MarkerFaceColor','g');
% errorbar(Bx, PLBsoftmcums_m, sqrt(PLBsoftmcums_v./PLBsoftmcums_n), '-^','Color',[0,0.5,0],'MarkerFaceColor',[0,0.5,0]);
errorbar(Bx, ALcums_m, sqrt(ALcums_v./ALcums_n), 'r-o','MarkerFaceColor','r'); hold on;
errorbar(Bx, PLRcums_m, sqrt(PLRcums_v./PLRcums_n), 'b-s','MarkerFaceColor','b'); hold off;

axis([0 25 0 1]);  %axis square;  %grid on; 
set(gca,'Ytick',0:0.2:1); ylabel('Cumulative BALD score', 'FontSize', 10);
set(gca,'Xtick',0:5:25);  xlabel('Revealing number', 'FontSize', 10);
ylabel('Cumulative BALD score', 'FontSize', 10);  %set(gca,'Ytick',0:10:100);
% legend('Hardmax', 'Hardmax *', 'Softmax', 'Softmax *', 'AL', 'random', 'Location','NorthWest'); legend('boxoff');
grid on;

% drawnow;
% title('Optimality measure: BALD percentile')

%%% BALD score Analysis part 5: entropy of whole system
m=ESCAL.mla; M = -log(0.5) +m.*log(m) +(1-m).*log(1-m); M= -M/log(0.5);
ALentr_m = nanmean(M);
ALentr_v = nanvar(M);
ALentr_n = sum(~isnan(M));

m=ESCPLR.mla; M = -log(0.5) +m.*log(m) +(1-m).*log(1-m); M= -M/log(0.5);
PLRentr_m = nanmean(M);
PLRentr_v = nanvar(M);
PLRentr_n = sum(~isnan(M));

% m=SCPLBhard.mla; M = -log(0.5) +m.*log(m) +(1-m).*log(1-m); M= -M/log(0.5);
% PLBentr_m = nanmean(M);
% PLBentr_v = nanvar(M);
% PLBentr_n = sum(~isnan(M));
% 
% m=SCPLBhardm.mla; M = -log(0.5) +m.*log(m) +(1-m).*log(1-m); M= -M/log(0.5);
% PLBmentr_m = nanmean(M);
% PLBmentr_v = nanvar(M);
% PLBmentr_n = sum(~isnan(M));
% 
% m=SCPLBsoft.mla; M = -log(0.5) +m.*log(m) +(1-m).*log(1-m); M= -M/log(0.5);
% PLBsoftentr_m = nanmean(M);
% PLBsoftentr_v = nanvar(M);
% PLBsoftentr_n = sum(~isnan(M));
% 
% m=SCPLBsoftm.mla; M = -log(0.5) +m.*log(m) +(1-m).*log(1-m); M= -M/log(0.5);
% PLBsoftmentr_m = nanmean(M);
% PLBsoftmentr_v = nanvar(M);
% PLBsoftmentr_n = sum(~isnan(M));

Bx=1:25;
hold on;% figure(3);
% errorbar(Bx, PLBentr_m, sqrt(PLBentr_v./PLBentr_n), 'k-d','MarkerFaceColor','k'); hold on;
% errorbar(Bx, PLBmentr_m, sqrt(PLBmentr_v./PLBmentr_n), 'd-','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5]);
% errorbar(Bx, PLBsoftentr_m, sqrt(PLBsoftentr_v./PLBsoftentr_n), 'g-^','MarkerFaceColor','g');
% errorbar(Bx, PLBsoftmentr_m, sqrt(PLBsoftmentr_v./PLBsoftmentr_n), '-^','Color',[0,0.5,0],'MarkerFaceColor',[0,0.5,0]);
errorbar(Bx, ALentr_m, sqrt(ALentr_v./ALentr_n), 'r-o','MarkerFaceColor','r'); hold on;
errorbar(Bx, PLRentr_m, sqrt(PLRentr_v./PLRentr_n), 'b-s','MarkerFaceColor','b'); hold off;

axis([0 25 0 1]);  %axis square;  %grid on;  
set(gca,'Ytick',0:0.2:1);  ylabel('Entropy gain', 'FontSize', 10);
set(gca,'Xtick',0:5:25);  xlabel('Revealing number', 'FontSize', 10);
ylabel('Entropy gain', 'FontSize', 10); %set(gca,'Ytick',0:0.1:0.7);
% legend('Hardmax', 'Hardmax *', 'Softmax', 'Softmax *', 'AL', 'random', 'Location','NorthWest'); legend('boxoff');
grid on;

% drawnow;
% title('Optimality measure: BALD percentile')

%% Revealing position Analysis Part 1
SCALpa=a_Get_ImageID(SCAL,20);
SCALsh=a_Get_ImageID(SCAL,6);
SCALsv=a_Get_ImageID(SCAL,30);

MAL = SCPLBsoft;
MALpa=a_Get_ImageID(SCPLBsoft,20);
MALsh=a_Get_ImageID(SCPLBsoft,6);
MALsv=a_Get_ImageID(SCPLBsoft,30);

plot(MALpa.RevealPosX,MALpa.RevealPosY,'g.'); hold on;
plot(MALsh.RevealPosX,MALsh.RevealPosY,'r.');
plot(MALsv.RevealPosX,MALsv.RevealPosY,'b.'); hold off;

% MALpa.Trials=400;
% MALpa.RevealPosX = BALDrevTpa.W_BALDx';
% MALpa.RevealPosY = BALDrevTpa.W_BALDy';
% MALpa.MaxRevealingTrial = ones(400,1)*25;
% 
% MALsh.Trials=400;
% MALsh.RevealPosX = BALDrevTsh.W_BALDx';
% MALsh.RevealPosY = BALDrevTsh.W_BALDy';
% MALsh.MaxRevealingTrial = ones(400,1)*25;
% 
% MALsv.Trials=400;
% MALsv.RevealPosX = BALDrevTsv.W_BALDx';
% MALsv.RevealPosY = BALDrevTsv.W_BALDy';
% MALsv.MaxRevealingTrial = ones(400,1)*25;
% 
% MAL = a_Combine_Data(MALpa,MALsh,MALsv);

%% Revealing position Analysis Part 2

MrevMap_all = a_Get_RevealMap(MAL);
MrevMap_pa  = a_Get_RevealMap(MALpa);
MrevMap_sh  = a_Get_RevealMap(MALsh);
MrevMap_sv  = a_Get_RevealMap(MALsv);

HrevMap_all = a_Get_RevealMap(SCAL);
HrevMap_pa  = a_Get_RevealMap(SCALpa);
HrevMap_sh  = a_Get_RevealMap(SCALsh);
HrevMap_sv  = a_Get_RevealMap(SCALsv);
%  Revealing position Analysis Part 3

dmrmpa=MrevMap_pa-MrevMap_all;  dmrmsh=MrevMap_sh-MrevMap_all; dmrmsv=MrevMap_sv-MrevMap_all;
dhrmpa=HrevMap_pa-HrevMap_all;  dhrmsh=HrevMap_sh-HrevMap_all; dhrmsv=HrevMap_sv-HrevMap_all;

vmax = max([max(max(dmrmpa)), max(max(dmrmsh)), max(max(dmrmsv)),...
            max(max(dhrmpa)), max(max(dhrmsh)), max(max(dhrmsv))]);
vmin = min([min(min(dmrmpa)), min(min(dmrmsh)), min(min(dmrmsv)),...
            min(min(dhrmpa)), min(min(dhrmsh)), min(min(dhrmsv))]);
dv = vmax-vmin;  % should be in the range 1 to length(get(gcf,'Colormap'))

%%  Revealing position Analysis Part 4

figure('Position',[100,200,1000,600]);
colormap(hot);  nc = length(get(gcf,'Colormap'));

subplot(6,9,[1:3, 9+1:9+3, 18+1:18+3]); image((dhrmpa-vmin)/dv*(nc-1)+1); axis square; axis off; title('Patchy');
subplot(6,9,[4:6, 9+4:9+6, 18+4:18+6]); image((dhrmsh-vmin)/dv*(nc-1)+1); axis square; axis off; title('Horizontal stripy');
subplot(6,9,[7:9, 9+7:9+9, 18+7:18+9]); image((dhrmsv-vmin)/dv*(nc-1)+1); axis square; axis off; title('Vertical stripy');

subplot(6,9,[27+1:27+3, 36+1:36+3, 45+1:45+3]); image((dmrmpa-vmin)/dv*(nc-1)+1); axis square; axis off;
subplot(6,9,[27+4:27+6, 36+4:36+6, 45+4:45+6]); image((dmrmsh-vmin)/dv*(nc-1)+1); axis square; axis off;
subplot(6,9,[27+7:27+9, 36+7:36+9, 45+7:45+9]); image((dmrmsv-vmin)/dv*(nc-1)+1); axis square; axis off;
% colorbar('Ytick', [1,32,64],'YTickLabel',{'-1','0','1'});
% text(0,0,'Subject');
% text(10,10,'Model');

%% Isoperimetric ratio analysis
AL=SCPLB;

ALpa=a_Get_ImageID(AL,20);
ALsh=a_Get_ImageID(AL,6);
ALsv=a_Get_ImageID(AL,30);

isoperiRpa=a_isoperiR(ALpa);
isoperiRsh=a_isoperiR(ALsh);
isoperiRsv=a_isoperiR(ALsv);

pam=nanmean(isoperiRpa);  pam(:,1:2)=[];
shm=nanmean(isoperiRsh);  shm(:,1:2)=[];
svm=nanmean(isoperiRsv);  svm(:,1:2)=[];

pav=nanvar(isoperiRpa)./sum((isnan(isoperiRpa)==0));  pav(:,1:2)=[];
shv=nanvar(isoperiRsh)./sum((isnan(isoperiRpa)==0));  shv(:,1:2)=[];
svv=nanvar(isoperiRsv)./sum((isnan(isoperiRpa)==0));  svv(:,1:2)=[];

clf;

f = [shm'+2*sqrt(shv'); flipdim(shm'-2*sqrt(shv'),1)];
z = (3:25)';  fill([z; flipdim(z,1)], f, [4 4 9]/10,'FaceAlpha', 0.6, 'LineStyle','None'); hold on; plot(z,shm, 'b-')

f = [svm'+2*sqrt(svv'); flipdim(svm'-2*sqrt(svv'),1)];
z = (3:25)';  fill([z; flipdim(z,1)], f, [9 4 4]/10,'FaceAlpha', 0.6, 'LineStyle','None'); hold on; plot(z, svm, 'r-')

f = [pam'+2*sqrt(pav'); flipdim(pam'-2*sqrt(pav'),1)];
z = (3:25)';  fill([z; flipdim(z,1)], f, [6 6 6]/10,'FaceAlpha', 0.6, 'LineStyle','None'); hold on; plot(z, pam, 'k-')

%% Histogramming as a function of revealing number

% BALDpa.RevealPosX=BALDpa.RevealPosX';
% BALDpa.RevealPosY=BALDpa.RevealPosY';
% BALDsh.RevealPosX=BALDsh.RevealPosX';
% BALDsh.RevealPosY=BALDsh.RevealPosY';
% BALDsv.RevealPosX=BALDsv.RevealPosX';
% BALDsv.RevealPosY=BALDsv.RevealPosY';

isoperiRpa=a_isoperiR(BALDpa);
isoperiRsh=a_isoperiR(BALDsh);
isoperiRsv=a_isoperiR(BALDsv);

bin = 0.1:0.1:0.9;

for i=1:5
    subplot(5,1,i)
    PAHist  = hist(isoperiRpa(:,i*5),bin);
    bar(bin, PAHist, 'r');    
    xlabel('Isoperimetric ratio', 'FontSize', 10);
    ylabel('Counts', 'FontSize', 10);    
end
figure(2)
for i=1:5
    subplot(5,1,i)
    Sist = hist([isoperiRsh(:,i*5);isoperiRsv(:,i*5)],bin);
    bar(bin, Sist, 'b');
    xlabel('Isoperimetric ratio', 'FontSize', 10);
    ylabel('Counts', 'FontSize', 10);
end


