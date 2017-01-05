%% co-optimize sigo & hyperparameter shift (2015-05-07)

AUX.Moments = [];
AUX.pvp = [];
AUX.pvd = [];
load('expDataWPix.mat');
load('lpx.mat'); %not really used....

%opt= [  lsigo,    phit,   lphis, phis0,     lbx,       lbd, lapse, prpa, window, prxpa, prCommon, shift];
opt =[log(0.1),   nan,     nan,   nan,    log(0),  log(100),     0,  1/2,      1,   1/3,      1/3,     5];

nincr1 = 10;
nincr2 = 10;
% sigoVec = linspace(0.5, 1.5, nincr1);
sigoVec = linspace(0.1, 1, nincr1);
shiftVec = linspace(0.1, 1, nincr2);
nsample = 500;

for parti = 1:3
    if(parti==1)
        AL = SCALwPix;  
        PLR = SCPLRwPix;  
        PLB = SCPLBwPix;  
        PLaB = SCPLaBwPix;
    elseif(parti==2)  
        AL = EPALwPix;  
        PLR = EPPLRwPix;  
        PLB = EPPLBwPix;  
        PLaB = EPPLaBwPix;  
    elseif(parti==3)  
        AL = BZALwPix;  
        PLR = BZPLRwPix;  
        PLB = BZPLBwPix;  
        PLaB = BZPLaBwPix;  
    end
    ProbAL = nan(AL.Trials, nincr1, nincr2, nsample);
    ProbPLR = nan(PLR.Trials, nincr1, nincr2, nsample);
    ProbPLB = nan(PLB.Trials, nincr1, nincr2, nsample);
    ProbPLaB = nan(PLaB.Trials, nincr1, nincr2, nsample);
    LprAL = nan(AL.Trials, nincr1, nincr2, nsample);
    LprPLR = nan(PLR.Trials, nincr1, nincr2, nsample);
    LprPLB = nan(PLB.Trials, nincr1, nincr2, nsample);
    LprPLaB = nan(PLaB.Trials, nincr1, nincr2, nsample);
    for j = 1:nincr1
        opt(1) = log(sigoVec(j));
        for k = 1:nincr2
            opt(12) = shiftVec(k);
            parfor isample = 1:nsample
                Gmap = zeros(AL.Trials,25,3);
                [~, ~, ProbAL(:,j,k,isample), LprAL(:,j,k,isample)] = a_Get_Model_Perf(AL,Gmap,AUX,opt);
                Gmap = lpx(parti,1).lik;
                [~, ~, ProbPLR(:,j,k,isample), LprPLR(:,j,k,isample)] = a_Get_Model_Perf(PLR,Gmap,AUX,opt);
                Gmap = lpx(parti,2).lik;
                [~, ~, ProbPLB(:,j,k,isample), LprPLB(:,j,k,isample)] = a_Get_Model_Perf(PLB,Gmap,AUX,opt);
                Gmap = lpx(parti,3).lik;
                [~, ~, ProbPLaB(:,j,k,isample), LprPLaB(:,j,k,isample)] = a_Get_Model_Perf(PLaB,Gmap,AUX,opt);
            end
            fprintf('parti=%d; sigo=%f; shift=%f\n',parti, exp(opt(1)), opt(12));        
        end
    end
    ProbPL = [ProbAL; ProbPLR; ProbPLB; ProbPLaB];   
    LprPL = [LprAL; LprPLR; LprPLB; LprPLaB];   
    if(parti==1)
        ProbSC = ProbPL;
        LprSC = LprPL;
    elseif(parti==2)  
        ProbEP = ProbPL;
        LprEP = LprPL;
    elseif(parti==3)
        ProbBZ = ProbPL;
        LprBZ = LprPL;
    end
end
%% diagnostic: look at -log(lik) with varying sigo and shifts
% nll=[];
% for parti=1:3
%     if (parti==1)
%         ind1 = 1;% + SCALwPix.Trials + SCPLRwPix.Trials + SCPLBwPix.Trials;
%         ind2 = SCALwPix.Trials + SCPLRwPix.Trials + SCPLBwPix.Trials + SCPLaBwPix.Trials;
%         Prob = ProbSC(ind1:ind2, :,:,:);
%     elseif (parti==2) 
%         ind1 = 1;% + EPALwPix.Trials + EPPLRwPix.Trials + EPPLBwPix.Trials;
%         ind2 = EPALwPix.Trials + EPPLRwPix.Trials + EPPLBwPix.Trials + EPPLaBwPix.Trials;
%         Prob = ProbEP(ind1:ind2, :,:,:);
%     elseif (parti==3) 
%         ind1 = 1;% + BZALwPix.Trials + BZPLRwPix.Trials + BZPLBwPix.Trials;
%         ind2 = BZALwPix.Trials + BZPLRwPix.Trials + BZPLBwPix.Trials + BZPLaBwPix.Trials;
%         Prob = ProbBZ(ind1:ind2, :,:,:);
%     end
%     % plot to see best para location
%     nllik = log(-nansum(max(log(1/501),log(nanmean(Prob,4))),1)); %want min
%     nll(:,:,parti) = nllik(1,:,:,1);
%     subplot(1,3,parti);
%     imagesc(nll(:,:,parti)); %colormap('rainbow');
%     hold on;
%     axis off;    
% end

%% diagnostic: individual -log(lik) separated by rev number and trial type
% nrMax = [5,10,15,20,25];
% pp = nan(600, 4);
% pp1 = nan(600, 4);
% counting = zeros(4,5,3);
% counting1 = zeros(4,5,3);
% clf;
% for parti=1:3
%     if (parti==1)
%         AL = SCALwPix;
%         PLR = SCPLRwPix;  
%         PLB = SCPLBwPix;  
%         PLaB = SCPLaBwPix;
%         p = max(log(1/501), log(nanmean(ProbSC(:,3,3,:),4))); %max lik para
%         p1 = max(log(1/501), log(nanmean(ProbSC(:,1,5,:),4))); %good-looking para
%     elseif (parti==2) 
%         AL = EPALwPix;
%         PLR = EPPLRwPix;  
%         PLB = EPPLBwPix;  
%         PLaB = EPPLaBwPix;
%         p = max(log(1/501), log(nanmean(ProbEP(:,3,3,:),4)));
%         p1 = max(log(1/501), log(nanmean(ProbEP(:,1,5,:),4)));
%     elseif (parti==3) 
%         AL = BZALwPix;
%         PLR = BZPLRwPix;  
%         PLB = BZPLBwPix;  
%         PLaB = BZPLaBwPix;
%         p = max(log(1/501), log(nanmean(ProbBZ(:,2,1,:),4)));
%         p1 = max(log(1/501), log(nanmean(ProbBZ(:,1,5,:),4)));
%     end
%     
% %     fprintf('size(p)=%d; size(p1)=%d.\n\n', numel(p), numel(p1));
% %     figure(2)
% %     subplot(1,3,parti);
% %     plot(p,p1,'.'); hold on;
% %     plot([-8, 0.1], [-8, 0.1], 'k--');
% %     axis([-7, 0.1, -7, 0.1]);
% 
%     
%     indAL = 1:AL.Trials;
%     indPLR = 1 + AL.Trials:AL.Trials + PLR.Trials;
%     indPLB = 1 + AL.Trials + PLR.Trials:AL.Trials + PLR.Trials + PLB.Trials;
%     indPLaB = 1 + AL.Trials + PLR.Trials + PLB.Trials:AL.Trials + PLR.Trials + PLB.Trials + PLaB.Trials;
%     nrMaxAL = AL.MaxRevealingTrial;
%     nrMaxPLR = PLR.MaxRevealingTrial;
%     nrMaxPLB = PLB.MaxRevealingTrial;
%     nrMaxPLaB = PLaB.MaxRevealingTrial;
%     
%     ProbAL = p(indAL);
%     ProbPLR = p(indPLR);
%     ProbPLB = p(indPLB);
%     ProbPLaB = p(indPLaB);
%     
%     ProbAL1 = p1(indAL);
%     ProbPLR1 = p1(indPLR);
%     ProbPLB1 = p1(indPLB);
%     ProbPLaB1 = p1(indPLaB);
% 
%     bm = zeros(5,1);
%     bp = zeros(5,1);
%     
%     figure(1); %scatter plots
%     for nrInd = 1:5              
%         %figure(1); % plot marginals
%         %subplot(3,5, (parti-1)*5 + nrInd);
%        
%         p = ProbAL(nrMaxAL==nrMax(nrInd));
%         pp(1:numel(p),1) = p(:);    counting(1,nrInd,parti) = numel(p);
%         p = ProbPLR(nrMaxPLR==nrMax(nrInd));
%         pp(1:numel(p),2) = p(:);    counting(2,nrInd,parti) = numel(p);
%         p = ProbPLB(nrMaxPLB==nrMax(nrInd));
%         pp(1:numel(p),3) = p(:);    counting(3,nrInd,parti) = numel(p);   
%         p = ProbPLaB(nrMaxPLaB==nrMax(nrInd));
%         pp(1:numel(p),4) = p(:);    counting(4,nrInd,parti) = numel(p);
%         %[counts,centers] = hist(pp, 0:0.1:1);
%         %[counts,centers] = hist(pp, -6.3:0.5:0);
%         %b = bar(centers, log(counts), 'stacked', 'EdgeColor', 'k');
%         %axis([-0.1, 1.1, 0, 100]);
%         %axis([-6.5, 0.1, 0, 20]);
%         %axis square;
%         
%         subplot(1,3,parti); hold on;
%         p1 = ProbAL1(nrMaxAL==nrMax(nrInd));
%         pp1(1:numel(p1),1) = p1(:);     counting1(1,nrInd,parti) = numel(p1);
%         p1 = ProbPLR1(nrMaxPLR==nrMax(nrInd));
%         pp1(1:numel(p1),2) = p1(:);     counting1(2,nrInd,parti) = numel(p1); 
%         p1 = ProbPLB1(nrMaxPLB==nrMax(nrInd));
%         pp1(1:numel(p1),3) = p1(:);     counting1(3,nrInd,parti) = numel(p1);  
%         p1 = ProbPLaB1(nrMaxPLaB==nrMax(nrInd));
%         pp1(1:numel(p1),4) = p1(:);     counting1(4,nrInd,parti) = numel(p1);
%         if (nrInd==1)       mt = 'o';
%         elseif (nrInd==2)   mt = 'v';
%         elseif (nrInd==3)   mt = '+';
%         elseif (nrInd==4)   mt = '^';
%         elseif (nrInd==5)   mt = 'x';
%         end
%         plot(pp(:,1), pp1(:,1), mt, 'color', [1,0,0], 'MarkerSize', 5);
%         plot(pp(:,2), pp1(:,2), mt, 'color', [0,0,1], 'MarkerSize', 5);
%         plot(pp(:,3), pp1(:,3), mt, 'color', [0,0,0], 'MarkerSize', 5);
%         plot(pp(:,4), pp1(:,4), mt, 'color', 0.5*[1,1,1], 'MarkerSize', 5);
%         plot([-8, 0.1], [-8, 0.1], 'k--');
%         axis([-7, 0.1, -7, 0.1]);
%         %axis square;
%         xlabel('log(lik) @ max-lik para');
%         ylabel('log(lik) @ good-looking para');
%         
%         bm(nrInd) = sum(sum(pp<pp1));
%         bp(nrInd) = sum(sum(pp>pp1));
%     end
%         
% end

%% diagnostic: find the naughty points
% ml = log(nanmean(ProbSC(:,3,3,:),4));
% gl = log(nanmean(ProbSC(:,1,5,:),4));
% ind = find( (ml>-1 & gl<-6) );
% %ind = 599;
% figure(2)
% temp = LprSC(ind,1,5,:); 
% hist( 1./(1+exp(-temp(:))), 0:0.01:1);
% axis([0, 1, 0, 200]);
% % hist( temp(:), -3:0.1:3 );
% % axis([-3, 3, 0, 200]);

%% check gradient

lpr = LprSC(:,2,5,:);

n = 20;
dval_analytic = zeros(n,1);
dval_numeric = zeros(n,1);
testId = 1;

for i = 1:n
    %lb = [lbd,     lapse];
    lb = [log(1e-2), 0.01];
    ub = [log(10),    0.5];
    x0 = rand(1,2).*(ub-lb)+lb;
    [val, dval] = optfunc_lpr(x0, lpr);
    dval_analytic(i) = dval(testId);
    
    dx = zeros(1,2);
    dx(testId) = 1e-5;
    x0u = x0 + dx;
    x0l = x0 - dx;
    [valu, ~] = optfunc_lpr(x0u, lpr);
    [vall, ~] = optfunc_lpr(x0l, lpr);
    dval_numeric(i) = (valu-vall)/2e-5;
    
    %disp(sprintf('i=%d;',i));
end

% plots
plot(dval_analytic, dval_numeric, 'r.');
hold on;
plot([-5e2,5e2],[-5e2,5e2]);
hold off;

%% optimize lbd and lapse

%load('probParalShifts-Lpr.mat');
%load('probScaleDown-Lpr-Sigo0d6-1d0.mat');

nsigo = 10;
nshift = 10;

seedn = 1; % CHECK: should match that inside opt_lpr
paran = 2; % lbd, lapse
nparti = 3;
opt_para = nan(seedn,paran,nsigo,nshift,nparti);
opt_val = nan(seedn,nsigo,nshift,nparti);
x0_para = nan(seedn,paran,nsigo,nshift,nparti);

% CHECK: a_Model has the right parameter shift/scale method 
for parti = 1:3
    if(parti==1)
        Lpr = LprSC;
    elseif(parti==2)  
        Lpr = LprEP;
    elseif(parti==3)  
        Lpr = LprBZ;
    end    
    for isigo = 1:nsigo
        for ishift = 1:nshift
            lpr = Lpr(:,isigo,ishift,:); %change participant here
            [opt_temp, x0_para(:,:,isigo,ishift,parti)] = opt_lpr(lpr);
            opt_para(:,:,isigo,ishift,parti) = opt_temp(:,1:2);
            opt_val(:,isigo,ishift,parti) = opt_temp(:,3);
        end
        fprintf('parti %d, isigo=%d; \n', parti, isigo); 
    end
end

%% diagnostic: should check to see other values of lbd and lapse are indeed not smaller
% lpr = LprSC(:,1,4,:);
% lpr = lpr(:);
% bd = exp(-0.523);
% lapse = 0.06;
% softlik = 1./(1+exp(-bd*lpr));
% lapse_softlik = (1-lapse)*softlik + 0.5*lapse;
% val = -sum(log(lapse_softlik));  % negative log posterior of all trials
% % bd=1, lapse=0: val = 476950
% % bd=exp(-0.523), lapse=0.06: val = 447460 -- indeed smaller

%% find optimal val and parameter

% first check opt_val to see if randon starting points produce same result
for parti = 1:nparti
    val = mean(opt_val(:,:,:,parti), 1);
    val = reshape(val, nsigo, nshift); %correct reshaping
    [minval, minInd] = min(val(:));    
    fprintf('min val = %f; \n', minval);
    [indrow, indcol] = ind2sub([nsigo, nshift], minInd);
    fprintf('min val ind: %d, %d; \n', indrow, indcol);    
    valm(:,:,parti) = val;
    
    parabd(:,:,parti) = reshape(opt_para(1,1,:,:,parti),nsigo, nshift);
    fprintf('opt lbd: %f; \n', parabd(indrow, indcol, parti));
    paralapse(:,:,parti) = reshape(opt_para(1,2,:,:,parti),nsigo, nshift);
    fprintf('opt lapse: %f; \n', paralapse(indrow, indcol, parti));
end

%% put the above in excel
parti = 3;
val = valm(:,:,parti);
pbd = parabd(:,:,parti);
pla = paralapse(:,:,parti);

%% aicbic

logL = [];
npara = [];
bic = [];

% name(1) = 'sigo + bd + lapse';
logL(1) = -855.5810102 - 881.1382713 - 692.7791787;
npara(1) = 3;

% name(2) = 'sigo + bd';
logL(2) = -855.5855874 - 882.2495886 - 706.4913414;
npara(2) = 2;

% name(3) = 'sigo + para shift + bd + lapse';
logL(3) = -835.392203 - 860.3825509 - 677.0477662;
npara(4) = 4;

% name(4) = 'sigo + para shift + bd';
logL(4) = -836.0894193 - 866.9951281 - 686.1979034;
npara(4) = 3;

% name(5) = 'sigo + normal shift + bd + lapse';
logL(5) = -827.8301684 - 852.0252496 - 682.8659725;
npara(5) = 4;

% name(6) = 'sigo + normal shift + bd';
logL(6) = -830.5615503 - 855.6668444 - 696.369523;
npara(6) = 3;

% name(7) = 'sigo + scale toward origin + bd + lapse';
logL(7) = -854.4183594 - 879.4838547 - 673.3100992;
npara(7) = 4;

% name(8) = 'sigo + scale toward origin + bd';
logL(8) = -854.7913311 - 883.0949887 - 674.9329051;
npara(8) = 3;

% name(9) = 'sigo + scale toward PA + bd + lapse';
logL(9) = -854.8020578 - 884.7875716 - 657.1588315;
npara(9) = 4;

% name(10) = 'sigo + scale toward PA + bd';
logL(10) = -854.8020592 - 883.0949887 - 657.438697;
npara(10) = 3;

nobs = 1502 + 1530 + 1376;
bic = -2*logL + npara*log(nobs);
bicDiff = bic - min(bic);
[v, ind] = min(bicDiff)

%% softmax curve
xx = 0:0.01:1;
lprxx = log(xx) - log(1-xx);
% SC
lapse = 0.0242;
bd = exp(0.3664);
yy = (1-lapse)./(1+exp(-bd*lprxx)) + lapse/2;
plot(xx, yy, 'r-'); hold on;
% EP
lapse = 0.0818;
bd = exp(1.5899);
yy = (1-lapse)./(1+exp(-bd*lprxx)) + lapse/2;
plot(xx, yy, 'b-');
% BZ
lapse = 0.0410;
bd = exp(1.0806);
yy = (1-lapse)./(1+exp(-bd*lprxx)) + lapse/2;
plot(xx, yy, 'k-');
xlabel('Posterior P(c=PA)');
ylabel('Softmax P(c=PA)');



