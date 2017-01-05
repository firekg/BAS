%% see if there are parameters that mimic the qualitative performance
% NOTE: a_Model has changed: so the following may not work
% a_Get_Model_Perf changed as well

[PerfAL, sdAL, ~, x] = a_Get_Human_Perf(SCAL);
[PerfPLR, sdPLR, ~, x] = a_Get_Human_Perf(SCPLRbt);
[PerfPLB, sdPLB, ~, x] = a_Get_Human_Perf(SCPLBbt);
[PerfPLaB, sdPLaB, ~, x] = a_Get_Human_Perf(SCPLaBbt);

%p = [lsigo,    phit,  phis,  phis0,    bx,   bd, lapse, prpa,  w];
opt= [log(0.1),    1,    10,      1,   0.1,    1,   0.1,  0.5,  1];
SCAL.RevealType = zeros(SCAL.Trials,1);

AUX.Moments=Moments; %Moments should be defined
AUX.pvp=[]; AUX.pvd=[];

[mPerfAL,msdAL,~] = a_Get_Model_Perf(SCAL,DIMSCAL,AUX,opt);
[mPerfPLR,msdPLR,~] = a_Get_Model_Perf(SCPLRbt,DIMSCPLR,AUX,opt);
[mPerfPLB,msdPLB,~] = a_Get_Model_Perf(SCPLBbt,DIMSCPLB,AUX,opt);
[mPerfPLaB,msdPLaB,~] = a_Get_Model_Perf(SCPLaBbt,DIMSCPLaB,AUX,opt);

% plots
subplot(1,2,1);  
errorbar(x,PerfPLB,sdPLB,'kd-','MarkerFaceColor','k');  hold on;
errorbar(x,PerfPLaB,sdPLaB,'^-','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5]);
errorbar(x,PerfPLR,sdPLR,'bs-','MarkerFaceColor','b');
errorbar(x,PerfAL,sdAL,'ro-','MarkerFaceColor','r'); hold off;
axis([0 30 0.4 1]);  xlabel('Number of revealing');  ylabel('Performance');
title('Subject')

subplot(1,2,2);
errorbar(x,mPerfPLB,msdPLB,'kd-','MarkerFaceColor','k');  hold on;
errorbar(x,mPerfPLaB,msdPLaB,'^-','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5]);
errorbar(x,mPerfPLR,msdPLR,'bs-','MarkerFaceColor','b');
errorbar(x,mPerfAL,msdAL,'ro-','MarkerFaceColor','r'); hold off;
axis([0 30 0.4 1]);  xlabel('Number of revealing');  ylabel('Performance');
title('Model')

%  Hard to find a parameter set that does well....the fit does about as
%  well as possible. May need other forms of fall off.

%% validation loop
pinM=zeros(10,8); poutM=pinM;
for ii=1:10
%%% generate data: input=D,DIM,perception parameters; output=answer choices
D=a_Combine_Data(SCPLRbt,SCPLBbt,SCPLaBbt);
DIM=[DIMSCPLR,DIMSCPLB,DIMSCPLaB];
yrM=randn(1000,25*9);
npix=770;  ImageSize=20;

%pvp=[w,lsigo,phit,phis,phis0,bx];
%pvd=[bd,lapse,prpa];
%p = [lsigo,     phit,     lphis,  phis0,        lbx,        lbd, lapse, prpa];
lb = [log(0.1),   0.5,  log(0.1),    0.5,  log(1e-3),  log(1e-2),  0.01,  0.5];
ub = [log(1.0),     5,   log(20),     10,     log(1),    log(10),   0.5,  0.5];
p_in = rand(1,8).*(ub-lb)+lb;
pvp= [1, p_in(1:5)];
pvd= p_in(6:8);
% also make sure Moments is OK

prpa=pvd(3);
lapse=pvd(2);
sigd=0.15;
for trial=1:D.Trials    
    yimage = DIM(trial).yimage;
    nv = D.MaxRevealingTrial(trial);
    xDx = 1+(npix-1)*(D.RevealPosX(trial,1:nv)/ImageSize+0.5);  % cm to pixel: inverting last two lines of a_GET_BALD
    xDy = 1-(npix-1)*(D.RevealPosY(trial,1:nv)/ImageSize-0.5);
    xD = [xDx',xDy'];
    yrin=yrM(trial,:)';
%     %conditioned on perceived value (bd doesn't play a role here)
%     [~,dlml,~]=a_Model(xD, yimage, pvp, Moments, pvd, yrin);
%     lpr = dlml(1)-log(exp(dlml(2))+exp(dlml(3))) + log(2*prpa)-log(1-prpa);
%     %simulated and overwrite AnswerChoice
%     if(lpr+sigd*randn(1)>0) 
%         if (rand(1)>lapse/2) %add in lapse rate
%             D.AnswerChoice(trial,1)=1;
%         else
%             D.AnswerChoice(trial,1)=2;
%         end
%     else
%         if (rand(1)>lapse/2)
%             D.AnswerChoice(trial,1)=2;
%         else
%             D.AnswerChoice(trial,1)=1;
%         end
%     end
    %conditioned on presented value and softmaxed
    % a_Model has changed, so this may not work now
    [lapse_softlik,~] = a_Model(xD, yimage, pvp, Moments, pvd, yrin, [], 1);    
    if(rand(1)<lapse_softlik) % use uniform random!
        D.AnswerChoice(trial,1)=1;
    else
        D.AnswerChoice(trial,1)=2;
    end
end

%%% fit answer choices and see if we get back the input parameters
n=1; w=1;

AUX.Moments=Moments; %Moments should be defined
AUX.pvp=[]; AUX.pvd=[];

opt_out=zeros(n,10);
x0_out=zeros(n,9);
for i=1:n
    %p = [lsigo,     phit,     lphis,  phis0,        lbx,        lbd, lapse, prpa];
    lb = [log(0.1),   0.5,  log(0.1),    0.5,  log(1e-3),  log(1e-2),  0.01,  0.5];
    ub = [log(1.0),     5,   log(20),     10,     log(1),    log(10),   0.5,  0.5];
    %x0 =[log(0.1),    0,     1,      2,     1,    1,  1e-4, 0.5];
    x0 = rand(1,8).*(ub-lb)+lb;
    f = @(x)a_AllTrial_SoftLikLapse(x,w,D,DIM,AUX);
    options = optimset('Algorithm','interior-point');
    options = optimset(options,'Display','iter');
    options = optimset(options,'GradObj','on');
    options = optimset(options,'TolX',1e-6,'TolFun',1e-6);
    [x,fval,exitflag,output] = fmincon(f,x0,[],[],[],[],lb,ub,[],options);
    opt_out(i,:)=[x,w,fval];
    x0_out(i,:)=[x0,w];
    %disp(sprintf('i=%d;',i));
end
%%%
pinM(ii,:)=p_in;
poutM(ii,:)=x;
disp(sprintf('ii=%d;',ii));
end

%% plot the fits
DR=SCPLRbt;  DRIM=DIMSCPLR;  c1=SCPLRbt.Trials;  DR.AnswerChoice=D.AnswerChoice(1:c1);
DB=SCPLBbt;  DBIM=DIMSCPLB;  c2=SCPLBbt.Trials;  DB.AnswerChoice=D.AnswerChoice(c1+1:c1+c2);
DaB=SCPLaBbt;  DaBIM=DIMSCPLaB;  c3=SCPLaBbt.Trials;  DaB.AnswerChoice=D.AnswerChoice(c1+c2+1:c1+c2+c3);

mx=5:5:25;

[PerfPLR, sdPLR, ~, ~] = a_Get_Human_Perf(DR);
[PerfPLB, sdPLB, ~, ~] = a_Get_Human_Perf(DB);
[PerfPLaB, sdPLaB, ~, ~] = a_Get_Human_Perf(DaB);

subplot(1,2,1);  
errorbar(mx,PerfPLB,sdPLB,'kd-','MarkerFaceColor','k');  hold on;
errorbar(mx,PerfPLaB,sdPLaB,'^-','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5]);
errorbar(mx,PerfPLR,sdPLR,'bs-','MarkerFaceColor','b'); hold off;
axis([0 30 0.4 1]);  xlabel('Number of revealing');  ylabel('Performance');
title('Input simulation')

opt = opt_out(1:9);
[mPerfPLR,msdPLR,~] = a_Get_Model_Perf(DR,DRIM,AUX,opt);
[mPerfPLB,msdPLB,~] = a_Get_Model_Perf(DB,DBIM,AUX,opt);
[mPerfPLaB,msdPLaB,~] = a_Get_Model_Perf(DaB,DaBIM,AUX,opt);

subplot(1,2,2);
errorbar(mx,mPerfPLB,msdPLB,'kd-','MarkerFaceColor','k');  hold on;
errorbar(mx,mPerfPLaB,msdPLaB,'^-','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5]);
errorbar(mx,mPerfPLR,msdPLR,'bs-','MarkerFaceColor','b'); hold off;
axis([0 30 0.4 1]);  xlabel('Number of revealing');  ylabel('Performance');
title('Fit')

%% Record sigma_p
D=a_Combine_Data(SCPLRbt,SCPLBbt,SCPLaBbt);
DIM=[DIMSCPLR,DIMSCPLB,DIMSCPLaB];
yrM=zeros(1000,25*9); % do not need yrM for recording sigma_p
npix=770;  ImageSize=20;  nexp=10;
Rsigp_simu=nan(D.Trials*nexp,25);
%Rsigp_fit=nan(D.Trials*nexp,25);

for i=1:nexp
    %p = [lsigo,    phit,  phis,  phis0,   bx,   bd, lapse, prpa];
    %p_r= [ -0.5,       5,    15,    0.1,  0.2,    0,   0.1,  0.5]; %test
    p_r = pinM(i,:); %pinM contains simu parameter values
    %p_r = poutM(i,:); %poutM contains fit parameter values
    pvp= [1, p_r(1:5)];
    pvd= p_r(6:8);
    % also make sure Moments is OK

    prpa=pvd(3);
    lapse=pvd(2);
    sigd=0.15;
    for trial=1:D.Trials    
        yimage = DIM(trial).yimage;
        nv = D.MaxRevealingTrial(trial);
        xDx = 1+(npix-1)*(D.RevealPosX(trial,1:nv)/ImageSize+0.5);  % cm to pixel: inverting last two lines of a_GET_BALD
        xDy = 1-(npix-1)*(D.RevealPosY(trial,1:nv)/ImageSize-0.5);
        xD = [xDx',xDy'];
        yrin=yrM(trial,:)';
        % a_Model has changed, so this may not work now
        [~,~,sigma_p]=a_Model(xD, yimage, pvp, Moments, pvd, yrin, [], 1);
        Rsigp_simu((i-1)*D.Trials+trial,1:nv)=sigma_p;
        %Rsigp_fit((i-1)*D.Trials+trial,1:nv)=sigma_p;
    end
    disp(sprintf('i=%d;',i));
end

%% Plot validation curves:

figure(1);
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
set(gcf,'Position',[0,0,20,8]);
set(gcf, 'PaperSize', [20,8]);
set(gcf, 'PaperPosition', [0,0,20,8]);
set(gcf, 'Color', [1,1,1]);

ncol=5;

subplot(2,ncol,ncol-2);
lo=-2.2;
up=0;
plot(pinM(:,1),poutM(:,1),'ro', [lo,up],[lo,up],'-');
axis([lo,up,lo,up]);  axis square;
set(gca,'FontSize',9);
xlabel('Input value');  ylabel('Fit value');
set(gca,'Xtick',-2:1:0,'XTickLabel',{'-2','-1','0'});
set(gca,'Ytick',-2:1:0,'YTickLabel',{'-2','-1','0'});
title('log(\sigma_o)')

subplot(2,ncol,2*ncol-2);
lo=3;
up=10;
plot(pinM(:,2),poutM(:,2),'ro', [lo,up],[lo,up],'-');
axis([lo,up,lo,up]);  axis square;
set(gca,'FontSize',9);
xlabel('Input value (s)');  ylabel('Fit value (s)');
set(gca,'Xtick',[5,10],'XTickLabel',{'5','10'});
set(gca,'Ytick',[5,10],'YTickLabel',{'5','10'});
title('\tau')

% subplot(2,ncol,ncol-1);
% lo=0;
% up=20;
% plot(pinM(:,3),poutM(:,3),'ro', [lo,up],[lo,up],'-');
% axis([lo,up,lo,up]);  axis square;
% set(gca,'FontSize',9);
% xlabel('Input value (cm^{-1})');  ylabel('Fit value (cm^{-1})');
% set(gca,'Xtick',0:10:20,'XTickLabel',{'0','10','20'});
% set(gca,'Ytick',0:10:20,'YTickLabel',{'0','10','20'});
% title('\lambda_p')

subplot(2,ncol,ncol-1);
lo=0;
up=2;
plot(pinM(:,6),poutM(:,6),'ro', [lo,up],[lo,up],'-');
axis([lo,up,lo,up]);  axis square;
set(gca,'FontSize',9);
xlabel('Input value');  ylabel('Fit value');
set(gca,'Xtick',0:1:2,'XTickLabel',{'0','1','2'});
set(gca,'Ytick',0:1:2,'YTickLabel',{'0','1','2'});
title('\beta_d')

subplot(2,ncol,2*ncol-1);
lo=0;
up=10;
plot(pinM(:,4),poutM(:,4),'ro', [lo,up],[lo,up],'-');
axis([lo,up,lo,up]);  axis square;
set(gca,'FontSize',9);
xlabel('Input value (cm)');  ylabel('Fit value (cm)');
set(gca,'Xtick',0:5:10,'XTickLabel',{'0','5','10'});
set(gca,'Ytick',0:5:10,'YTickLabel',{'0','5','10'});
title('\Delta_o')

subplot(2,ncol,ncol);
lo=0;
up=1;
plot(pinM(:,5),poutM(:,5),'ro', [lo,up],[lo,up],'-');
axis([lo,up,lo,up]);  axis square;
set(gca,'FontSize',9);
xlabel('Input value');  ylabel('Fit value');
set(gca,'Xtick',0:0.5:1,'XTickLabel',{'0','0.5','1'});
set(gca,'Ytick',0:0.5:1,'YTickLabel',{'0','0.5','1'});
title('\beta_x')

subplot(2,ncol,ncol*2);
lo=0;
up=0.31;
plot(pinM(:,7),poutM(:,7),'ro', [lo,up],[lo,up],'-');
axis([lo,up,lo,up]);  axis square;
set(gca,'FontSize',9);
xlabel('Input value');  ylabel('Fit value');
set(gca,'Xtick',0:0.1:0.3,'XTickLabel',{'0','','0.2',''});
set(gca,'Ytick',0:0.1:0.3,'YTickLabel',{'0','','0.2',''});
title('\kappa_d')

subplot(2,ncol,[1:2,ncol+1:ncol+2]);
lo=-4;
up=4;
for i=1:10
    plot(log(Rsigp_simu((i-1)*902+1:i*902,:)),log(Rsigp_fit((i-1)*902+1:i*902,:)),'.','Color',rand(1,3), 'markersize', 2); hold on;
end
plot([lo,up],[lo,up],'-'); hold off
axis([lo,up,lo,up]);  axis square;
set(gca,'FontSize',9);
xlabel('log( simulation value )');  ylabel('log( fit value )');
set(gca,'Xtick',-4:2:6,'XTickLabel',{'-4','-2','0','2','4','6'});
set(gca,'Ytick',-4:2:6,'YTickLabel',{'-4','-2','0','2','4','6'});
title('\sigma_p')

%%
export_fig('d_val2.pdf');
%saveas(gcf,'d_val','pdf')

%% 2015-05-07 Validate bx weight parameter

[PROP] = prop_sim([]);

for parti = 1:3
    for cond = 1:4
        mla = PROP(parti, cond).mla;
        ind0 = mla == 0.5; % actually 0 elements
        ind1 = mla > 0.5;
        ind2 = mla < 0.5;
        [nrow, ncol] = size(mla);
        PROP(parti, cond).mAnswerChoice = nan(nrow, ncol);        
        m0 = round(rand(nrow, ncol));
        PROP(parti, cond).mAnswerChoice(ind0) = m0(ind0);
        PROP(parti, cond).mAnswerChoice(ind1) = 1;
        PROP(parti, cond).mAnswerChoice(ind2) = 2;
    end
end

% changes code in a_optimize_likelihood & a_Get_Model_Perf with these new answers
% validated -- see plot lpxWeight-valid


