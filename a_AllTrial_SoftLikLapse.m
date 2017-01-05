function [val, dval]=a_AllTrial_SoftLikLapse(x,w,D,DIM,AUX)
% For optimization: [val, dval]=a_AllTrial_SoftLikLapse(x,w,D,DIM,AUX)
% UNFINISHED: function [val, dval]=a_AllTrial_SoftLikLapse(x,w,D,DIM,yrM,Moments, exPerfM, exPerfMstd, cond)

%% define stuff
lsigo=x(1);  sigo=exp(lsigo);
phit=x(2);
lphis=x(3);  phis=exp(lphis);
phis0=x(4);
lbx=x(5);  bx=exp(lbx);
pvp=[w,sigo,phit,phis,phis0,bx];  % parameters vector for perception

lbd=x(6);  bd=exp(lbd);
lapse=x(7);
prpa=x(8);
pvd=[bd,lapse,prpa]; % parameters vector for decision

AUX.pvp=pvp;
AUX.pvd=pvd;

%% sum of trial-by-trial likelihood
[LlikRvec,LpriRvec,dpriFvec,dlsigoFvec,dphitFvec,dphisFvec,dbxFvec] = a_AllTrial_Likelihood(DIM,D,AUX);
lpr = LlikRvec+LpriRvec;  % log posterior ratio

% make sure value not <Inf and >-Inf
lpr = max(-1000,lpr); lpr=min(1000,lpr);
dlsigoFvec = max(-1000,dlsigoFvec); dlsigoFvec=min(1000,dlsigoFvec);
dphitFvec = max(-1000,dphitFvec); dphitFvec=min(1000,dphitFvec);
dphisFvec = max(-1000,dphisFvec); dphisFvec=min(1000,dphisFvec);
dbxFvec = max(-1000,dbxFvec); dbxFvec=min(1000,dbxFvec);

softlik = 1./(1+exp(-bd*lpr));  % softmax of lpr with bd
lapse_softlik = (1-lapse)*softlik + 0.5*lapse;  % add lapse rate to softlik

% precompute a common factor:
temp=(1-lapse)./lapse_softlik.*softlik./(1+exp(bd*lpr));

%gradients
dval_lsigo = bd*temp.*dlsigoFvec; % gradient wrt lsigo
dval_phit = bd*temp.*dphitFvec; % gradient wrt phit
dval_lphis = bd*temp.*dphisFvec(:,1); % gradient wrt phis
dval_phis0 = bd*temp.*dphisFvec(:,2); % gradient wrt phis0
dval_lbx = bd*temp.*dbxFvec; % gradient wrt bx
dval_lbd = temp.*lpr; % gradient wrt bd
dval_lapse = (0.5-softlik)./lapse_softlik; % gradient wrt lapse
dval_prpa = temp.*(bd*dpriFvec); % gradient wrt ppa: same as for wrt to b1; just replace lpr --> b1/ppa/(1-ppa)

% output1
val = -sum(log(lapse_softlik));  % negative log posterior of all trials
% output2
dval(1)=-sigo*sum(dval_lsigo);
dval(2)=-sum(dval_phit);
dval(3)=-phis*sum(dval_lphis);
dval(4)=-sum(dval_phis0);
dval(5)=-bx*sum(dval_lbx);
dval(6)=-bd*sum(dval_lbd);
dval(7)=-sum(dval_lapse);
dval(8)=-0*sum(dval_prpa);

%% UNFINISHED
% % for sum of condition-by-condition (condition averaged) likelihood
% if(nargin>npbase) %cond is a n-by-2 matrix, stores the start and end index of each condition
%     [nc,nrn]=size(exPerfM);
%     PerfM=zeros(nc,nrn);  dPerfM=zeros(nc*nrn,8);
%     for c=1:size(cond,1)
%         c1=cond(c,1); c2=cond(c,2);
%         MaxReveal = D.MaxRevealingTrial(c1:c2);
%         Prob = lapse_softlik(c1:c2);
%         der=zeros(c2-c1+1,8);
%         der(:,1)=dval_lsigo(c1:c2);
%         der(:,2)=dval_phit(c1:c2);
%         der(:,3)=dval_phis(c1:c2);
%         der(:,4)=0;
%         der(:,5)=0;
%         der(:,6)=dval_bd(c1:c2);
%         der(:,7)=dval_lapse(c1:c2);
%         der(:,8)=0;
%     
%         RevN=MaxReveal;
%         P = Prob;
%         uni_revn = unique(RevN);
%         rn =length(uni_revn);
%         perf=zeros(rn,1);
%         dperf=zeros(rn,8); %nrn-by-8
%         for i=1:rn
%            q=P(RevN==uni_revn(i));
%            dq=der(RevN==uni_revn(i),:);
%            perf(i) = mean(q);
%            dperf(i,:) = mean(dq,1);
%         end
%         
%         PerfM(c,:)=perf'; %total 3-by-5: #cond-by-#rn
%         dPerfM( (c-1)*nrn+1:c*nrn, :)=dperf;        
%     end
%     
%     PerfM =reshape(PerfM',nc*nrn,1);
%     exPerfM =reshape(exPerfM',nc*nrn,1);
%     exPerfMstd =reshape(exPerfMstd,nc*nrn,1);
%     
%     val=sum( (PerfM-exPerfM).^2 ./ exPerfMstd.^2 /2);     
%     dval=sum( repmat( (PerfM-exPerfM) ./ exPerfMstd.^2 , 1,8).*dPerfM, 1);    
% end


end

    
