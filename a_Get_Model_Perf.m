function [perf, sd, Prob, LPR]=a_Get_Model_Perf(D,Gmap,AUX,pv)

% npix=770;  ImageSize=20;
% w=3; pin=[ npix/20, npix/20, npix/6, npix/30, -0.2867];
% b1=0.2687; lapse=0.1015;
% ppa=0.2248;

npix=770;
ImageSize=20;

w=pv(9);
lsigo=pv(1);  sigo = exp(lsigo);
phit=pv(2);
phis=pv(3);
phis0=pv(4);
lbx=pv(5);  bx=exp(lbx);
shift = pv(12);
pvp=[w,sigo,phit,phis,phis0,bx,shift];  % parameters vector for perception

lbd=pv(6);  bd = exp(lbd);
lapse=pv(7);
prpa=pv(8);
prxpa = pv(10);
prCommon = pv(11);
pvd=[bd,lapse,prpa, prxpa, prCommon];  % parameters vector for decision

% Correct = zeros(D.Trials,1);  %also a probability
Prob = zeros(D.Trials,1);
LPR = zeros(D.Trials,1);

% ProbAB = zeros(D.Trials,2);
% lmlM = zeros(D.Trials,3);

for trial=1:D.Trials    
    %yimage = DIM(trial).yimage;

    nv = D.MaxRevealingTrial(trial);
    xDx = 1+(npix-1)*(D.RevealPosX(trial,1:nv)/ImageSize+0.5);  % cm to pixel
    xDy = 1-(npix-1)*(D.RevealPosY(trial,1:nv)/ImageSize-0.5);
    xD = [xDx',xDy'];
    yrin = randn(nv,1);  % consider fixing random seed?
    yD = D.RevealZ(trial, 1:nv)';
    AUX.Moments.lpx = Gmap(trial,nv,:);    
    
    %[Answer,lapse_softlik,lmlpa,lmlsh,lmlsv]=a_Get_Model_Answer(xD, yimage, pin, w, prpa, b1, lapse, phit, phis);
    %Correct(trial) = Answer==D.AnswerReal(trial);
    %ProbAB(trial,1)=lapse_softlik;
    %ProbAB(trial,2)=1-lapse_softlik;
    %lmlM(trial,1)=lmlpa;
    %lmlM(trial,2)=lmlsh;
    %lmlM(trial,3)=lmlsv;

    [lapse_softlik, lpr, ~] = a_Model(xD, yD, pvp, AUX.Moments, pvd, yrin);
    
%     % for getting correct
%     if(D.AnswerReal(trial)==1)
%         Correct(trial) = lapse_softlik;        
%     else
%         Correct(trial) = 1 - lapse_softlik;        
%     end
    
    % for getting Prob
    if(lapse_softlik==0.5)
        lapse_softlik = rand();
    end    
    if(D.AnswerChoice(trial)==1)
        Prob(trial) = (lapse_softlik > 0.5);
        LPR(trial) = lpr;
    else
        Prob(trial) = 1 - (lapse_softlik > 0.5);
        LPR(trial) = -lpr;
    end

%     % for validation
%     if(D.mAnswerChoice(trial,nv)==1)
%         Prob(trial) = (lapse_softlik > 0.5);        
%     else
%         Prob(trial) = 1 - (lapse_softlik > 0.5);        
%     end
    
    %fprintf('sigo=%f, bx=%f, softlik =%f;\n',pvp(2),pvp(6),lapse_softlik);    
    %disp(sprintf('Prob=%f;',lapse_softlik));    
    %disp(sprintf('trial=%d;',trial));
end

RevN = D.MaxRevealingTrial;
uni_revn = unique(RevN);
rn = length(uni_revn);
perf = zeros(rn,1);
sd = zeros(rn,1);
% P = Correct;
% for i=1:rn
%    q = P(RevN==uni_revn(i)); 
%    perf(i) = mean(q);
%    total = length(q);
%    sd(i) = sqrt(var(q)/total);
% end

% x=uni_revn;

end
%% Replaced by [Answer,lapse_softlik]=a_Model(xD, yimage, pin, w, phit, phis, b1, lapse, prpa);

% function [Answer,lapse_softlik,lmlpa,lmlsh,lmlsv] = a_Get_Model_Answer(xD, yimage, pin, w, prpa, b1, lapse, phit, phis)
% 
% % get data with window
% xn=size(xD,1);  center=floor(w/2);
% xDw=[];  sDw=[];  tDw=[];
% maxnr=25;  % maximum revealing
% tD=((1:xn)'-xn)/maxnr;  % timed by revealing number
% 
% cdeg=20/770*sin(10/35)/10*180/pi;  % convert pixel to degree
% edeg=4;
% sD = cdeg*min(sqrt(sq_dist(xD',xD'))+1e100*diag(ones(1,xn)))';
% sD(sD<edeg)=1; sD(sD>=edeg)=0;
% 
% for i=0:w-1
% for j=0:w-1
%     xDw=[xDw; xD+repmat([i,j]-center,xn,1)];
%     tDw=[tDw; tD];
%     sDw=[sDw; sD];
% end
% end
% 
% xDw=round(xDw);
% npix=770; xDw(xDw<1)=1; xDw(xDw>npix)=npix;
% 
% xDn=size(xDw,1);  yD=zeros(xDn,1);
% for i=1:xDn
%     yD(i) = yimage(xDw(i,2),xDw(i,1));  % x is column of image, y is row of image;
% end
% 
% % initialize guess models
% paL1=pin(1);
% paL2=pin(2);
% sL1=pin(3);
% sL2=pin(4);
% lsigo=pin(5);
% 
% pahyp.cov = [log(paL1), log(paL2), log(1)];  pahyp.lik = lsigo;
% shhyp.cov = [log(sL1), log(sL2), log(1)];  shhyp.lik = lsigo;  psh=(1-prpa)/2;
% svhyp.cov = [log(sL2), log(sL1), log(1)];  svhyp.lik = lsigo;  psv=(1-prpa)/2;
% 
% covfunc = {@covSEard};  %likfunc = {@likGauss};
% 
% %
% [nn, ~] = size(xDw);  %[nt, ~] = size(tDw);
% 
% sn2 = exp(2*lsigo);
% phitM=diag(phit.^(tDw))+phis*diag(sDw);
% 
% K = feval(covfunc{:}, pahyp.cov, xDw);  % evaluate covariance matrix
% L = chol(K+sn2*phitM);  % Cholesky factor of covariance with noise
% alpha = solve_chol(L,yD);
% lmlpa = -((yD)'*alpha/2 + sum(log(diag(L))) + nn*log(2*pi)/2);  % log marg lik
% 
% K = feval(covfunc{:}, shhyp.cov, xDw);
% L = chol(K+sn2*phitM);
% alpha = solve_chol(L,yD);
% lmlsh = -((yD)'*alpha/2 + sum(log(diag(L))) + nn*log(2*pi)/2);
% 
% K = feval(covfunc{:}, svhyp.cov, xDw);
% L = chol(K+sn2*phitM);
% alpha = solve_chol(L,yD);
% lmlsv = -((yD)'*alpha/2 + sum(log(diag(L))) + nn*log(2*pi)/2);
% 
% lml = max([lmlpa,lmlsh,lmlsv]);
% dlmlpa = lmlpa-lml;
% dlmlsh = lmlsh-lml;
% dlmlsv = lmlsv-lml;
% 
% lnor = log(prpa*exp(dlmlpa)+psh*exp(dlmlsh)+psv*exp(dlmlsv));
% 
% lwpa = log(prpa)+dlmlpa-lnor;
% lwsh = log(psh)+dlmlsh-lnor;
% lwsv = log(psv)+dlmlsv-lnor;
% 
% lpr = lwpa -log(exp(lwsh)+exp(lwsv));
% softlik = 1/(1+exp(-b1*lpr));  % softmax of lpr with b1
% lapse_softlik = (1-lapse)*softlik + 0.5*lapse;  % add lapse rate to softlik
% 
% if ( rand(1)<=lapse_softlik )
%     Answer=1;
% else
%     Answer=2;
% end
% 
% % % the optimal decision making answer
% % if (isnan(mla) || isnan(mlb))
% %     Answer=-1;
% % elseif (mla==mlb)
% %     Answer = ceil(rand(1)*2);
% % else
% %     [~,Answer]=max([mla,mlb]);
% % end
% 
% 
% end

%% OLD CODE using gp

% % get data with window
% xn=size(xD,1);  xDw=[];  center=floor(w/2);
% for i=0:w-1
% for j=0:w-1
%     xDw=[xDw; xD+repmat([i,j]-center,xn,1)];
% end
% end
% 
% xDw=round(xDw);
% npix=770; xDw(xDw<1)=1; xDw(xDw>npix)=npix;
% 
% xDn=size(xDw,1);  yD=zeros(xDn,1);
% for i=1:xDn
%     yD(i) = yimage(xDw(i,2),xDw(i,1)); % x is column in image, y is row in image;
% end
% 
% % initialize guess models
% MaL1=pin(1);
% MaL2=pin(2);
% MbL1=pin(3);
% MbL2=pin(4);
% Mabn=pin(5);
% 
% pahyp.cov = [log(MaL1), log(MaL2), log(1)];  pahyp.lik = Mabn;  %ppa=0.5;
% shhyp.cov = [log(MbL1), log(MbL2), log(1)];  shhyp.lik = Mabn;  psh=(1-prpa)/2;
% svhyp.cov = [log(MbL2), log(MbL1), log(1)];  svhyp.lik = Mabn;  psv=(1-prpa)/2;
% 
% covfunc = {@covSEard};  likfunc = {@likGauss};
% 
% % marginal likelihood for p(model|D, theta)
% lmlpa = -gp(pahyp, @infExact, [], covfunc, likfunc, xDw, yD);
% lmlsh = -gp(shhyp, @infExact, [], covfunc, likfunc, xDw, yD);
% lmlsv = -gp(svhyp, @infExact, [], covfunc, likfunc, xDw, yD);
% 
% % % final marginal likelihoods
% % nor = ppa*exp(-nlmlpa) + psh*exp(-nlmlsh) + psv*exp(-nlmlsv);
% % mla = ppa*exp(-nlmlpa)/nor;
% % mlb = (psh*exp(-nlmlsh) + psv*exp(-nlmlsv))/nor;    
    
