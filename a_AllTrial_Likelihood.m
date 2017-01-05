function [LlikRvec,LpriRvec,dpriFvec,dlsigoFvec,dphitFvec,dphisFvec,dbxFvec]=a_AllTrial_Likelihood(DIM,D,AUX)
%% define stuff
npix=770;  ImageSize=20;
prpa=AUX.pvd(3);

LlikRvec=zeros(D.Trials,1); % log likelihood ratio vector
LpriRvec=zeros(D.Trials,1); % log prior ratio vector
dpriFvec=zeros(D.Trials,1); % derivative of prior factor vector
dlsigoFvec=zeros(D.Trials,1); % derivative of observation noise factor vector
dphitFvec=zeros(D.Trials,1); % derivative of temporal forgetting factor vector
dphisFvec=zeros(D.Trials,2); % derivative of spatial accuity rate vector
dbxFvec=zeros(D.Trials,1); % derivative of revealing location bias vector

%%
for trial=1:D.Trials
    
    %Type = D.ImageID(trial,1);
    %X    = D.ImageID(trial,2);    
    %Y    = D.ImageID(trial,3);
    %Num  = D.ImageID(trial,4);
    yimage = DIM(trial).yimage;
        
    nv = D.MaxRevealingTrial(trial);
    xDx = 1+(npix-1)*(D.RevealPosX(trial,1:nv)/ImageSize+0.5);  % cm to pixel: inverting last two lines of a_GET_BALD
    xDy = 1-(npix-1)*(D.RevealPosY(trial,1:nv)/ImageSize-0.5);
    xD = [xDx',xDy'];
    
    [dlml,TrQ,TrphitQ,TrphisQ,lpxzvec]=a_Model(xD, yimage, AUX.pvp, AUX.Moments);

    llr = dlml(1)-log(exp(dlml(2))+exp(dlml(3))); %log likelihood ratio

    %for derivatives
    lwbsh= dlml(2)-log(exp(dlml(2))+exp(dlml(3)));
    lwbsv= dlml(3)-log(exp(dlml(2))+exp(dlml(3)));    
    dlsigoF= 0.5*(TrQ(1) - TrQ(2)*exp(lwbsh) - TrQ(3)*exp(lwbsv));
    dphitF= 0.5*(TrphitQ(1) - TrphitQ(2)*exp(lwbsh) - TrphitQ(3)*exp(lwbsv));
    dphisF= 0.5*(TrphisQ(:,1) - TrphisQ(:,2)*exp(lwbsh) - TrphisQ(:,3)*exp(lwbsv));
    dbxF = lpxzvec(1) - lpxzvec(2)*exp(lwbsh) - lpxzvec(3)*exp(lwbsv);
    %disp(sprintf('v=%f;',dphisF(2))); %debug

    if (D.AnswerChoice(trial)==1)
        LlikRvec(trial) = llr;
        LpriRvec(trial) = log(2*prpa) -log(1-prpa);
        dpriFvec(trial) = 1/prpa/(1-prpa);
        dlsigoFvec(trial) = dlsigoF;
        dphitFvec(trial) = dphitF;
        dphisFvec(trial,:) = dphisF';
        dbxFvec(trial) = dbxF;
    elseif (D.AnswerChoice(trial)==2) % just the negative of the above
        LlikRvec(trial) = -llr;
        LpriRvec(trial) = log(1-prpa) -log(2*prpa);
        dpriFvec(trial) = -1/prpa/(1-prpa);
        dlsigoFvec(trial) = -dlsigoF;
        dphitFvec(trial) = -dphitF;
        dphisFvec(trial,:) = -dphisF';
        dbxFvec(trial) = -dbxF;
    else % didn't make it in time to make a choice
        LlikRvec(trial) = 0;
        LpriRvec(trial) = 0;
        dpriFvec(trial) = 0;
        dlsigoFvec(trial) = 0;
        dphitFvec(trial) = 0;
        dphisFvec(trial,:) = [0,0,0];
        dbxFvec(trial) = 0;
    end
        %disp(sprintf('trial=%d;',trial));
end

end
%% Replaced by [dlml,TrQ,TrphitQ,TrphisQ] = a_Model(xD, yimage, pin, w, phit, phis)

% function [dlml,TrQ,TrphitQ,TrphisQ] = a_AllTrial_Likelihood_sub(xD, yimage, pin, w, phit, phis)
% % w: revealing window width; w=1 recovers the single pixel case
% 
% % get data with window
% xn=size(xD,1);  center=floor(w/2);
% xDw=[];  sDw=[];  tDw=[];
% maxnr=25;  % maximum revealing
% tD=((1:xn)'-xn)/maxnr;  % timed by revealing number
% 
% % % sD model 1
% % cdeg=20/770*sin(10/35)/10*180/pi;  % convert pixel to degree
% % edeg=1/2;  % accuity fall off rate in degree from Najemnik05
% % sD=-1./sum(exp(-sqrt(sq_dist(xD',xD'))*cdeg/edeg),2);
% 
% % sD model: nearest neighbor effect
% % cdeg=20/770*sin(10/35)/10*180/pi;  % convert pixel to degree
% % sD = cdeg*min(sqrt(sq_dist(xD',xD'))+1e100*diag(ones(1,xn)))';
% 
% % sD model: step function noise increas
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
% shhyp.cov = [log(sL1), log(sL2), log(1)];  shhyp.lik = lsigo;
% svhyp.cov = [log(sL2), log(sL1), log(1)];  svhyp.lik = lsigo;
% 
% covfunc = {@covSEard};  %likfunc = {@likGauss};
% 
% %
% [nn, ~] = size(xDw);  %[nt, ~] = size(tDw);
% sn2 = exp(2*lsigo);
% 
% % with sD model 1
% % phitM=diag(phit.^(tDw)).*diag(phis.^(sDw));  % phit=1 is the no forgetting case
% % dphitM=diag((tDw).*phit.^(tDw-1)).*diag(phis.^(sDw));
% % dphisM=diag(phit.^(tDw)).*diag(sDw.*phis.^(sDw-1));
% 
% % with sD model n.n. effect
% phitM=diag(phit.^(tDw))+phis*diag(sDw);
% dphitM=diag((tDw).*phit.^(tDw-1));
% dphisM=diag(sDw);
% 
% K = feval(covfunc{:}, pahyp.cov, xDw);  % evaluate covariance matrix
% L = chol(K+sn2*phitM);  % Cholesky factor of covariance with noise
% alpha = solve_chol(L,yD);
% lmlpa = -((yD)'*alpha/2 + sum(log(diag(L))) + nn*log(2*pi)/2);  % log marg lik
% TrQpa = trace((alpha*alpha')*phitM - solve_chol(L,phitM));  % der wrt sigo
% TrphitQpa = sn2*trace((alpha*alpha')*dphitM - solve_chol(L,dphitM)); %der wrt phit
% TrphisQpa = sn2*trace((alpha*alpha')*dphisM - solve_chol(L,dphisM)); %der wrt phis
% 
% K = feval(covfunc{:}, shhyp.cov, xDw);
% L = chol(K+sn2*phitM);
% alpha = solve_chol(L,yD);
% lmlsh = -((yD)'*alpha/2 + sum(log(diag(L))) + nn*log(2*pi)/2);
% TrQsh = trace((alpha*alpha')*phitM - solve_chol(L,phitM));
% TrphitQsh = sn2*trace((alpha*alpha')*dphitM - solve_chol(L,dphitM));
% TrphisQsh = sn2*trace((alpha*alpha')*dphisM - solve_chol(L,dphisM)); %der wrt phis
% 
% K = feval(covfunc{:}, svhyp.cov, xDw);
% L = chol(K+sn2*phitM);
% alpha = solve_chol(L,yD);
% lmlsv = -((yD)'*alpha/2 + sum(log(diag(L))) + nn*log(2*pi)/2);
% TrQsv = trace((alpha*alpha')*phitM - solve_chol(L,phitM));
% TrphitQsv = sn2*trace((alpha*alpha')*dphitM - solve_chol(L,dphitM));
% TrphisQsv = sn2*trace((alpha*alpha')*dphisM - solve_chol(L,dphisM)); %der wrt phis
% 
% lml = max([lmlpa,lmlsh,lmlsv]);
% dlmlpa = lmlpa-lml;
% dlmlsh = lmlsh-lml;
% dlmlsv = lmlsv-lml;
% 
% dlml=[dlmlpa,dlmlsh,dlmlsv];
% TrQ=[TrQpa,TrQsh,TrQsv];
% TrphitQ=[TrphitQpa,TrphitQsh,TrphitQsv];
% TrphisQ=[TrphisQpa,TrphisQsh,TrphisQsv];
% 
% end

%% OLD CODE using gp

% marginal likelihood for p(model|D, theta)
% lmlpa = -gp(pahyp, @infExact, [], covfunc, likfunc, xDw, yD);
% lmlsh = -gp(shhyp, @infExact, [], covfunc, likfunc, xDw, yD);
% lmlsv = -gp(svhyp, @infExact, [], covfunc, likfunc, xDw, yD);

%% Carl's code for finding log marginal likelihood and its derivative

% K = feval(covfunc{:}, svhyp.cov, xDw);
% L = chol(K/sn2+eye(nn));
% alpha = solve_chol(L,yD)/sn2;
% lmlsv = -((yD)'*alpha/2 + sum(log(diag(L))) + nn*log(2*pi*sn2)/2);
% TrQsv = trace(alpha*alpha' - solve_chol(L,eye(nn))/sn2);

% moving the sn2 to its proper place as is in Carl's book
% K = feval(covfunc{:}, shhyp.cov, xDw);
% L = chol(K+sn2*eye(nn));
% alpha = solve_chol(L,yD);
% lmlsh = -((yD)'*alpha/2 + sum(log(diag(L))) + nn*log(2*pi)/2);
% TrQsh = trace(alpha*alpha' - solve_chol(L,eye(nn)));
    
