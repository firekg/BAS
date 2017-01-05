function [varargout] = a_Model(xD, yimage, pvp, Moments, pvd, yrin, z)
% For optimization  : [dlml,TrQ,TrphitQ,TrphisQ] = a_Model(xD, yimage, pvp, Moments)
% For performance +1: [lapse_softlik,lw] = a_Model(xD, yimage, pvp, Moments, pvd)
% For validation  +2: [lapse_softlik] = a_Model(xD, yimage, pvp, Moments, pvd, yrin)
% For BAS score   +3: [mla, matH_BAS, xsH_BAS] = a_Model(xD, yimage, pvp, Moments, pvd, yrin, z)

%% Deficiencies:
% sq_dist in PerceptionNoise() and subKs=[...,feval()] is being called more than needed?
% PerceptionNoise is not coded with actual fixation position and times yet

%% PART 1: setup parameters

npmin=4; %minimum number of parameters

w=pvp(1);
sigo=pvp(2);  sigo2=sigo^2; %perception noise variance level
% phit=pvp(3);
% phis=pvp(4);
% phis0=pvp(5);
bx=pvp(6);

if(nargin>npmin)
    bd=pvd(1);
    lapse=pvd(2);
    prpa=pvd(3); %should belong to pvp...
end


%% PART 2: get data 
% position data with window
xn=size(xD,1);  center=floor(w/2);
xDw1=[];  xDw2=[];
for i=0:w-1
for j=0:w-1
    xDw1=[xDw1, xD(:,1)+(i-center)*ones(xn,1)]; %xn-by-w^2
    xDw2=[xDw2, xD(:,2)+(j-center)*ones(xn,1)]; %xn-by-w^2
end
end

xDw = [reshape(xDw1',xn*w*w,1), reshape(xDw2',xn*w*w,1)];
xDw=round(xDw);
npix=770; xDw(xDw<1)=1; xDw(xDw>npix)=npix;

% color data with window
xDn=size(xDw,1);  yD=zeros(xDn,1);
for i=1:xDn
    yD(i) = yimage(xDw(i,2),xDw(i,1));  % x is column of image, y is row of image;
end


%% PART 3: initialize things & set up probe points z in pixel 
lmlv=zeros(1,3);
lpxzvec=zeros(1,3);

if(nargin==npmin)
    TrQ=lmlv; TrphitQ=lmlv; TrphisQ=zeros(2,3);
end

if(nargin==npmin+3)
    nz=size(z,1);
    matlml=zeros(xn,3);  matmu=zeros(nz+1,3,xn-1);  matvar=matmu;
    matH_BAS=zeros(nz,xn-1);  xsH_BAS=zeros(1,xn-1);
    mla=NaN(1,xn);
end

%% PART 4: building perception noise model

if (nargin<npmin+3)
    % get data: revealing and fixation positions and times
    rp = xDw; % revealing positions
    fp = [rp;rp]+1e-3;  % retrace revealed position; add offset to help phis derivative
    % fp = [rp;flipud(rp)];  % retrace backward revealed position
    % fp = [rp; repmat(mean(rp,1),xDn,1)];  % stay on most crowded place
    nTd = size(fp,1); % decision time in number of points
    rt = reshape(repmat((0:xn-1),w*w,1),xn*w*w,1); % revealing time
    ft1 = reshape(repmat((0:nTd-1),w*w,1),nTd*w*w,1); % fixation start time
    ft2 = ft1+1; % fixation end time

    if(nargin==npmin)
        [phiM,dphitM,dphisM]=PerceptionNoise(pvp,rp,fp,rt,ft1,ft2,nTd,2);
    else
        [phiM]=PerceptionNoise(pvp,rp,fp,rt,ft1,ft2,nTd,1);
    end
end


%% PART 5: initialize model
pahyp.cov = [log(npix/20), log(npix/20), log(1)];  %pahyp.lik = lsigo; %now a complicated observation noise
shhyp.cov = [log(npix/6), log(npix/30), log(1)];  %shhyp.lik = lsigo;
svhyp.cov = [log(npix/30), log(npix/6), log(1)];  %svhyp.lik = lsigo;
covfunc = {@covSEard};  %likfunc = {@likGauss};


%% PART 6: GP calculation
hyper(1,:)=pahyp.cov;
hyper(2,:)=shhyp.cov;
hyper(3,:)=svhyp.cov;

snI2 = 0.01^2; %image noise variance; hard-wired

% add perception nosie to presented color
if (nargin==npmin+2)
    phiV = diag(phiM);
    yD = yD+yrin(1:xDn).*sqrt(sigo2*phiV); %yD already has snI in it
end
% a bad trick: pre-computed yr makes optimization deterministic

if (nargin<npmin+3)
for iht=1:3 %image hyperparameter type
    K = feval(covfunc{:}, hyper(iht,:), xDw);  % evaluate covariance matrix K
    L = chol(K + snI2*eye(xDn) + sigo2*phiM);  % Cholesky factor of covariance with image and perception noise
    alpha = solve_chol(L,yD);
    lmlv(iht) = -((yD)'*alpha/2 + sum(log(diag(L))) + xDn*log(2*pi)/2);  % log marginal lik vector
    if (nargin==npmin)
        TrQ(iht) = 2*sigo*trace((alpha*alpha')*phiM - solve_chol(L,phiM));  % deri wrt sigo
        TrphitQ(iht) = sigo2*trace( (alpha*alpha')*dphitM-solve_chol(L,dphitM) ); %deri wrt phit
        TrphisQ(1,iht) = sigo2*trace( (alpha*alpha')*dphisM(:,:,1)-solve_chol(L,dphisM(:,:,1)) ); %deri wrt phis
        TrphisQ(2,iht) = sigo2*trace( (alpha*alpha')*dphisM(:,:,2)-solve_chol(L,dphisM(:,:,2)) ); %deri wrt phis0
    end
end
end

%% PART 6.5: GP calculation for BAS

if (nargin==npmin+3)
for iht=1:3
    K = feval(covfunc{:}, hyper(iht,:), xDw);  % evaluate covariance matrix    
    kss = ones(nz+1,1);
    Ks  = feval(covfunc{:}, hyper(iht,:), xDw, z);  %Ks=xDn-by-nz
    for nr=1:xn-1
        subi=nr*w*w;
        
        % noise model
        subxDw = xDw(1:subi,:);
        rp = [subxDw; 1e3*zeros(1,2)]; %assume noise on z ~ that on a far away revealing
        fp = rp+1e-3; %rp(end)=fp(end)-->fixation on z
        nTd = nr+1;
        rt = reshape(repmat((0:nTd-1),w*w,1),nTd*w*w,1); % revealing time
        ft1 = rt; % fixation start time
        ft2 = ft1+1; % fixation end time
        [phiM]=PerceptionNoise(pvp,rp,fp,rt,ft1,ft2,nTd,1);
        stmz = phiM(end,end);
        phiM(:,end)=[];  phiM(end,:)=[];
        
        % posterior at t: previous opening > t > coming opening 
        subK = K(1:subi,1:subi);  %checked numerically, submatrix of L = chol(submatrix of K)
        subL = chol(subK + snI2*eye(subi) + sigo2*phiM);  % Cholesky factor of covariance with noise
        suby = yD(1:subi,:);
        
        % simu perceived values
        % phiV = diag(phiM);
        % suby = suby+yrin(1:subi).*sqrt(sigo2*phiV);
        
        suba = solve_chol(subL,suby); %has to recompute this
        matlml(nr,iht) = -((suby)'*suba/2 + sum(log(diag(subL))) + subi*log(2*pi)/2);
               
        subKs = [Ks(1:subi,:), feval(covfunc{:}, hyper(iht,:), xDw(1:subi,:), xD(nr+1,:))]; % add xs to subKs
        matmu(:,iht,nr) = subKs'*suba;
        subV = subL'\subKs;
        fs2 = kss - sum(subV.*subV,1)' + sigo2*stmz; %assuming noise on xs = that on z
        matvar(:,iht,nr) = max(fs2,0); %predictive variance on z and xs

    end
end
end

%% PART 7: P(z|x) effect
% this is only good for window=1

%adding log[P(z|x)] to lmlv
if(nargin<npmin+3)
    for i=2:xDn
        xy = xDw(i,:);
        lpxzvec(1) = lpxzvec(1)+ 0.5*( -log(2*pi) -log(Moments.dcovpa(i)) -(xy-Moments.mVpa(i,:))*(Moments.icovMpa(:,:,i)*(xy-Moments.mVpa(i,:))'));
        lpxzvec(2) = lpxzvec(2)+ 0.5*( -log(2*pi) -log(Moments.dcovsh(i)) -(xy-Moments.mVsh(i,:))*(Moments.icovMsh(:,:,i)*(xy-Moments.mVsh(i,:))'));
        lpxzvec(3) = lpxzvec(3)+ 0.5*( -log(2*pi) -log(Moments.dcovsv(i)) -(xy-Moments.mVsv(i,:))*(Moments.icovMsv(:,:,i)*(xy-Moments.mVsv(i,:))'));
    end
    lmlv = lmlv + bx*lpxzvec;
end

%% PART 8: posterior
if(nargin<npmin+3)
% difference log marginal likelihood (for numerical precision)
lml = max(lmlv);
dlmlpa = lmlv(1)-lml;
dlmlsh = lmlv(2)-lml;
dlmlsv = lmlv(3)-lml;
dlml=[dlmlpa,dlmlsh,dlmlsv];

% add softmax and lapse rate
if(nargin>npmin)
    psh=(1-prpa)/2;  psv=(1-prpa)/2;
    lnor = log(prpa*exp(dlmlpa)+psh*exp(dlmlsh)+psv*exp(dlmlsv));

    lwpa = log(prpa)+dlmlpa-lnor;
    lwsh = log(psh)+dlmlsh-lnor;
    lwsv = log(psv)+dlmlsv-lnor;
    lws = log(exp(lwsh)+exp(lwsv)); %log weight stripy

    lw = [lwpa, lwsh, lwsv];
    lpr = lwpa -lws;  % for choosing Type 1 (P)
    softlik = 1/(1+exp(-bd*lpr));  % softmax of lpr with bd
    lapse_softlik = (1-lapse)*softlik + 0.5*lapse;  % add lapse rate to softlik    
end
end

%% PART 8.5: BAS output prep
if(nargin==npmin+3)
for nr=1:xn-1
    lml=max(matlml(nr,:));
    dlmlpa = matlml(nr,1)-lml;
    dlmlsh = matlml(nr,2)-lml;
    dlmlsv = matlml(nr,3)-lml;

    psh=(1-prpa)/2;  psv=(1-prpa)/2;
    lnor = log(prpa*exp(dlmlpa)+psh*exp(dlmlsh)+psv*exp(dlmlsv));

    lwpa = log(prpa)+dlmlpa-lnor;
    lwsh = log(psh)+dlmlsh-lnor;
    lwsv = log(psv)+dlmlsv-lnor;

    lw = [lwpa, lwsh, lwsv];
    mu=reshape(matmu(:,:,nr),(nz+1)*3,1);
    var=reshape(matvar(:,:,nr),(nz+1)*3,1);

    H_BAS = a_GetBALDGM_mex(lw, mu, var);
        
    H_BAS = max(0, H_BAS);
    matH_BAS(:,nr)=H_BAS(1:end-1);
    xsH_BAS(1,nr)=H_BAS(end);
    mla(nr)=exp(lw(1))/sum(exp(lw));
end
end


%% PART 9: Final output
if(nargin==npmin)
    varargout={dlml,TrQ,TrphitQ,TrphisQ,lpxzvec};
elseif(nargin==npmin+1)
    varargout={lapse_softlik,lw};
elseif(nargin==npmin+2)
    sigma_p=sqrt(sigo2./stV);
    varargout={lapse_softlik,dlml,sigma_p};
elseif(nargin==npmin+3)
    varargout={mla, matH_BAS, xsH_BAS};
end

end

%% ########################################################################
%% Perception noise function
function [varargout]=PerceptionNoise(pvp,rp,fp,rt,ft1,ft2,nTd,mode)
% mode 1: just gives phiM
% mode 2: also gives derivatives

w=pvp(1);
phit=pvp(3);
phis=pvp(4);
phis0=pvp(5);
nr = size(rp,1);

    %space (better to convert to degree, see notes 2014-01-09 on a way to do that)
    pix2cm = 20/770; % convert pixel to cm
    drf = pix2cm*sqrt(sq_dist(rp',fp')); % distance between revealing and fixation positions in cm

    %time
    rtM = repmat(rt,1,nTd*w*w);  %the size is not nTd in general!
    ft1M = repmat(ft1',nr*w*w,1);
    ft2M = repmat(ft2',nr*w*w,1);
    t1rf = ft1M.*(ft1M>=rtM); % time fixation start thresholded by revealing time
    t2rf = ft2M.*(ft2M>=rtM);

    %perception noise modulation
    %sV = 1-1./(1+exp(-phis*(drf-phis0))); %spatial modulation vector -- sigmoid
    sVe=log(2)*(drf/phis0).^(phis);  sV=exp(-sVe); %spatial modulation vector
    tV = phit*exp(-nTd/phit)*(exp(t2rf/phit)-exp(t1rf/phit)); %temporal modulation vector
    stV = sum(sV.*tV,2); % spatiotemporal modulation vector
    phiM = diag(1./stV); % final observation noise modulation without overall noise level
    % disp(sprintf('max v=%f; min v=%f.',max(1./stV),min(1./stV))); %debug
    
    if(mode==2)
        % derivative wrt phit
        nV = sum( -sV.*( exp(-nTd/phit)*( (1+nTd/phit-t2rf/phit).*exp(t2rf/phit) - (1+nTd/phit-t1rf/phit).*exp(t1rf/phit) )),2); %a numerator vector
        dphitM = diag( nV./(stV.^2) );
        % derivative wrt phis
        %nV = sum( ((drf-phis0).*exp(-phis*(drf-phis0)).*((1-sV).^2) ).*tV ,2); %for sigmoid
        nV = sum( log(drf/phis0).*sVe.*sV.*tV ,2);
        dphisM(:,:,1)= diag( nV./(stV.^2) );
        % derivative wrt phis0
        %nV = sum( -(phis*exp(-phis*(drf-phis0)).*((1-sV).^2)).*tV ,2); %for sigmoid
        nV = sum( -phis/phis0.*sVe.*sV.*tV ,2);
        dphisM(:,:,2)= diag( nV./(stV.^2) );
    end
   
    if(mode==1)
        varargout={phiM};
    else
        varargout={phiM,dphitM,dphisM};
    end
end
