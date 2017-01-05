% replaced by a_Model
function [mla, matH_BALD, xsH_BALD] = a_BALDscore(xD, yimage, pvp, pvd, n, yrin)
%% PART 1: setup parameters
w=pvp(1);
lsigo=pvp(2);
phit=pvp(3);
phis=pvp(4);
phis0=pvp(5);
lphismax=pvp(6);

if(nargin>3)
    bd=pvd(1);
    lapse=pvd(2);
    prpa=pvd(3);
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
xx = linspace(1,npix,n);  %1:npix
yy = linspace(1,npix,n);  %1:npix
[x1,x2]=meshgrid(xx,yy);
z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % in pixel x,y coordinate
nz=size(z,1);
matlml=zeros(xn,3);  matmu=zeros(nz+1,3,xn-1);  matvar=matmu;
matH_BALD=zeros(nz,xn-1);  xsH_BALD=zeros(1,xn-1);
mla=NaN(1,xn);


%% PART 4: initialize model
pahyp.cov = [log(npix/20), log(npix/20), log(1)];  %pahyp.lik = lsigo;
shhyp.cov = [log(npix/6), log(npix/30), log(1)];  %shhyp.lik = lsigo;
svhyp.cov = [log(npix/30), log(npix/6), log(1)];  %svhyp.lik = lsigo;
covfunc = {@covSEard};  %likfunc = {@likGauss};


%% PART 5: GP calculation

sn2 = exp(2*lsigo);
cdeg = 20/770*sin(10/42)/10*180/pi; % convert pixel to degree

hyper(1,:)=pahyp.cov;
hyper(2,:)=shhyp.cov;
hyper(3,:)=svhyp.cov;

for iht=1:3
    K = feval(covfunc{:}, hyper(iht,:), xDw);  % evaluate covariance matrix    
    kss = ones(nz+1,1); %feval(covfunc{:}, hyper(iht,:), z, 'diag');  %kss=nz-by-nz
    Ks  = feval(covfunc{:}, hyper(iht,:), xDw, z);  %Ks=xDn-by-nz
    for nr=1:xn-1
        subi=nr*w*w;
        
        subxDw = xDw(1:subi,:);
        rp = [subxDw; zeros(1,2)];
        fp = rp; %rp(end)=fp(end)-->fixation on z
        nTd = nr+1;
        rt = reshape(repmat((0:nTd-1),w*w,1),nTd*w*w,1); % revealing time
        ft1 = rt; % fixation start time
        ft2 = ft1+1; % fixation end time

        drf = cdeg*sqrt(sq_dist(rp',fp')); % distance between revealing and fixation positions
        t1rf = repmat(ft1',nTd*w*w,1)-repmat(rt,1,nTd*w*w); % time fixation start relative to revealing
        t1rf(t1rf<0)=0;
        t2rf = repmat(ft2',nTd*w*w,1)-repmat(rt,1,nTd*w*w); % time fixation end relative to revealing
        t2rf(t2rf<0)=0;

        stM = (exp(phit*t2rf)-exp(phit*t1rf)).*exp(-phis*drf); % spatiotemporal modulation matrix
        stV = sum(stM,2); % spatiotemporal modulation vector
        phiM = diag( phit*exp(phit*nTd)./stV ); % final observation noise modulation
        stmz = phiM(end,end);
        phiM(:,end)=[];  phiM(end,:)=[];
        

        % posterior at t: previous opening > t > coming opening 
        subK = K(1:subi,1:subi);  %checked numerically, submatrix of L = chol(submatrix of K)
        subL = chol(subK+sn2*phiM);  % Cholesky factor of covariance with noise
        suby = yD(1:subi,:);
        phiV = diag(phiM);
        suby = suby+yrin(1:subi).*sqrt(sn2*phiV);
        
        suba = solve_chol(subL,suby); %has to recompute this
        matlml(nr,iht) = -((suby)'*suba/2 + sum(log(diag(subL))) + subi*log(2*pi)/2);
       
        subKs = [Ks(1:subi,:), feval(covfunc{:}, hyper(iht,:), xDw(1:subi,:), xD(nr+1,:))]; % add xs to subKs
        matmu(:,iht,nr) = subKs'*suba;
        subV = subL'\subKs;
        fs2 = kss - sum(subV.*subV,1)' + sn2*stmz;
        matvar(:,iht,nr) = max(fs2,0);

    end
end


%% PART 6b: posterior and BALD for BALD score property
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
    %lws = log(exp(lwsh)+exp(lwsv));

    % Mate likes soft lik NOT passed to movement model, only to decision model    
    %lpr = lwpa -lws;
    %softlik = 1/(1+exp(-bd*lpr));
    %lapse_softlik = (1-lapse)*softlik + 0.5*lapse;
    %llslpa = log(lapse_softlik);
    %llsls = log(1-lapse_softlik);
    %llslsh = llsls + lwsh - lws;
    %llslsv = llsls + lwsv - lws;
    %lw  = [llslpa,llslsh,llslsv];

    lw = [lwpa, lwsh, lwsv];
    mu=reshape(matmu(:,:,nr),(nz+1)*3,1);
    var=reshape(matvar(:,:,nr),(nz+1)*3,1);

    H_BALD = a_GetBALDGM_mex(lw, mu, var);
    H_BALD = max(0, H_BALD);
    matH_BALD(:,nr)=H_BALD(1:end-1);
    xsH_BALD(1,nr)=H_BALD(end);
    mla(nr)=exp(lw(1))/sum(exp(lw));
end


end
