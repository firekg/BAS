function [mla, mlb, H_BALD, z] = Get_BALDmap_window(xD, yimage, pin, w)
% w: revealing window width; w=1 recovers the single pixel case
% %% for testing
% y=imread('image/Type1X20Y20Num6.bmp');  y=cast(y,'double');
% yimage=(y(:,:,1))/255*10-5;  %pixel value max=255; mode of extreme value~5;
% npix=770;  ImageSize=20;  w=1;
% pin=[ npix/20, npix/20, npix/6, npix/30, log(1e-2)];
% xD = [0,0]*0.2*npix/ImageSize + npix/2 + 1;
% % xD = [1,1];

global Config

% set up probe points z in pixel
n = 20;  npix=770;
xx = linspace(1,npix,n);  [x1,x2]=meshgrid(xx,xx);  z=[reshape(x1,n*n,1),reshape(x2,n*n,1)];
xn=size(xD,1);  xDw=[];  center=floor(w/2);
for i=0:w-1
for j=0:w-1
    xDw=[xDw; xD+repmat([i,j]-center,xn,1)];
end
end

xDw=round(xDw);
npix=770; xDw(xDw<1)=1; xDw(xDw>npix)=npix;

xDn=size(xDw,1);  yD=zeros(xDn,1);
for i=1:xDn
    yD(i) = yimage(xDw(i,2),xDw(i,1));
end


% initialize guess models
MaL1=pin(1);
MaL2=pin(2);
MbL1=pin(3);
MbL2=pin(4);
Mabn=pin(5);

Config.GuessAcovfunc = {@covSEard};
Config.GuessAlikfunc = {@likGauss};
Config.GuessAhyp.cov = [log(MaL1), log(MaL2), log(1)];
Config.GuessAhyp.lik = Mabn;

Config.GuessBcovfunc = {@covSEard};
Config.GuessBlikfunc = {@likGauss};
Config.GuessBhyp.cov = [log(MbL1), log(MbL2), log(1)];
Config.GuessBhyp.lik = Mabn;

likfunca = Config.GuessAlikfunc;  covfunca = Config.GuessAcovfunc;  hypa = Config.GuessAhyp;
likfuncb = Config.GuessBlikfunc;  covfuncb = Config.GuessBcovfunc;  hypb = Config.GuessBhyp;

Config.ModelTheta = [0,90]/180*pi;
ntheta = length(Config.ModelTheta);
ptheta = ones(1,ntheta)/ntheta;


% initialize things
nn = n*n;
ma = zeros(nn,ntheta);  vara = zeros(nn,ntheta);  nlmla = zeros(1,ntheta);
mb = zeros(nn,ntheta);  varb = zeros(nn,ntheta);  nlmlb = zeros(1,ntheta);
   
    
    for i=1:ntheta

        RotA = [cos(Config.ModelTheta(i)), -sin(Config.ModelTheta(i)); sin(Config.ModelTheta(i)), cos(Config.ModelTheta(i))];
        xDr = (RotA*xDw')';  zr = (RotA*z')';

        % predictive distribution with data hypa or hypb
        [matheta, varatheta] = gp(hypa, @infExact, [], covfunca, likfunca, xDr, yD, zr);
        [mbtheta, varbtheta] = gp(hypb, @infExact, [], covfuncb, likfuncb, xDr, yD, zr);
        ma(:,i) = matheta;  vara(:,i) = varatheta;
        mb(:,i) = mbtheta;  varb(:,i) = varbtheta;

        % marginal likelihood for p(model|D, theta)
        nlmlatheta = gp(hypa, @infExact, [], covfunca, likfunca, xDr, yD);
        nlmlbtheta = gp(hypb, @infExact, [], covfuncb, likfuncb, xDr, yD);
        nlmla(i) = nlmlatheta;  nlmlb(i) = nlmlbtheta;

    end


    % weights of the gaussians for HE:
    nor = sum(ptheta.*exp(-nlmla)) + sum(ptheta.*exp(-nlmlb));
    wa = ptheta.*exp(-nlmla)/nor;
    wb = ptheta.*exp(-nlmlb)/nor;

    % rescale weights to avoid numerical error
    wa=max(1e-15,wa);  wb=max(1e-15,wb);
    nor = sum(wa)+sum(wb);
    wa=wa/nor;  wb=wb/nor;
    
    
    % first term: H(y|x,D)
    H_ME = zeros(nn,1);    
    for i=1:nn
        H_ME(i,1) = GetEntropyGM_mex( [wa,wb], [ma(i,:),mb(i,:)], [vara(i,:),varb(i,:)] );
    end

%     % weights of the gaussians for BALD2
%     nora = sum(ptheta.*exp(-nlmla));
%     norb = sum(ptheta.*exp(-nlmlb));
%     wa2 = ptheta.*exp(-nlmla)/nora;
%     wb2 = ptheta.*exp(-nlmlb)/norb;
%     % make weights non-zero
% %     wa2=max(1e-10,wa2);  wb2=max(1e-10,wb2);
% 
%     % marginal likelihood for p(model|D)
%     mla = sum(ptheta.*exp(-nlmla))/nor;
%     mlb = sum(ptheta.*exp(-nlmlb))/nor;

    mla=sum(wa);  mlb=sum(wb);
    wa2 = wa/mla;  wb2=wb/mlb;

    % second term H of p(y|model, x, D)
    H_BALD2 = zeros(nn,1);    
    for i=1:nn
        H_BALD2(i,1) = mla*GetEntropyGM_mex( wa2, ma(i,:), vara(i,:) ) + mlb*GetEntropyGM_mex( wb2, mb(i,:), varb(i,:) );
    end
    
    H_BALD = H_ME - H_BALD2;
    H_BALD = max(1e-15, H_BALD);  %numerical error ~ at 2e-16
%     H_BALD = (H_BALD-min(H_BALD))/(max(H_BALD)-min(H_BALD))+0.001*(rand(n*n,1));


%     % avoid getting stuck level 1: extreme punishment of visited sites
%     Ind=zeros(xDn,1);
%     for i=1:xDn
%         [~,Ind(i)] = min(sum((z-repmat(xDw(i,:),n*n,1)).^2,2));
%     end
%     H_BALD(Ind) = 0; 

    
%     H_BALD2a = zeros(nn,1);    
%     H_BALD2b = zeros(nn,1);    
%     for i=1:nn
%         H_BALD2a(i,1) = mla*GetEntropyGM_mex( wa2, ma(i,:), vara(i,:) );
%         H_BALD2b(i,1) = mlb*GetEntropyGM_mex( wb2, mb(i,:), varb(i,:) );
%     end
%     H_BALD2=H_BALD2a+H_BALD2b;
    
% YR=[];    
% YR.ma=ma;    
% YR.mb=mb;    
% YR.vara=vara;
% YR.varb=varb;
% YR.wa=wa;
% YR.wb=wb;
% YR.wa2=wa2;
% YR.wb2=wb2;
% YR.mla=mla;
% YR.mlb=mlb;
% YR.H_ME=H_ME;
% YR.H_BALD2a=H_BALD2a;
% YR.H_BALD2b=H_BALD2b;

end

