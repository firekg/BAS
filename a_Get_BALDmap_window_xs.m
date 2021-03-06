function [llslpa, llsls, lH_BALD, lH_BALDxs, z] = a_Get_BALDmap_window_xs(xDxs, yimage, pin, w, ppa, b1, lapse, n, zxlim, zylim)
% n is number rows and cols of probe-point matrix z
% xD in pixel coordinate
% yimage in pixel coordinate

% split xDxs into xD and x*
xDxsn=size(xDxs,1);
xD=xDxs(1:xDxsn-1,:);
xs=xDxs(xDxsn,:);

% set up probe points z in pixel
xx = linspace(zxlim(1),zxlim(2),n);  %1:npix
yy = linspace(zylim(1),zylim(2),n);  %1:npix
[x1,x2]=meshgrid(xx,yy);
z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % in pixel x,y coordinate
z=[z;xs];

% get data with window
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
    yD(i) = yimage(xDw(i,2),xDw(i,1)); % in pixel row,col coordinate: x is column in image, y is row in image;
end

% initialize guess models
MaL1=pin(1);
MaL2=pin(2);
MbL1=pin(3);
MbL2=pin(4);
Mabn=pin(5);

pahyp.cov = [log(MaL1), log(MaL2), log(1)];  pahyp.lik = Mabn;  %ppa=0.5;
shhyp.cov = [log(MbL1), log(MbL2), log(1)];  shhyp.lik = Mabn;  psh=(1-ppa)/2;
svhyp.cov = [log(MbL2), log(MbL1), log(1)];  svhyp.lik = Mabn;  psv=(1-ppa)/2;

covfunc = {@covSEard};  likfunc = {@likGauss};

% marginal likelihood for p(model|D, theta)
lmlpa = -gp(pahyp, @infExact, [], covfunc, likfunc, xDw, yD);
lmlsh = -gp(shhyp, @infExact, [], covfunc, likfunc, xDw, yD);
lmlsv = -gp(svhyp, @infExact, [], covfunc, likfunc, xDw, yD);

[mpa, vpa] = gp(pahyp, @infExact, [], covfunc, likfunc, xDw, yD, z);
[msh, vsh] = gp(shhyp, @infExact, [], covfunc, likfunc, xDw, yD, z);
[msv, vsv] = gp(svhyp, @infExact, [], covfunc, likfunc, xDw, yD, z);

lml = max([lmlpa,lmlsh,lmlsv]);
dlmlpa = lmlpa-lml;
dlmlsh = lmlsh-lml;
dlmlsv = lmlsv-lml;

lnor = log(ppa*exp(dlmlpa)+psh*exp(dlmlsh)+psv*exp(dlmlsv));

lwpa = log(ppa)+dlmlpa-lnor;
lwsh = log(psh)+dlmlsh-lnor;
lwsv = log(psv)+dlmlsv-lnor; %exp(lwpa)+exp(lwsh)+exp(lwsv)
lws = log(exp(lwsh)+exp(lwsv)); %log weight stripy

lpr = lwpa -lws; %log posterior ratio
softlik = 1/(1+exp(-b1*lpr));  % softmax of lpr with b1
lapse_softlik = (1-lapse)*softlik + 0.5*lapse;

llslpa = log(lapse_softlik); % log lapse soft likelihood for pa
llsls = log(1-lapse_softlik); % log lapse soft likelihood for stripy
llslsh = llsls + lwsh - lws;
llslsv = llsls + lwsv - lws;

lw  = [llslpa,llslsh,llslsv];  %sum(exp(lw))
mu  = [mpa;msh;msv];
var = [vpa;vsh;vsv];

H_BALD = a_GetBALDGM_mex( lw, mu, var );
lH_BALD = log(max(0, H_BALD(1:end-1)));  % numerical error ~ at 2e-16
lH_BALDxs = log(max(0, H_BALD(end)));

end

