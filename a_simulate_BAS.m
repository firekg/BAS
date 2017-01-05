load('SCdata.mat');
DIMSCAL=a_DIM(SCAL);
DIMSCPLRbt=a_DIM(SCPLRbt);
yrM=0*rand(1000,25);

%% initialize parameters
% Just the simplest perception noise
%pv= [ lsigo,  phit,   lphis, phis0,        lbx,        lbd, lapse, prpa, w];
pv=[log(0.3),   100,  log(1),   100,  log(1e-5),  log(1e-5),     0,  0.5, 1];
%pw=[lbm, km,    nm, bias];
pw=[log(1),  0,  0.01,    0];

%% make z
npix=770;  ImageSize=20;  n=120;
xx = linspace(1,npix,n);
yy = linspace(1,npix,n);
[x1,x2]=meshgrid(xx,yy);
z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
dxx = xx(2)-xx(1);

%% create Gmap
nm=pw(3); 
mu = [npix, npix]./2;
covar = [nm^2, 0; 0, nm^2];
F = mvnpdf([x1(:) x2(:)],mu,covar); F = reshape(F,length(x2),length(x1));
ng=min(floor(n/2),ceil(nm/dxx*5));
Gmap = F(round(n/2-ng):round(n/2+ng),round(n/2-ng):round(n/2+ng));
ss=sum(sum(Gmap));
if(ss==0)
    Gmap=1;
else
    Gmap = Gmap/sum(sum(Gmap)); %normalize
end
%% simulate BAS and mBAS (BAS with motor noise) on SCAL trials
ALbas = a_Get_BALDscoreProp(SCAL,DIMSCAL,pv,pw,yrM,Gmap,z,3);
ALmbas = a_Get_BALDscoreProp(SCAL,DIMSCAL,pv,pw,yrM,Gmap,z,4);
%% computing bias map
% get REV from draft_fig_rev first
mapsize=120;  GaussFilter=a_GaussFilter(120);  zz=REV(4,1,1);  
PbiasM=zeros(mapsize^2,25);  biasw=zeros(25,1); 
for i=1:25
    Pparti=a_Get_RevealMap_mod2(zz,1:i,1:zz.Trials,GaussFilter,mapsize,2);
    Pnbas=a_Get_RevealMap_mod2(ALmbas,1:i,1:ALmbas.Trials,GaussFilter,mapsize,2);
    q=max(max((Pnbas./Pparti)));
    Pbias=q*Pparti-Pnbas;  Pbias=Pbias/sum(Pbias(:));
    Pbias=reshape(Pbias',mapsize^2,1);  %to match the convention of probe points z
    PbiasM(:,i)=Pbias;  biasw(i)=(q-1)/q;
    disp(sprintf('i=%d;',i));
end
%pw(4)=(q-1)/q;  pw(4)=0.5;
%subplot(2,2,1); imagesc(Pnbas);  subplot(2,2,2); imagesc(Pparti);
%subplot(2,2,3); imagesc(Pnbas./Pparti);  subplot(2,2,4); imagesc(Pbias);
%% simulate nBAS (noisy BAS with bias map)
biasw=zeros(25,1);
SCALnbas = a_Get_BALDscoreProp(SCAL,DIMSCAL,pv,pw,yrM,Gmap,z,4,PbiasM,biasw);  % motor noise + bias

%% compute BAS prop on SCAL
SCALprop = a_Get_BALDscoreProp(SCAL,DIMSCAL,pv,pw,yrM,Gmap,z,2);
SCPLRprop = a_Get_BALDscoreProp(SCPLRbt,DIMSCPLRbt,pv,pw,yrM,Gmap,z,2);

%% BALD score Analysis part 5: entropy of whole system
m=SCALprop.mla;
M = -log(0.5) +m.*log(m) +(1-m).*log(1-m); M= -M/log(0.5);
ALentr_m = nanmean(M);  ALentr_v = nanvar(M);  ALentr_n = sum(~isnan(M));

m=SCPLRprop.mla;
M = -log(0.5) +m.*log(m) +(1-m).*log(1-m); M= -M/log(0.5);
PLRentr_m = nanmean(M);  PLRentr_v = nanvar(M);  PLRentr_n = sum(~isnan(M));

m=ALbas.mla;
M = -log(0.5) +m.*log(m) +(1-m).*log(1-m); M= -M/log(0.5);
BASentr_m = nanmean(M);  BASentr_v = nanvar(M);  BASentr_n = sum(~isnan(M));

m=SCALnbas.mla;
M = -log(0.5) +m.*log(m) +(1-m).*log(1-m); M= -M/log(0.5);
nBASentr_m = nanmean(M);  nBASentr_v = nanvar(M);  nBASentr_n = sum(~isnan(M));

Bx=1:25;
errorbar(Bx, BASentr_m, sqrt(BASentr_v./BASentr_n), 'k-d','MarkerFaceColor','k'); hold on;
errorbar(Bx, nBASentr_m, sqrt(nBASentr_v./nBASentr_n), 'd-','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5]);
errorbar(Bx, ALentr_m, sqrt(ALentr_v./ALentr_n), 'r-o','MarkerFaceColor','r'); hold on;
errorbar(Bx, PLRentr_m, sqrt(PLRentr_v./PLRentr_n), 'b-s','MarkerFaceColor','b'); hold off;
axis([0 25 0 1]);  grid on; %axis square;  
set(gca,'Ytick',0:0.2:1);  ylabel('Entropy reduction', 'FontSize', 10);
set(gca,'Xtick',0:5:25);  xlabel('Revealing number', 'FontSize', 10);
leg=legend('BAS', 'biased nBAS', 'free-scan', 'random', 'Location','SouthEast'); 
set(leg,'Color','w');

%%
set(gcf, 'Color', [1,1,1]);
export_fig('EntropyReduction.pdf')


