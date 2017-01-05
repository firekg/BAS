%% I just need one BAS map to start with
% 2015-09-09
% copied from prop_sim

load('expDataWPix.mat', 'SCALwPix');

%pv= [  lsigo,    phit,   lphis, phis0,        lbx,        lbd, lapse, prpa, window, prxpa, prCommon, shift];
pv =[log(1.0),   nan,     nan,   nan,       log(0),        nan,     0,  1/2,      1,   1/3,      1/3,    12];
%pw=[lbm,  km,  nm, bias];
pw=[   0, nan, nan,    0]; % none of the last 3 are used here

% make probe points z
npix=770;
%ImageSize=20;
n=110;
xx = linspace(1,npix,n);
yy = linspace(1,npix,n);
[x1,x2]=meshgrid(xx,yy);
z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
%dxx = xx(2)-xx(1);

AL = a_Dtrials(SCALwPix, 100);
AL.RevealZ = SCALwPix.RevealZ(100,:); % a hack to include RevealZ in a_Dtrials

pv(1) = log(0.5);
pv(6) = 0.34901995;
pv(7) = 0.04471412;
pv(12) = 16;

Gmap = zeros(AL.Trials,25,3);
% BASmap = a_Get_BALDscoreProp(AL,[],pv,pw,randn(1000,25), Gmap, z, 2.6);
% saved as BASmap-2015-09-09

AL.RevealZ(1) = 0.8;
BASmap_0d8 = a_Get_BALDscoreProp(AL,[],pv,pw,zeros(1000,25), Gmap, z, 2.6);
AL.RevealZ(1) = 2.0;
BASmap_2d0 = a_Get_BALDscoreProp(AL,[],pv,pw,zeros(1000,25), Gmap, z, 2.6);
% saved as BASmap-2015-09-14


%% quick analysis
% load('BASmap-2015-09-09.mat');
load('BASmap-2015-09-14.mat');

npix=770;
n=110;
xx = linspace(1,npix,n);
yy = linspace(1,npix,n);
[x1,x2] = meshgrid(xx,yy);
z=[reshape(x1,n*n,1),reshape(x2,n*n,1)];
xdeg = linspace(-13.9, 13.9, n);

pixx = -6:0.1:6;
t = 2;

indm = n*(n/2) + n/2;
indi = rem(indm,n);
indf = indi + n*(n-1);
indset = indi:n:indf;

% length scales
lenS = n/30;
lenM = n/20;
lenL = n/6;
xboxS = [n/2; n/2-lenS; n/2-lenS; n/2];
xboxM = [n/2; n/2-lenM; n/2-lenM; n/2];
xboxL = [n/2; n/2-lenL; n/2-lenL; n/2];
ybox = [0; 0; 1.5; 1.5];

% colormap
clf; 
graymin=0;
graymax=1;
grayres=256;
graywarp = 5;%1/4;
gf=@(graylevel) graylevel.^graywarp;
graylevel=linspace(graymin,graymax,grayres); 
cm=gf(graylevel)'*ones(1,3);  % scale color map
set(gcf,'ColorMap',cm,'PaperPositionMode','auto');
maxdata = 1.5;

% plot BASmap_0d8
BASmap = BASmap_0d8;
score = BASmap.BALDmap(:,:,t);
revx = BASmap.RevealPosX(:);
revy = BASmap.RevealPosY(:);

subplot(2,2,1); hold on;
% colormap(gray)
hold on;
D = score(:);
D2=interp1(linspace(0,1,grayres),(1:grayres),D/maxdata,'nearest');    
D2 = reshape(D2, n, n);
image([1, npix], [1,npix], D2, 'CDataMapping', 'direct');
% plot(revx(1:t-1), revy(1:t-1), 'gx');
plot(z(indset,1), z(indset,2), 'r');
axis([1,npix,1,npix]);
axis off;
axis square;
title('Percieved pixel value = 0.8');

subplot(2,2,3); hold on;
w =  BASmap.ml(:,t,:);  w = w(:);
mu = BASmap.pred_mean(indset,:,1,t);
var = BASmap.pred_var(indset,:,1,t);
nvar = numel(indset);
unexpl_var = zeros(nvar,1);
expl_var = zeros(nvar,1);
% not sure if right... looks alright
for type = 1:3
    wt = w(type);
    unexpl_var = unexpl_var + wt*var(:,type);
    expl_var = expl_var + wt*(1-wt)*mu(:,type).^2;
end
for itype = 1:3
    wi = w(itype);
    if(itype>1)
        for jtype = 1:itype-1            
            wj = w(jtype);
            expl_var = expl_var - 2*wi*wj*mu(:,itype).*mu(:,jtype);            
        end
    end
end
fill(xboxL, ybox, 0.8*[1,1,1], 'FaceAlpha', 1, 'LineStyle','None');
fill(xboxM, ybox, 0.6*[1,1,1], 'FaceAlpha', 1, 'LineStyle','None');
fill(xboxS, ybox, 0.5*[1,1,1], 'FaceAlpha', 1, 'LineStyle','None');
plot(unexpl_var, 'b-');
plot(expl_var, 'r-');
plot(unexpl_var + expl_var, 'k-');
axis([1, n, 0, inf]);
xlabel('Degree');
ylabel('Variance');
mm = (n-1)/27.8;
bb = n - mm*13.9;
set(gca, 'XTick', linspace(-10*mm+bb, 10*mm+bb, 5));
set(gca, 'XTickLabel', {'-10','-5','0','5','10'});
%legend('unexpl.', 'expl.', 'total', 'location', 'east');
axis square

% plot BASmap_2d0 - a copy of the above
BASmap = BASmap_2d0;
score = BASmap.BALDmap(:,:,t);
revx = BASmap.RevealPosX(:);
revy = BASmap.RevealPosY(:);

subplot(2,2,2); hold on;
% colormap(gray)
hold on;
D = score(:);
D2=interp1(linspace(0,1,grayres),(1:grayres),D/maxdata,'nearest');    
D2 = reshape(D2, n, n);
image([1, npix], [1,npix], D2, 'CDataMapping', 'direct');
% plot(revx(1:t-1), revy(1:t-1), 'gx');
plot(z(indset,1), z(indset,2), 'r');
axis([1,npix,1,npix]);
axis off;
axis square;
title('Percieved pixel value = 2.0');
h=gca;
pos=get(h,'Position');
b = colorbar;
set(b,'YTick',[1 grayres],'YTickLabel',{'0',sprintf('%0.3f',maxdata/log(2))});
set(h,'Position',pos);

subplot(2,2,4); hold on;
w =  BASmap.ml(:,t,:);  w = w(:);
mu = BASmap.pred_mean(indset,:,1,t);
var = BASmap.pred_var(indset,:,1,t);

nvar = numel(indset);
unexpl_var = zeros(nvar,1);
expl_var = zeros(nvar,1);
% not sure if right... looks alright
for type = 1:3
    wt = w(type);
    unexpl_var = unexpl_var + wt*var(:,type);
    expl_var = expl_var + wt*(1-wt)*mu(:,type).^2;
end
for itype = 1:3
    wi = w(itype);
    if(itype>1)
        for jtype = 1:itype-1            
            wj = w(jtype);
            expl_var = expl_var - 2*wi*wj*mu(:,itype).*mu(:,jtype);            
        end
    end
end
plot(unexpl_var, 'b-');
plot(expl_var, 'r-');
plot(unexpl_var + expl_var, 'k-');
fill(xboxL, ybox, 0.8*[1,1,1], 'FaceAlpha', 1, 'LineStyle','None');
fill(xboxM, ybox, 0.6*[1,1,1], 'FaceAlpha', 1, 'LineStyle','None');
fill(xboxS, ybox, 0.5*[1,1,1], 'FaceAlpha', 1, 'LineStyle','None');
plot(unexpl_var, 'b-');
plot(expl_var, 'r-');
plot(unexpl_var + expl_var, 'k-');
axis([1, n, 0, inf]);
xlabel('Degree');
ylabel('Variance');
mm = (n-1)/27.8;
bb = n - mm*13.9;
set(gca, 'XTick', linspace(-10*mm+bb, 10*mm+bb, 5));
set(gca, 'XTickLabel', {'-10','-5','0','5','10'});
legend('unexplained', 'explained', 'total', 'location', 'east');
axis square


set(gcf,'Color',[1,1,1]);

%% older version
% function [BASmap] =BASmap()
% %% BASmap has 2 layers: participants x3 (SC,EP,BZ); MC samp x20
%     %pv= [ lsigo,  phit,   lphis, phis0,        lbx,        lbd, lapse, prpa, w];
%     pv =[log(1),   nan,     nan,   nan,        nan,        nan,     0,  0.5, 1];
%     %pw=[lbm,  km,  nm, bias];
%     pw=[   0, nan, nan,    0]; % none of these are used here
%     % make probe points z
%     npix=770;  ImageSize=20;  n=110;  xx = linspace(1,npix,n);  yy = linspace(1,npix,n);
%     [x1,x2]=meshgrid(xx,yy);  z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
% 
%     for i=1:1
%         if(i==1)      
%             load('SCdata.mat');  
%             AL=SCAL;
%             AL=a_Dtrials(AL,100); %just pick out one trial, to check approximation
%             DIMAL=a_DIM(SCAL);  
%             pv(1) = log(1.1);  % from max likelihood fit
%         elseif(i==2)  
%             load('EPdata.mat');  
%             AL=EPAL;  
%             DIMAL=a_DIM(EPAL);  
%             pv(1) = log(1.2);  % from max likelihood fit
%         elseif(i==3)  
%             load('BZdata.mat');  
%             AL=BZAL;  
%             DIMAL=a_DIM(BZAL);  
%             pv(1) = log(0.9); % from max likelihood fit
%         end
%         for isamp = 1:1
%             BASmap(i,isamp) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,2.5);  %BAS
%             %BASmap(i,isamp) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,2.6);  %Max Ent
%             disp(sprintf('Done parti %d samp %d.', i, isamp));
%         end
%     end
% 
% end