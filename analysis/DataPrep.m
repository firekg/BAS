%% simulate to get BAS and pnBAS
% SIM has 2 layers: participants x3 (SC,EP,BZ); conditions x3 (BAS, mBAS, pnBAS)
% PROP has 2 layers: participants x3 (SC,EP,BZ); conditions x2 (ALprop, PLRprop)
%yrM=randn(1000,25);  
Gmap=1;  biasw=zeros(26,1);
%pv= [ lsigo,  phit,   lphis, phis0,        lbx,        lbd, lapse, prpa, w];
%pv=[log(0.3),   100,  log(1),   100,  log(1e-5),  log(1e-5),     0,  0.5, 1]; % v2: simplest perception noise
%pv=[log(2.2),   100,  log(1),   100,  log(1e-5),  log(1e-5),     0,  0.5, 1]; % v3
pv =[log(1),   nan,     nan,   nan,        nan,        nan,     0,  0.5, 1]; % v4
%pw=[lbm,  km,  nm, bias];
pw=[   0, nan, nan,    0]; %none of the last 3 are used in current version
% make probe points z
npix=770;  ImageSize=20;  n=110;  xx = linspace(1,npix,n);  yy = linspace(1,npix,n);
[x1,x2]=meshgrid(xx,yy);  z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
dxx = xx(2)-xx(1);

for i=1:3
    if(i==1)      load('SCdata.mat');  AL=SCAL;  PLR=SCPLRbt;  DIMAL=a_DIM(SCAL);  DIMPLR=a_DIM(SCPLRbt);
    elseif(i==2)  load('EPdata.mat');  AL=EPAL;  PLR=EPPLRbt;  DIMAL=a_DIM(EPAL);  DIMPLR=a_DIM(EPPLRbt);
    elseif(i==3)  load('BZdata.mat');  AL=BZAL;  PLR=BZPLRbt;  DIMAL=a_DIM(BZAL);  DIMPLR=a_DIM(BZPLRbt);
    end
    pw(1)=inf;     SIM(i,1) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,3);  %BAS
    pw(1)=inf;     SIM(i,2) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,4,[],biasw);  %mBAS = BAS + motor noise
    pw(1)=log(1);  SIM(i,3) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),[],z,4,[],biasw);  %nBAS = BAS + motor noise + selection noise
    disp(sprintf('parti %d SIM done.',i));
    PROP(i,1) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %AL prop
    PROP(i,2) = a_Get_BALDscoreProp(PLR,DIMPLR,pv,pw,randn(1000,25),Gmap,z,2); %PLR prop
    disp(sprintf('parti %d PROP done.',i));
end
save('SIM10.mat','SIM');
save('PROP10.mat','PROP');

%% make a matrix of structure to store all data
% REV has 3 layers: 
% patterns x4 (PA,SH,SV,ALL); 
% participants x4 (SC,EP,BZ,AVG);
% conditions x3 (AL, BAS, pnBAS)
load('SIMf.mat');

load('SCdata.mat');  
REV(1,1,1)=a_DimageID(SCAL,20);       REV(2,1,1)=a_DimageID(SCAL,6);       REV(3,1,1)=a_DimageID(SCAL,30);  
REV(1,1,2)=a_DimageID(SIM(1,1),20);   REV(2,1,2)=a_DimageID(SIM(1,1),6);   REV(3,1,2)=a_DimageID(SIM(1,1),30);  
REV(1,1,3)=a_DimageID(SIM(1,3),20);   REV(2,1,3)=a_DimageID(SIM(1,3),6);   REV(3,1,3)=a_DimageID(SIM(1,3),30);  
load('EPdata.mat');
REV(1,2,1)=a_DimageID(EPAL,20);       REV(2,2,1)=a_DimageID(EPAL,6);       REV(3,2,1)=a_DimageID(EPAL,30);  
REV(1,2,2)=a_DimageID(SIM(2,1),20);   REV(2,2,2)=a_DimageID(SIM(2,1),6);   REV(3,2,2)=a_DimageID(SIM(2,1),30);  
REV(1,2,3)=a_DimageID(SIM(2,3),20);   REV(2,2,3)=a_DimageID(SIM(2,3),6);   REV(3,2,3)=a_DimageID(SIM(2,3),30);  
load('BZdata.mat');
REV(1,3,1)=a_DimageID(BZAL,20);       REV(2,3,1)=a_DimageID(BZAL,6);       REV(3,3,1)=a_DimageID(BZAL,30);  
REV(1,3,2)=a_DimageID(SIM(3,1),20);   REV(2,3,2)=a_DimageID(SIM(3,1),6);   REV(3,3,2)=a_DimageID(SIM(3,1),30);  
REV(1,3,3)=a_DimageID(SIM(3,3),20);   REV(2,3,3)=a_DimageID(SIM(3,3),6);   REV(3,3,3)=a_DimageID(SIM(3,3),30);  

% combining pattern as the 4th pattern
for j=1:3
for k=1:3
    zz1=REV(1,j,k);    zz2=REV(2,j,k);    zz3=REV(3,j,k);
    REV(4,j,k)=a_Combine_Data(zz1,zz2,zz3);    
end
end
% combining participants, the 4th participant
for i=1:4
for k=1:3
    REV(i,4,k)=a_Combine_Data(REV(i,1,k),REV(i,2,k),REV(i,3,k));    
end
end
save('REVf.mat','REV');

%% caluclate revealing densities (difference and mean)
% define gaussian filter
npix=770; ImageSize=20; bin=1;
GaussFilter=a_GaussFilter(npix);

revMap=zeros(770,770,4,4,3);
for k=1:3
for j=1:4
for i=1:4
    zz=REV(i,j,k);
    revMap(:,:,i,j,k) = a_Get_RevealMap_mod2(zz,1:25,1:zz.Trials,GaussFilter,770,2);
%     imagesc(revMap(:,:,i,j,k));
%     pause
end
end
end

% caluclate difference densities
drevMap=zeros(770,770,3,4,3);
for k=1:3
for j=1:4
for i=1:3
    drevMap(:,:,i,j,k) = revMap(:,:,i,j,k)-revMap(:,:,4,j,k);
    %imagesc(drevMap(:,:,i,j,k));
    %pause
end
end
end

% scale maps to [-1,1]
drevV=abs(reshape(drevMap,numel(drevMap),1));
drevMap=drevMap/max(drevV);
revV=abs(reshape(revMap,numel(revMap),1));
revMap=revMap/max(revV);

%% caluclate revealing densities with revealing number (difference and mean)
% define gaussian filter
npix=770/2; ImageSize=20; bin=1;
GaussFilter=a_GaussFilter(npix);

RrevMap=zeros(25,npix,npix,4,4,3);
for r=1:25
for k=1:3
for j=1:4
for i=1:4
    zz=REV(i,j,k);
    RrevMap(r,:,:,i,j,k) = a_Get_RevealMap_mod2(zz,1:r,1:zz.Trials,GaussFilter,npix,2);
%     imagesc(revMap(:,:,i,j,k));
%     pause
end
end
end
end

% caluclate difference densities
RdrevMap=zeros(25,npix,npix,3,4,3);
for r=1:25
for k=1:3
for j=1:4
for i=1:3
    RdrevMap(r,:,:,i,j,k) = RrevMap(r,:,:,i,j,k)-RrevMap(r,:,:,4,j,k);
    %imagesc(drevMap(:,:,i,j,k));
    %pause
end
end
end
end

% scale maps to [-1,1]
RdrevV=abs(reshape(RdrevMap,numel(RdrevMap),1));
v=max(RdrevV);  RdrevMap=RdrevMap/v;
RrevV=abs(reshape(RrevMap,numel(RrevMap),1));
v=max(RrevV);  RrevMap=RrevMap/max(RrevV);

%% calculate similarity measure (LONG!!)
load('REVC.mat');

% define Gaussian window again
npix=770; ImageSize=20; bin=10;
var=(20/bin)^2; %var=(npix/ImageSize/2)^2;
n=sqrt(var)*6;
[x1,x2]=meshgrid(1:n,1:n);
z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % in pixel coordinate
GaussFilter = zeros(n*n,1);
for i=1:n*n
    GaussFilter(i) = mvnpdf(z(i,:),[1,1]*n/2,[1,0;0,1]*var);
end
GaussFilter = reshape(GaussFilter,n,n);

% start bootstrap loop
nmap=npix/bin; %P and Q map size
nrev=25;  %rev=1:25;  %nrev & rev should be defined together: if nrev>1 --> rev=scalar; if nrev=1, rev can be vector
nboot=10; %number of bootstrap samples
SMboot=zeros(nrev,nboot,6,3,4);

for parti=1:4 %in total: 1:4
for pq=1:3 %cases: human-human, human-optimal, human-noisy_BAS    
    if(pq==1)
        % Q1=P1;  Q2=P2;  Q3=P3;  %Qav=Pav;
        %avoid correlation at 1st revealing
        P1=REV(1,parti,1);  P2=REV(2,parti,1);  P3=REV(3,parti,1);  %Pav=REV(4,parti,1);
        Q1=a_Dtrials(P1, 1:2:P1.Trials);
        Q2=a_Dtrials(P2, 1:2:P2.Trials);
        Q3=a_Dtrials(P3, 1:2:P3.Trials);
        P1=a_Dtrials(P1, 2:2:P1.Trials);
        P2=a_Dtrials(P2, 2:2:P2.Trials);
        P3=a_Dtrials(P3, 2:2:P3.Trials);
    elseif(pq==2)
        P1=REV(1,parti,1);  P2=REV(2,parti,1);  P3=REV(3,parti,1);  %Pav=REV(4,parti,1);
        Q1=REV(1,parti,2);  Q2=REV(2,parti,2);  Q3=REV(3,parti,2);  %Qav=REV(4,parti,2);  %BAS
    elseif(pq==3)
        P1=REV(1,parti,1);  P2=REV(2,parti,1);  P3=REV(3,parti,1);  %Pav=REV(4,parti,1);
        Q1=REV(1,parti,3);  Q2=REV(2,parti,3);  Q3=REV(3,parti,3);  %Qav=REV(4,parti,3);  %noisy BAS   
    end
    for zz=1:6
        SMpq=zeros(nrev,nboot);
        for irev=1:nrev %number of revealings
            rev=1:irev;
            for iboot=1:nboot %number of bootstraps
                % bootstrapped sample
                P1map = a_Get_RevealMap_mod2(P1,rev,1:P1.Trials,GaussFilter,nmap,1);
                P2map = a_Get_RevealMap_mod2(P2,rev,1:P2.Trials,GaussFilter,nmap,1);
                P3map = a_Get_RevealMap_mod2(P3,rev,1:P3.Trials,GaussFilter,nmap,1);
                Pavmap = (P1map+P2map+P3map)/3;
                %Pavmap = a_Get_RevealMap_mod2(Pav,rev,1:Pav.Trials,GaussFilter,nmap,1);
                Q1map = a_Get_RevealMap_mod2(Q1,rev,1:Q1.Trials,GaussFilter,nmap,1);                
                Q2map = a_Get_RevealMap_mod2(Q2,rev,1:Q2.Trials,GaussFilter,nmap,1);                
                Q3map = a_Get_RevealMap_mod2(Q3,rev,1:Q3.Trials,GaussFilter,nmap,1);
                Qavmap = (Q1map+Q2map+Q3map)/3;
                %Qavmap = a_Get_RevealMap_mod2(Qav,rev,1:Qav.Trials,GaussFilter,nmap,1);                
                % mean-corrected
                P1map=P1map-Pavmap;  P2map=P2map-Pavmap;  P3map=P3map-Pavmap;
                Q1map=Q1map-Qavmap;  Q2map=Q2map-Qavmap;  Q3map=Q3map-Qavmap;
                % P-Q 6 cases: PA-PA, SH-SH, SV-SV, PA-SH, PA-SV, SH-SV
                if(zz==1)      P=P1map;  Q=Q1map;
                elseif(zz==2)  P=P2map;  Q=Q2map;
                elseif(zz==3)  P=P3map;  Q=Q3map;
                elseif(zz==4)  P=P1map;  Q=Q2map;
                elseif(zz==5)  P=P1map;  Q=Q3map;
                elseif(zz==6)  P=P2map;  Q=Q3map;
                end                
                P=reshape(P,numel(P),1);  Q=reshape(Q,numel(Q),1);
                % calculate JS divergence
                    %A=0.5*P+0.5*Q;
                    %SMpq(irev,iboot) = 0.5*KLDiv(P,A)+0.5*KLDiv(Q,A);
                % Or mean square error
                    %SMpq(irev,iboot) = mean((P-Q).^2);
                % Or variance explained
                    %dP=P-mean(P);  dQ=Q-mean(Q);
                    %SMpq(irev,iboot) = sum(dP.*dQ)^2/sum(dP.^2)/sum(dQ.^2);
                % Or correlation
                    c=corrcoef(P,Q); %2-by-2
                    SMpq(irev,iboot) =c(1,2);
            end
        end
    %SMpqm=mean(SMpq,1); %averaging out revealing number
    SMboot(:,:,zz,pq,parti)=SMpq;
    disp(sprintf('parti=%d; pq=%d; zz=%d.',parti,pq,zz));    
    end
end
end
save('SMbootC.mat','SMboot');
% SMboot has 5 layers:
% revealing number 25x: 1:25
% bootstrap samples bootx: 1:nboot
% pattern (mis)match 6x: PA-PA, SH-SH, SV-SV, PA-SH, PA-SV, SH-SV
% conditions 3x: self-self, self-BAS, self-pnBAS
% participants 4x: SC, EP, BZ, AVG

%% Compute performance
% PERF has 4 layers: revealing number x5(1:5:25); parti x4; condiitons x3 (AL,PLR,PLB); moments x2 (mean,SD)
PERF=zeros(5,4,3,2);
for i=1:4
    if(i==1)      load('SCdata.mat');  AL=SCAL;  PLR=SCPLRbt;  PLB=SCPLBbt;
    elseif(i==2)  load('EPdata.mat');  AL=EPAL;  PLR=EPPLRbt;  PLB=EPPLBbt;
    elseif(i==3)  load('BZdata.mat');  AL=BZAL;  PLR=BZPLRbt;  PLB=BZPLBbt;
    elseif(i==4)  
        AL=a_Combine_Data(SCAL,EPAL,BZAL);
        PLR=a_Combine_Data(SCPLRbt,EPPLRbt,BZPLRbt);
        PLB=a_Combine_Data(SCPLBbt,EPPLBbt,BZPLBbt);
    end    
    [PERF(:,i,1,1), PERF(:,i,1,2), ~, ~] = a_Get_Human_Perf(AL);
    [PERF(:,i,2,1), PERF(:,i,2,2), ~, ~] = a_Get_Human_Perf(PLR);
    [PERF(:,i,3,1), PERF(:,i,3,2), ~, ~] = a_Get_Human_Perf(PLB);
end
save('PERF.mat','PERF');

%% Compute information gain
% SIM has 2 layers: participants x3 (SC,EP,BZ); conditions x3 (BAS, mBAS, pnBAS)
% PROP has 2 layers: participants x3 (SC,EP,BZ); conditions x2 (ALprop, PLRprop)
% INFO has 4 layers: revealing number x25 (1:25); parti x4; condiitons x5 (AL,PLR,BAS,pnBAS,mBAS); moments x2 (mean,SD)
load('SIM.mat');
load('PROP.mat');
INFO=zeros(25,4,4,2);
for i=1:3
    m=PROP(i,1).mla;  M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  M= -M/log(0.5);
    INFO(:,i,1,1)=nanmean(M);  INFO(:,i,1,2)= sqrt(nanvar(M)./sum(~isnan(M)));
    m=PROP(i,2).mla;  M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  M= -M/log(0.5);
    INFO(:,i,2,1)=nanmean(M);  INFO(:,i,2,2)= sqrt(nanvar(M)./sum(~isnan(M)));
    m=SIM(i,1).mla;  M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  M= -M/log(0.5);
    INFO(:,i,3,1)=nanmean(M);  INFO(:,i,3,2)= sqrt(nanvar(M)./sum(~isnan(M)));
    m=SIM(i,3).mla;  M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  M= -M/log(0.5);
    INFO(:,i,4,1)=nanmean(M);  INFO(:,i,4,2)= sqrt(nanvar(M)./sum(~isnan(M)));
    m=SIM(i,2).mla;  M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  M= -M/log(0.5);
    INFO(:,i,5,1)=nanmean(M);  INFO(:,i,5,2)= sqrt(nanvar(M)./sum(~isnan(M)));
    disp(sprintf('i=%d;',i));
end
for i=4
    m=[PROP(1,1).mla; PROP(2,1).mla; PROP(3,1).mla];     
    M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  M= -M/log(0.5);
    INFO(:,i,1,1)=nanmean(M);  INFO(:,i,1,2)= sqrt(nanvar(M)./sum(~isnan(M)));
    m=[PROP(1,2).mla; PROP(2,2).mla; PROP(3,2).mla];     
    M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  M= -M/log(0.5);
    INFO(:,i,2,1)=nanmean(M);  INFO(:,i,2,2)= sqrt(nanvar(M)./sum(~isnan(M)));
    m=[SIM(1,1).mla; SIM(2,1).mla; SIM(3,1).mla];     
    M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  M= -M/log(0.5);
    INFO(:,i,3,1)=nanmean(M);  INFO(:,i,3,2)= sqrt(nanvar(M)./sum(~isnan(M)));
    m=[SIM(1,3).mla; SIM(2,3).mla; SIM(3,3).mla];     
    M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  M= -M/log(0.5);
    INFO(:,i,4,1)=nanmean(M);  INFO(:,i,4,2)= sqrt(nanvar(M)./sum(~isnan(M)));
    m=[SIM(1,2).mla; SIM(2,2).mla; SIM(3,2).mla];     
    M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  M= -M/log(0.5);
    INFO(:,i,5,1)=nanmean(M);  INFO(:,i,5,2)= sqrt(nanvar(M)./sum(~isnan(M)));
    disp(sprintf('i=%d;',i));
end
save('INFOf.mat','INFO');

%% Compute Correlation curves
% CORR has 5 layers:
% revealing number 25x: 1:25
% moments x3: (mean, upper 95%, lower 95%)
% pattern (mis)match 2x: same, different
% conditions 3x: self-self, self-BAS, self-pnBAS
% participants 4x: SC, EP, BZ, AVG
load('SMboot4.mat');
nboot=100;
temp_same=mean(SMboot(:,:,1:3,:,:),3);  temp_same = sort(temp_same,2);
temp_diff=mean(SMboot(:,:,4:6,:,:),3);  temp_diff = sort(temp_diff,2);
CORR(:,1,1,:,:) = mean(temp_same,2);
CORR(:,2,1,:,:) = temp_same(:,round(nboot*0.95),1,:,:);
CORR(:,3,1,:,:) = temp_same(:,round(nboot*0.05),1,:,:);
CORR(:,1,2,:,:) = mean(temp_diff,2);
CORR(:,2,2,:,:) = temp_diff(:,round(nboot*0.95),1,:,:);
CORR(:,3,2,:,:) = temp_diff(:,round(nboot*0.05),1,:,:);
save('CORR4.mat','CORR');





