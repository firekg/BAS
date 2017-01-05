function plot_BASmap(BASmap)
% files in maps.mat

nsamp = 1;
%find max rev for each trial
MaxRev = sum(BASmap(1,1).RevealPosX~=0, 2);

% find best matching BASmap with PTILE
%Meas3 = nan(600,25,nsamp);
for isamp=1:nsamp
    Meas3(:,:,isamp) = BASmap(1,isamp).PTILExs;
    %Meas3(:,:,isamp) = BASmap(1,isamp).opt_xs;
end
Meas2 = nanmean(Meas3,3);
Meas1 = nanmean(Meas2,2);  %plot(Ptile1)
Meas1(MaxRev<=5) = 0; % only consider more than 5 moves
[~,ind] = max(Meas1);

% choose a random ind
trInd=100; %a good one, rev num=15, SV
mr = MaxRev(trInd);

% get the right BASmaps
BALDmap4 = nan(1,12100,mr,nsamp);
%BALDmap4diff = nan(1,12100,mr,nsamp);
ML4 = nan(1,mr,3,nsamp); %this one doesn't make sense with MC averaging
for isamp=1:nsamp
    BALDmap4(1,:,:,isamp) = BASmap(1,isamp).BALDmap(trInd,:,1:mr);
    %BALDmap4diff(1,:,:,isamp) = BASmap(1,1).BALDmap(ind,:,1:mr) - BASmap(1,isamp).BALDmap(ind,:,1:mr);
    ML4(1,:,:,isamp) = BASmap(1,isamp).ml(trInd,1:mr,:);
end
BALDmap3 = nanmean(BALDmap4,4);
%BALDmap3diff = nanmean(BALDmap4diff,4);
ML3 = nanmean(ML4,4);

% plotting the map
npix=770;  ImageSize=20;  n=110;  xx = linspace(1,npix,n);  yy = linspace(1,npix,n);
[x1,x2]=meshgrid(xx,yy);  z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate

px = linspace(1,npix,n);

xDx = BASmap(1,1).RevealPosX(trInd,1:mr+1);
xDy = BASmap(1,1).RevealPosY(trInd,1:mr+1);

%load underlying image
load('SCdata');
D = SCAL;
Type = D.ImageID(trInd,1);
X    = D.ImageID(trInd,2);    
Y    = D.ImageID(trInd,3);
Num  = D.ImageID(trInd,4);
%s=['im\Type1X20Y20Num',num2str(jj)];
ims=strcat('im\Type',num2str(Type),'X',num2str(X),'Y',num2str(Y),'Num', num2str(Num));

figure(4)
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
set(gcf,'Position',[1,2,20,20]);
set(gcf, 'PaperSize', [20,20]);
set(gcf, 'PaperPosition', [1,2,20,20]);
% colormap(gray);

% Mate's grayscale
graymin=0;
graymax=1;
grayres=100;
% graywarp=1/4; %BASmap
% gf=@(D) (D/0.0117).^graywarp; % BASmap
graywarp=1; %MaxEntMap
gf=@(D) (D/1.71).^graywarp; %MaxEntMap
cms=linspace(graymin,graymax,grayres); 
cm=cms'*ones(1,3);
set(gcf,'ColorMap',cm);

pd=0.005; %distance between plots
ps=0.15; %plot size
lm=0.02; %left margin
ls1=lm+2*pd; %1st column left start position
ts1=0.8; %1st row top start position

for jj=1:3
for ii=1:5
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]
    rev = (jj-1)*5 + ii;
    if(rev==1)
        img=imread(ims,'bmp');
        image(img);
    else
        Bmap = BALDmap3(1,:,rev);
        
        D2=interp1(linspace(0,1,grayres),(1:grayres),gf(Bmap),'nearest'); % scale score
        D2=reshape(D2,n,n);
        image(px,px,D2,'CDataMapping','direct');
        
        %imagesc(px, px, reshape(log(Bmap),n,n));  %for BAS
        %imagesc(px, px, reshape(Bmap,n,n));  %for MAx Ent
        
        hold on;
%         % plot fixations
%         for r=1:rev-1
%             %cl = (max(-4,(r-rev))+5)/5;  %this is the fading
%             plot( xDx(r), xDy(r), '.', 'MarkerEdgeColor',[0,0.6,0]);
%         end
%         plot( xDx(rev), xDy(rev), 'y.');
%         %plot( xDx(rev-1:rev), xDy(rev-1:rev), 'y-'); %line distracting        
%         % plot best location
%         maxB = max(Bmap);
%         indmaxB = find(Bmap == maxB);
%         if(numel(indmaxB)<10)
%             for r = 1:numel(indmaxB)
%                 indz = indmaxB(r); 
%                 plot(z(indz,1), z(indz,2), 'x', 'MarkerSize', 3, 'MarkerEdgeColor', [0,0,1]);
%             end
%         end
%         title(sprintf('%d%%', floor(Meas2(trInd,rev))),'Color', 'k')
        axis square;
    end
    set(gca,'Xtick',[],'Ytick',[]);
end
end

%legned bar mean
jj=4; ii=1;
ls=ls1+(ii-1)*ps+(ii-1)*pd;
ts=ts1-(jj-1)*ps-(jj-1)*pd;
subplot('Position',[ls, ts+4*(ps-pd/2)/5, ps-pd/2, (ps-pd/2)/5]); %[jj,ii]    
imagesc(gf(linspace(0,1,1000))); 
set(gca,'Xtick',[],'Ytick',[]);

% restart figure spacing for likelihood
pd=0.01; %distance between plots
ps=0.05; %plot size
lm=0.03; %left margin
ls1=lm+2*pd; %1st column left start position
ts1=0.2; %1st row top start position

% for jj=1:3
% for ii=1:5
%     ls=ls1+(ii-1)*ps+(ii-1)*pd;
%     ts=ts1-(jj-1)*ps-(jj-1)*pd;
%     subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]
%     rev = (jj-1)*5 + ii;
%     bar(1, ML3(1,rev,1), 'FaceColor', 'r'); hold on;
%     bar(2, ML3(1,rev,2), 'FaceColor', 'g');
%     bar(3, ML3(1,rev,3), 'FaceColor', 'b');
%     plot(0:4, 0.5*ones(5,1), 'k:');
%     axis([0 4 0 1]);  %axis square;  %grid on;
%     set(gca,'Xtick',[],'Ytick',[],'box','off');
% end
% end
    
% % save figure
set(gcf,'Color',[1,1,1]);
% export_fig('BASMaps.pdf')
export_fig('MaxEntMaps.pdf')


end