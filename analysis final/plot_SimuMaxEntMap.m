function plot_MaxEntmap(D)

nsamp = 1;
%find max rev for each trial
MaxRev = sum(D(1,1).RevealPosX~=0, 2);

% choose a random ind
trInd=100; %a good one, rev num=15, SV
mr = MaxRev(trInd);

% get the right BASmaps
map4 = nan(1,12100,mr,nsamp);
ML4 = nan(1,mr,3,nsamp); %this one doesn't make sense with MC averaging
for isamp=1:nsamp
    map4(1,:,:,isamp) = D(1,isamp).map(trInd,:,1:mr);
    ML4(1,:,:,isamp) = D(1,isamp).ml(trInd,1:mr,:);
end
map3 = nanmean(map4,4);
ML3 = nanmean(ML4,4);

% plotting the map
npix=770;  ImageSize=20;  n=110;  xx = linspace(1,npix,n);  yy = linspace(1,npix,n);
[x1,x2]=meshgrid(xx,yy);  z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate

px = linspace(1,npix,n);

xDx = D(1,1).RevealPosX(trInd,1:mr);
xDy = D(1,1).RevealPosY(trInd,1:mr);

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
        Bmap = map3(1,:,rev-1);
        colormap(gray);
        imagesc(px, px, reshape(Bmap,n,n));
        hold on;
        % plot fixations
        for r=1:rev-1
            %cl = (max(-4,(r-rev))+5)/5;  %this is the fading
            plot( xDx(r), xDy(r), '.', 'MarkerEdgeColor',[0,0.6,0]);
        end
        plot( xDx(rev), xDy(rev), 'y.');
        %plot( xDx(rev-1:rev), xDy(rev-1:rev), 'y-'); %line distracting        
        % plot best location
        axis square;
        title(sprintf('%d', rev-1),'Color','k')
    end
    set(gca,'Xtick',[],'Ytick',[]);
end
end

% %legned bar mean
% jj=4; ii=1;
% ls=ls1+(ii-1)*ps+(ii-1)*pd;
% ts=ts1-(jj-1)*ps-(jj-1)*pd;
% subplot('Position',[ls, ts+4*(ps-pd/2)/5, ps-pd/2, (ps-pd/2)/5]); %[jj,ii]    
% imagesc(1:1000); 
% set(gca,'Xtick',[],'Ytick',[]);

% restart figure spacing for likelihood
pd=0.01; %distance between plots
ps=0.05; %plot size
lm=0.03; %left margin
ls1=lm+2*pd; %1st column left start position
ts1=0.2; %1st row top start position

for jj=1:3
for ii=1:5
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]
    rev = (jj-1)*5 + ii;
    bar(1, ML3(1,rev,1), 'FaceColor', 'r'); hold on;
    bar(2, ML3(1,rev,2), 'FaceColor', 'g');
    bar(3, ML3(1,rev,3), 'FaceColor', 'b');
    plot(0:4, 0.5*ones(5,1), 'k:');
    axis([0 4 0 1]);  %axis square;  %grid on;
    set(gca,'Xtick',[],'Ytick',[],'box','off');
end
end
    
% save figure
set(gcf,'Color',[1,1,1]);
export_fig('MaxEntMaps.pdf')



end