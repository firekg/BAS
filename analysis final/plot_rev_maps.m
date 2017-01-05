function plot_rev_maps(revMap, drevMap)
%% Plot drevMap, revMap (nominally Fig. 2)
figure(2);
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
set(gcf,'Position',[1,2,20,20]);
set(gcf, 'PaperSize', [20,20]);
set(gcf, 'PaperPosition', [1,2,20,20]);

pd=0.005; %distance between plots
ps=0.1; %plot size
lm=0.05; %left margin
ls1=lm+2*pd; %1st column left start position
ts1=0.85; %1st row top start position

% %define image color functions
% rf = @(M) min(1, 0.55+(1+M).^7.5);
% gf = @(M) (1 - abs(M)).^1.0;
% bf = @(M) min(1, (1-M).^7.5);

img=zeros(770,770,3); % initialize image
nleg=200;  %number of entry in legend
leg=zeros(1,nleg,3); % initialized legend

nparti = 4;
showRevn = 25;

% differences of participants
for ii=1:3
for jj=1:nparti
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]
    img=revMapColor(drevMap(:,:,ii,jj,1,showRevn));
    image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
    %imagesc(drevMap(:,:,ii,jj,1)); colormap(hot);  freezeColors;  axis off;    
end
end
% mean of participants
ii=4;
for jj=1:nparti
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]    
    img=revMapColor(revMap(:,:,ii,jj,1,showRevn));
    img(:,:,1) = img(:,:,2);
    img(:,:,3) = img(:,:,2); %make it black and white
    image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
    %imagesc(revMap(:,:,ii,jj,1)); colormap(hot);  freezeColors;  axis off;    
end

% % differences of nnBAS
% nnBAS = 3;
% jj=5;
% for ii=1:3
%     ls=ls1+(ii-1)*ps+(ii-1)*pd;
%     ts=ts1-(jj-1)*ps-(jj-1)*pd;
%     subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]    
%     img=revMapColor(drevMap(:,:,ii,4,nnBAS));
%     image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
%     %imagesc(drevMap(:,:,ii,4,2)); colormap(hot);  freezeColors;  axis off;    
% end
% % mean of nn BAS
% ii=4;
%     ls=ls1+(ii-1)*ps+(ii-1)*pd;
%     ts=ts1-(jj-1)*ps-(jj-1)*pd;
%     subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]    
%     img=revMapColor(revMap(:,:,ii,4,nnBAS));
%     img(:,:,1) = img(:,:,2);
%     img(:,:,3) = img(:,:,2); %make it black and white
%     image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
%     %imagesc(revMap(:,:,ii,4,2)); colormap(hot);  freezeColors;  axis off;    

% differences of noisy BAS
mnBAS = 2;
jj = nparti+1;
for ii=1:3
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]    
    img=revMapColor(drevMap(:,:,ii,4,mnBAS,showRevn));
    image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
    %imagesc(drevMap(:,:,ii,4,3)); colormap(hot);  freezeColors;  axis off;    
end
% mean of noisy BAS
ii=4;
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]    
    img=revMapColor(revMap(:,:,ii,4,mnBAS,showRevn));
    img(:,:,1) = img(:,:,2);
    img(:,:,3) = img(:,:,2); %make it black and white
    image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
    %imagesc(revMap(:,:,ii,4,3)); colormap(hot);  freezeColors;  axis off;    

%legned bar differece
jj = nparti+2;
ii=1;
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts+4*(ps-pd/2)/5, ps-pd/2, (ps-pd/2)/5]); %[jj,ii]    
    leg=revMapColor(linspace(-1,1,nleg));
    image(leg);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);

%legned bar mean
jj = nparti+2;
ii=4;
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts+4*(ps-pd/2)/5, ps-pd/2, (ps-pd/2)/5]); %[jj,ii]    
    leg=revMapColor(linspace(0,1,nleg));
    leg(:,:,1) = leg(:,:,2);
    leg(:,:,3) = leg(:,:,2); %make it black and white
    image(leg);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);

set(gcf,'Color',[1,1,1]);

% save figure
% export_fig('revMap_newNoise.pdf')

end