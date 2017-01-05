function plot_rev_maps_order(RrevMap, RdrevMap)
%% Plot drevMap, revMap of avg participant with revealing order: still messy
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
set(gcf,'Position',[1,2,22,22]);
set(gcf, 'PaperSize', [22,22]);
set(gcf, 'PaperPosition', [1,2,22,22]);

pd=0.001; %distance between plots
ps=0.03; %plot size
lm=0.01; %left margin
ls1=lm+2*pd; %1st column left start position
ts1=0.95; %1st row top start position

% %define image color functions
% rf = @(M) min(1, 0.55+(1+M).^7.5);
% gf = @(M) (1 - abs(M)).^1.0;
% bf = @(M) min(1, (1-M).^7.5);

img = zeros(770,770,3); % initialize image
nleg = 200;  %number of entry in legend
leg = zeros(1,nleg,3); % initialized legend

parti = 4; %avg parti
revVec = 1:25;
cond = 2; %BAS
incr = 10;

% difference maps
for col=1:3
for row=1:25
    rev = revVec(row);
    ls=ls1+(col-1)*ps+(col-1)*pd;
    ts=ts1-(row-1)*ps-(row-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]
    temp = RdrevMap(:,:,col,parti,cond,rev);
    img = revMapColor(subSampImg(temp, incr));
    image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
    %imagesc(drevMap(:,:,ii,jj,1)); colormap(hot);  freezeColors;  axis off;    
end
end
% avg maps
for col = 1:4;
for row = 1:25
    rev = revVec(row);
    colsh = col + 3;
    ls=ls1+(colsh-1)*ps+(colsh-1)*pd;
    ts=ts1-(row-1)*ps-(row-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]
    temp = RrevMap(:,:,col,parti,cond,rev);
    img = revMapColor(subSampImg(temp, incr));
    img(:,:,1) = img(:,:,2);
    img(:,:,3) = img(:,:,2); %make it black and white
    image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
    %imagesc(revMap(:,:,ii,jj,1)); colormap(hot);  freezeColors;  axis off;    
end
end

% %legned bar differece
% row=8; col=1;
%     ls=ls1+(col-1)*ps+(col-1)*pd;
%     ts=ts1-(row-1)*ps-(row-1)*pd;
%     subplot('Position',[ls, ts+4*(ps-pd/2)/5, ps-pd/2, (ps-pd/2)/5]); %[jj,ii]    
%     leg=revMapColor(linspace(-1,1,nleg));
%     image(leg);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
% 
% %legned bar mean
% row=8; col=4;
%     ls=ls1+(col-1)*ps+(col-1)*pd;
%     ts=ts1-(row-1)*ps-(row-1)*pd;
%     subplot('Position',[ls, ts+4*(ps-pd/2)/5, ps-pd/2, (ps-pd/2)/5]); %[jj,ii]    
%     leg=revMapColor(linspace(0,1,nleg));
%     leg(:,:,1) = leg(:,:,2);
%     leg(:,:,3) = leg(:,:,2); %make it black and white
%     image(leg);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);

set(gcf,'Color',[1,1,1]);

% save figure
% export_fig('revMap_order.pdf')

end
%%
function subSampImg = subSampImg(img, incr)
    subSampImg = img(1:incr:end,:);
    subSampImg = subSampImg(:,1:incr:end);
end