%% Plot: sample images, PERF, INFO (nominally Fig. 1)
load('PERF.mat');
load('INFOf.mat');

figure(1);
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
set(gcf,'Position',[1,2,20,20]);
set(gcf, 'PaperSize', [20,20]);
set(gcf, 'PaperPosition', [1,2,20,20]);

%plot sample iamges
pd=0.005; %distance between plots
ps=0.07; %plot size
lm=0.05; %left margin
ls1=lm+2*pd; %1st column left start position
ts1=0.85; %1st row top start position
%ii for columns; jj for rows (bad convention...)
for ii=1:3
for jj=1:5
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]
    if(ii==1)      s=['im\Type1X20Y20Num',num2str(jj)];
    elseif(ii==2)  s=['im\Type2X6Y30Num',num2str(jj)];
    elseif(ii==3)  s=['im\Type2X30Y6Num',num2str(jj)];
    end
    img=imread(s,'bmp');
    image(img);  set(gca,'Xtick',[],'Ytick',[]);
end
end

%plot PERF
pd=0.04; %distance between plots
pws=0.2; %plot width
phs=0.16; %plot height
lm=0.00; %left margin
ls1=lm+2*pd; %1st column left start position
ts2=0.25; %1st row top start position
x=5:5:25;
% PERF has 4 layers: revealing number x5(1:5:25); parti x4; condiitons x3 (AL,PLR,PLB); moments x2 (mean,SD)
for ii=1:4
    ls=ls1+(ii-1)*pws+(ii-1)*pd;
    ts=ts2;
    subplot('Position',[ls, ts, pws-pd/2, phs-pd/2]);
    plot([0,30],[0.5,0.5],':','Color',[0.5,0.5,0.5]);  hold on;  plot([0,30],[1,1],':','Color',[0.5,0.5,0.5]);
    errorbar(x,PERF(:,ii,3,1),PERF(:,ii,3,2),'k^-','MarkerFaceColor',[0,0,0],'MarkerSize',3);
    errorbar(x,PERF(:,ii,2,1),PERF(:,ii,2,2),'bs-','MarkerFaceColor',[0,0,1],'MarkerSize',3);
    errorbar(x,PERF(:,ii,1,1),PERF(:,ii,1,2),'ro-','MarkerFaceColor',[1,0,0],'MarkerSize',3);    
%     if(ii==1)
%         set(gca,'Ytick',0.3:0.1:1,'YTickLabel',{'','','0.5','','','','','1'}, 'TickLength',[0.01,0]);    
%         ylabel('P(correct)');  % ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [3, 0, 0]);
%     else
%         set(gca,'Ytick',0.3:0.1:1,'YTickLabel',{'','','','','','','',''}, 'TickLength',[0.01,0]);    
%     end
    ylabel('P(correct)');  % ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [3, 0, 0]);  set(gca,'FontSize',10);
    set(gca,'Ytick',0.3:0.1:1,'YTickLabel',{'','','0.5','','','','','1'}, 'TickLength',[0.01,0]);    
    xlabel('Revealing number');
    set(gca,'Xtick',0:5:25,'XTickLabel',{'0','5','10','15','20','25'}, 'TickLength',[0.01,0]);
    axis([0 27.5 0.3 1.05]);  box off;
    title(['P',num2str(ii)]);
    set(gca,'FontSize',10);
end

%plot INFO
pd=0.04; %distance between plots
pws=0.2; %plot width
phs=0.16; %plot height
lm=0.00; %left margin
ls1=lm+2*pd; %1st column left start position
ts2=0.05; %1st row top start position
x=1:25;
% INFO has 4 layers: revealing number x25 (1:25); parti x3; condiitons x4 (AL,PLR,BAS,nBAS); moments x2 (mean,SD)
for ii=1:4
    ls=ls1+(ii-1)*pws+(ii-1)*pd;
    ts=ts2;
    subplot('Position',[ls, ts, pws-pd/2, phs-pd/2]);
    %plot([0,30],[1,1],':','Color',[0.5,0.5,0.5]);   %line on 1
    %BAS->nBAS->PLR->AL
    cs=[0.7,0.7,0.7];  a_ErrorShade(x,INFO(:,ii,3,1)-INFO(:,ii,3,2),INFO(:,ii,3,1)+INFO(:,ii,3,2),cs,1);  hold on;
    cs=[0.7,0.7,1.0];  a_ErrorShade(x,INFO(:,ii,2,1)-INFO(:,ii,2,2),INFO(:,ii,2,1)+INFO(:,ii,2,2),cs,1);
    cs=[1.0,0.7,0.7];  a_ErrorShade(x,INFO(:,ii,1,1)-INFO(:,ii,1,2),INFO(:,ii,1,1)+INFO(:,ii,1,2),cs,1);
    cs=[0.7,1.0,0.7];  a_ErrorShade(x,INFO(:,ii,4,1)-INFO(:,ii,4,2),INFO(:,ii,4,1)+INFO(:,ii,4,2),cs,1);
    cs=[0.2,0.5,0.2];  a_ErrorShade(x,INFO(:,ii,5,1)-INFO(:,ii,5,2),INFO(:,ii,5,1)+INFO(:,ii,5,2),cs,1);
    cl=[0.0,0.0,0.0];  plot(x,INFO(:,ii,3,1),'-','Color',cl, 'MarkerFaceColor',cl);
    cl=[0.0,0.0,1.0];  plot(x,INFO(:,ii,2,1),'-','Color',cl, 'MarkerFaceColor',cl);    
    cl=[1.0,0.0,0.0];  plot(x,INFO(:,ii,1,1),'-','Color',cl, 'MarkerFaceColor',cl);
    cl=[0.0,0.8,0.0];  plot(x,INFO(:,ii,4,1),'-','Color',cl, 'MarkerFaceColor',cl);
    cl=[0.2,0.4,0.2];  plot(x,INFO(:,ii,5,1),'-','Color',cl, 'MarkerFaceColor',cl);
%     if(ii==1)
%         set(gca,'Ytick',0:0.1:1,'YTickLabel',{'0','','','','','0.5','','','','','1'}, 'TickLength',[0.01,0]);    
%         ylabel('Information gain');  % ylabh = get(gca,'YLabel'); set(ylabh,'Position',get(ylabh,'Position') + [3, 0, 0]);  set(gca,'FontSize',10);
%     else
%         set(gca,'Ytick',0:0.1:1,'YTickLabel',{'','','','','','','','','','',''}, 'TickLength',[0.01,0]);    
%     end
    ylabel('Information (bit)');
    set(gca,'Ytick',0:0.1:1,'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','','','','','1'}, 'TickLength',[0.01,0]);    
    %set(gca,'Ytick',0:0.1:1,'YTickLabel',{'0','','','','','0.5','','','','','1'}, 'TickLength',[0.01,0]);    
    xlabel('Revealing number');
    set(gca,'Xtick',0:5:25,'XTickLabel',{'0','5','10','15','20','25'}, 'TickLength',[0.01,0]);
    axis([0 27.5 0 0.5]);  box off;  
    %axis([0 27.5 0 1.05]);  box off;  
    set(gca,'FontSize',10);
end
%leg=legend('BAS', 'biased nBAS', 'free-scan', 'random', 'Location','SouthEast'); 
%set(leg,'Color','w');

% % save figure
% set(gcf,'Color',[1,1,1]);
% export_fig('PerfInfoSigo1.pdf')

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

% differences participants
for ii=1:3
for jj=1:4
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]
    img=revMapColor(drevMap(:,:,ii,jj,1));
    image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
    %imagesc(drevMap(:,:,ii,jj,1)); colormap(hot);  freezeColors;  axis off;    
end
end
% mean participants
ii=4;
for jj=1:4
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]    
    img=revMapColor(revMap(:,:,ii,jj,1));
    image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
    %imagesc(revMap(:,:,ii,jj,1)); colormap(hot);  freezeColors;  axis off;    
end

% differences BAS
jj=5;
for ii=1:3
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]    
    img=revMapColor(drevMap(:,:,ii,4,2));
    image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
    %imagesc(drevMap(:,:,ii,4,2)); colormap(hot);  freezeColors;  axis off;    
end
% mean BAS
ii=4;
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]    
    img=revMapColor(revMap(:,:,ii,4,2));
    image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
    %imagesc(revMap(:,:,ii,4,2)); colormap(hot);  freezeColors;  axis off;    

% differences noisy BAS
jj=6;
for ii=1:3
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]    
    img=revMapColor(drevMap(:,:,ii,4,3));
    image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
    %imagesc(drevMap(:,:,ii,4,3)); colormap(hot);  freezeColors;  axis off;    
end
% mean noisy BAS
ii=4;
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]    
    img=revMapColor(revMap(:,:,ii,4,3));
    image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
    %imagesc(revMap(:,:,ii,4,3)); colormap(hot);  freezeColors;  axis off;    

    
%legned bar differece
jj=7; ii=1;
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts+4*(ps-pd/2)/5, ps-pd/2, (ps-pd/2)/5]); %[jj,ii]    
    leg=revMapColor(linspace(-1,1,nleg));
    image(leg);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);

%legned bar mean
jj=7; ii=4;
    ls=ls1+(ii-1)*ps+(ii-1)*pd;
    ts=ts1-(jj-1)*ps-(jj-1)*pd;
    subplot('Position',[ls, ts+4*(ps-pd/2)/5, ps-pd/2, (ps-pd/2)/5]); %[jj,ii]    
    leg=revMapColor(linspace(0,1,nleg));
    image(leg);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);

% % save figure
% set(gcf,'Color',[1,1,1]);
% export_fig('revMapCsigo1.pdf')

%% Plot SMboot--correlation (nominally Fig. 3)
figure(3);
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
set(gcf,'Position',[1,2,20,20]);
set(gcf, 'PaperSize', [20,20]);
set(gcf, 'PaperPosition', [1,2,20,20]);

load('CORR4.mat');

%plot CORR
pd=0.04; %distance between plots
pws=0.18; %plot width
phs=0.15; %plot height
lm=0.08; %left margin
ls1=lm; %1st column left start position
ts1=0.70; %1st row top start position
x=1:25;
% CORR has 5 layers:
% revealing number 25x: 1:25
% moments x3: (mean, upper 95%, lower 95%)
% pattern (mis)match 2x: same, different
% conditions 3x: self-self, self-BAS, self-pnBAS
% participants 4x: SC, EP, BZ, AVG
for ii=4
for jj=1:3
    ls=ls1;
    ts=ts1-(jj-1)*phs-(jj-1)*pd;
    subplot('Position',[ls, ts, pws-pd/2, phs-pd/2]);    
    plot([0,30],[0,0],'-','Color',[0.0,0.0,0.0]); hold on;   %line on 0
    cl=[1.0,0.5,0.0];  
    plot(x,CORR(:,1,1,jj,ii),'-','Color',cl, 'MarkerFaceColor',cl);
    plot(x,CORR(:,2,1,jj,ii),':','Color',cl, 'MarkerFaceColor',cl);
    plot(x,CORR(:,3,1,jj,ii),':','Color',cl, 'MarkerFaceColor',cl);
    cl=[0.5,0.0,1.0];  
    plot(x,CORR(:,1,2,jj,ii),'-','Color',cl, 'MarkerFaceColor',cl);
    plot(x,CORR(:,2,2,jj,ii),':','Color',cl, 'MarkerFaceColor',cl);
    plot(x,CORR(:,3,2,jj,ii),':','Color',cl, 'MarkerFaceColor',cl);
    ylabel('Correlation');
    set(gca,'Ytick',-1:0.5:1,'YTickLabel',{'-1','','0','','1'}, 'TickLength',[0.01,0]);    
    xlabel('Revealing number');
    set(gca,'Xtick',0:5:25,'XTickLabel',{'0','5','10','15','20','25'}, 'TickLength',[0.01,0]);
    axis([0 25 -1.05 1.05]);  box off;  
    set(gca,'FontSize',10);
end
end

% barplots of the 25th CORR
pwd=0.02; %distance width between plots
phd=0.04; %distance height between plots
pws=0.08; %plot width
phs=0.15; %plot height
lm=0.08; %left margin
ls1=0.3; %1st column left start position
ts1=0.70; %1st row top start position
for ii=1:4
for jj=1:3
    ls=ls1+(ii-1)*pws+(ii-1)*pwd;
    ts=ts1-(jj-1)*phs-(jj-1)*phd;
    subplot('Position',[ls, ts, pws-pwd/2, phs-phd/2]);
    %same
    cb=[1.0,0.5,0.0]; 
    bar(1, CORR(25,1,1,jj,ii), 'FaceColor', cb); hold on;    
    plot([1,1],[CORR(25,2,1,jj,ii),CORR(25,3,1,jj,ii)],'k-','LineWidth',1);
    %diff
    cb=[0.5,0.0,1.0];
    bar(2, CORR(25,1,2,jj,ii), 'FaceColor', cb);
    plot([2,2],[CORR(25,2,2,jj,ii),CORR(25,3,2,jj,ii)],'k-','LineWidth',1);
    axis([0 3 -1 1]);  %axis square;  %grid on;
    set(gca,'Ytick',-1:0.5:1,'YTickLabel',{'-1','','0','','1'}, 'TickLength',[0.01,0]);    
    set(gca,'Xtick',[],'box','off');
end
end
% save figure
set(gcf,'Color',[1,1,1]);
export_fig('Corr4.pdf')

%% Plot SMboot--correlation for each participants (nominally Fig. 4)
figure(4);
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
set(gcf,'Position',[1,2,20,20]);
set(gcf, 'PaperSize', [20,20]);
set(gcf, 'PaperPosition', [1,2,20,20]);

%plot CORR
load('CORR4.mat');
pd=0.04; %distance between plots
pws=0.18; %plot width
phs=0.15; %plot height
lm=0.08; %left margin
ls1=lm; %1st column left start position
ts1=0.75; %1st row top start position
x=1:25;
for ii=1:4
for jj=1:3
    ls=ls1+(ii-1)*pws+(ii-1)*pd;
    ts=ts1-(jj-1)*phs-(jj-1)*pd;
    subplot('Position',[ls, ts, pws-pd/2, phs-pd/2]);
    plot([0,30],[0,0],'-','Color',[0.0,0.0,0.0]); hold on;   %line on 0
    cl=[1.0,0.5,0.0];  
    plot(x,CORR(:,1,1,jj,ii),'-','Color',cl, 'MarkerFaceColor',cl);
    plot(x,CORR(:,2,1,jj,ii),':','Color',cl, 'MarkerFaceColor',cl);
    plot(x,CORR(:,3,1,jj,ii),':','Color',cl, 'MarkerFaceColor',cl);
    cl=[0.5,0.0,1.0];  
    plot(x,CORR(:,1,2,jj,ii),'-','Color',cl, 'MarkerFaceColor',cl);
    plot(x,CORR(:,2,2,jj,ii),':','Color',cl, 'MarkerFaceColor',cl);
    plot(x,CORR(:,3,2,jj,ii),':','Color',cl, 'MarkerFaceColor',cl);
    ylabel('Correlation');
    set(gca,'Ytick',-1:0.5:1,'YTickLabel',{'-1','','0','','1'}, 'TickLength',[0.01,0]);    
    xlabel('Revealing number');
    set(gca,'Xtick',0:5:25,'XTickLabel',{'0','5','10','15','20','25'}, 'TickLength',[0.01,0]);
    axis([0 25 -1.0 1.0]);  box off;  
    set(gca,'FontSize',10);
end
end
% save figure
set(gcf,'Color',[1,1,1]);
export_fig('CorrParticipants4.pdf')


%% Plot revMap with revealing number progression (nominally Fig. 5)
figure(5);
set(gcf,'Units','centimeters');
set(gcf,'Toolbar','None');
% set(gcf,'MenuBar','None');
set(gcf,'Position',[1,2,20,20]);
set(gcf, 'PaperSize', [20,20]);
set(gcf, 'PaperPosition', [1,2,20,20]);

pd=0.001; %distance between plots
ps=0.038; %plot size
lm=0.04; %left margin
ls1=lm; %1st column left start position
ts1=0.95; %1st row top start position

%define image color functions
rf = @(M) min(1, (1+M).^2.5);
gf = @(M) (1 - abs(M)).^2.5;
bf = @(M) min(1, (1-M).^2.5);

img=zeros(770/2,770/2,3); % initialize image
nleg=200;  %number of entry in legend
leg=zeros(1,nleg,3); % initialized legend

for parti=1:3
    % differences participants
    for ii=1:3
    for jj=1:25
        ls=ls1+(ii-1)*ps+(ii-1)*pd + (parti-1)*0.2;
        ts=ts1-(jj-1)*ps-(jj-1)*pd;
        subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]
        img(:,:,1) = feval(rf,RdrevMap(jj,:,:,ii,parti,1)); %r    
        img(:,:,2) = feval(gf,RdrevMap(jj,:,:,ii,parti,1)); %g    
        img(:,:,3) = feval(bf,RdrevMap(jj,:,:,ii,parti,1)); %b    
        image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
        %imagesc(drevMap(:,:,ii,jj,1)); colormap(hot);  freezeColors;  axis off;    
    end
    end
    % mean participants
    ii=4;
    for jj=1:25
        ls=ls1+(ii-1)*ps+(ii-1)*pd  + (parti-1)*0.2;
        ts=ts1-(jj-1)*ps-(jj-1)*pd;
        subplot('Position',[ls, ts, ps-pd/2, ps-pd/2]); %[jj,ii]    
        img(:,:,1) = feval(rf,RrevMap(jj,:,:,ii,parti,1)); %r    
        img(:,:,2) = feval(gf,RrevMap(jj,:,:,ii,parti,1)); %g    
        img(:,:,3) = feval(bf,RrevMap(jj,:,:,ii,parti,1)); %b    
        image(img);  freezeColors;  set(gca,'Xtick',[],'Ytick',[]);
        %imagesc(revMap(:,:,ii,jj,1)); colormap(hot);  freezeColors;  axis off;    
    end
end
