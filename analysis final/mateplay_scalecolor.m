clear all;
load BAS_data_Mate;

%%

figure(1); 

clf; 

graymin=0;
graymax=1;
grayres=256;
graywarp=1/4;

gf=@(graylevel) graylevel.^graywarp;

graylevel=linspace(graymin,graymax,grayres); 
cm=gf(graylevel)'*ones(1,3);  % scale color map


ww=[1/10 1/5 1/4 1/3 1/2 1];

[ix,iy,p,p0]=deal(zeros(size(x)));

set(gcf,'ColorMap',cm,'PaperPositionMode','auto');

for i=1:size(data,3); 
    subplot(1,9,i); 
    
    D=data(:,:,i);
    D=D(:);
    D2=interp1(linspace(0,1,grayres),(1:grayres),D/max(D),'nearest');
    
    D2=reshape(D2,size(data(:,:,i)));

    image(px,px,D2,'CDataMapping','direct');
    hold on;
    plot(x(i),y(i),'o','MarkerSize',5,'MarkerEdgeColor',[1 1 0],'LineWidth',2);
    hold off;
        
    axis square; 
    set(gca,'Visible','off'); 
    
    ix(i)=interp1(px,(1:size(data,2)),x(i),'nearest');
    iy(i)=interp1(px,(1:size(data,1)),y(i),'nearest');
    d=data(iy(i),ix(i),i);
    p(i)=sum(D<=d)/numel(D)*100;   
    
    h=gca;
    pos=get(h,'Position');
    
end

b=colorbar;
set(b,'YTick',[1 grayres],'YTickLabel',{'0','max'});
set(h,'Position',pos);

saveas(1,'mateplot_scalecolor','epsc2');
disp([90 91 77 72 75 85 96 94 98; p]);
