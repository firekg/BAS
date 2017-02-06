clear;
load for_mate;

%%

offdiag=@(A) A(~eye(size(A)));

%% making sure we use the same number of trials per image type

imagetypes=unique(D.image_type);
Nimtyp=length(imagetypes);

Ntrials=zeros(1,Nimtyp);

for j=1:Nimtyp
    Ntrials(j)=length(unique(D.trial(D.image_type==imagetypes(j))));
end

Ntrials=min(Ntrials);

% ttix=false(size(D.trial));
% for j=1:Nimtyp
%     iix=(D.image_type==imagetypes(j));
%     trials=unique(D.trial(iix));
%     tix=randperm(length(trials));
%     ttix=ttix | (iix & ismember(D.trial,trials(tix(1:Ntrials))));
% end
% 
% fn=fieldnames(D);
% for i=1:length(fn)
%     X=D.(fn{i});
%     D.(fn{i})=X(ttix);
% end

%% maps using all the data

xgrid=linspace(-10,10,11);
ygrid=linspace(-10,10,11);

M=createmap(D.x,D.y,xgrid,ygrid);

imagetypes=unique(D.image_type);
Nimtyp=length(imagetypes);

[Mit,Mitc]=deal(zeros(length(xgrid),length(ygrid),Nimtyp));
for j=1:Nimtyp
    ix=(D.image_type==imagetypes(j));
    Mit(:,:,j)=createmap(D.x(ix),D.y(ix),xgrid,ygrid);
    Mitc(:,:,j)=Mit(:,:,j)-M;
end

%         M1(:,:,j)=createmap(D.x(ix1),D.y(ix1),xgrid,ygrid);
%         M2(:,:,j)=createmap(D.x(ix2),D.y(ix2),xgrid,ygrid);
% 
%         M1c(:,:,j)=M1(:,:,j)-M;
%         M2c(:,:,j)=M2(:,:,j)-M;

%%

figure(1);
clf;

clims_map=[0 0.4];
clims_diffmap=[-0.1 0.1];

subplot(3,3,2);
imagesc(M,clims_map);
set(gca,'FontSize',16,'XTick',[],'YTick',[],'XDir','normal','YDir','normal');
axis square;
colorbar;

for j=1:Nimtyp    
    subplot(3,3,3+j);
    imagesc(Mit(:,:,j),clims_map);
    set(gca,'FontSize',16,'XTick',[],'YTick',[],'XDir','normal','YDir','normal');
    axis square;
    colorbar;
    
    subplot(3,3,6+j);
    imagesc(Mitc(:,:,j),clims_diffmap/8);
    set(gca,'FontSize',16,'XTick',[],'YTick',[],'XDir','normal','YDir','normal');
    axis square;
    colorbar;

end

%%

Nsplits=10;

[M1,M2,M1c,M2c]=deal(zeros(length(xgrid),length(ygrid),Nimtyp));
[r,rc]=deal(nan(Nimtyp,Nimtyp,Nsplits));

splitsize=@(N) floor(Ntrials/2);

for i=1:Nsplits
    
    for j=1:Nimtyp
        
        disp([i j]);

        iix=(D.image_type==imagetypes(j));
        
        %% splitting by trials
        
        trials=unique(D.trial(iix));
        trials=trials(randperm(length(trials)));
        
        Ntrials=splitsize(length(trials));
        
        ix1=(iix & ismember(D.trial,trials(1:Ntrials)));
        ix2=(iix & ismember(D.trial,trials(Ntrials+1:2*Ntrials)));
                
        %% splitting by revealings
        
%         revs=find(iix);
%         revs=revs(randperm(length(revs)));
%         
%         Nrevs=splitsize(length(revs));
%         
%         ix1=revs(1:Nrevs);
%         ix2=revs(Nrevs+1:2*Nrevs);
        
        %% computing maps and correlations
        M1(:,:,j)=createmap(D.x(ix1),D.y(ix1),xgrid,ygrid);
        M2(:,:,j)=createmap(D.x(ix2),D.y(ix2),xgrid,ygrid);

        M1c(:,:,j)=M1(:,:,j)-M;
        M2c(:,:,j)=M2(:,:,j)-M;

        r(j,j,i)=corr2(M1(:,:,j),M2(:,:,j));
        rc(j,j,i)=corr2(M1c(:,:,j),M2c(:,:,j));
        if j>1
            for k=1:j-1
                r(j,k,i)=corr2(M1(:,:,j),M2(:,:,k));
                r(k,j,i)=corr2(M1(:,:,k),M2(:,:,j));
                rc(j,k,i)=corr2(M1c(:,:,j),M2c(:,:,k));
                rc(k,j,i)=corr2(M1c(:,:,k),M2c(:,:,j));
            end
        end

        %%
        
        figure(2);
        
        subplot(2,3,j);
        imagesc(M1(:,:,j),clims_map);
        set(gca,'FontSize',16,'XTick',[],'YTick',[],'XDir','normal','YDir','normal');
        axis square;
        colorbar;
        
        subplot(2,3,3+j);
        imagesc(M2(:,:,j),clims_map);
        set(gca,'FontSize',16,'XTick',[],'YTick',[],'XDir','normal','YDir','normal');
        axis square;
        colorbar;
        
        drawnow;
                        
        figure(3);
        
        subplot(2,3,j);
        imagesc(M1c(:,:,j),clims_diffmap);
        set(gca,'FontSize',16,'XTick',[],'YTick',[],'XDir','normal','YDir','normal');
        axis square;
        colorbar;
        
        subplot(2,3,3+j);
        imagesc(M2c(:,:,j),clims_diffmap);
        set(gca,'FontSize',16,'XTick',[],'YTick',[],'XDir','normal','YDir','normal');
        axis square;
        colorbar;
        
        drawnow;
        
    end
end

%%

figure(4);

for j=1:Nimtyp
    for k=1:Nimtyp
        subplot(Nimtyp,Nimtyp,(j-1)*Nimtyp+k);
        histogram(r(j,k,:));        
    end
end

figure(5);

for j=1:Nimtyp
    for k=1:Nimtyp
        subplot(Nimtyp,Nimtyp,(j-1)*Nimtyp+k);
        histogram(rc(j,k,:));        
    end
end

%%

rf=atanh(r);

figure(6);

for j=1:Nimtyp
    for k=1:Nimtyp
        subplot(Nimtyp,Nimtyp,(j-1)*Nimtyp+k);
        histogram(rf(j,k,:));        
    end
end

%%

figure(10); 
clf; 

hold on; 
plot([0 0],[-8 6],'-k',[-2 7],[0 0],'-k');
ix=(D.image_type==1)&(D.rev==1); %&(rand(size(D.rev))<0.5); 
p1=plot(D.x(ix),D.y(ix),'.r'); 
ix=(D.image_type==2)&(D.rev==1); 
p2=plot(D.x(ix),D.y(ix),'.b'); 
ix=(D.image_type==3)&(D.rev==1); 
p3=plot(D.x(ix),D.y(ix),'.g'); hold off;
set(gca,'FontSize',16,'box','on','LineWidth',2,'XLim',[-2 7],'YLim',[-8 6]); 
xlabel('fixation x'); 
ylabel('fixation y');
legend([p1 p2 p3],{'image type 1','image type 2','image type 3'}); 
legend boxoff;

figure(11); 
clf; 

subplot(3,1,1);

XX=logspace(-3,1,30); %[0.001:0.001:0.01 0.02:0.01:0.1 0.2:0.1:1 2:10]; 

ix=(D.image_type==1)&(D.rev==1);
[N1,X1]=hist(sqrt(D.x(ix).^2+D.y(ix).^2),XX); 
ix=(D.image_type==2)&(D.rev==1); 
[N2,X2]=hist(sqrt(D.x(ix).^2+D.y(ix).^2),XX); 
ix=(D.image_type==3)&(D.rev==1); 
[N3,X3]=hist(sqrt(D.x(ix).^2+D.y(ix).^2),XX); 

plot(X1,N1/sum(N1),'-r',X2,N2/sum(N2),'-b',X3,N3/sum(N3),'-g','LineWidth',2); 
set(gca,'XScale','log','FontSize',16,'box','off','LineWidth',2); 
xlabel('distance of first fixation from origin'); 
ylabel('#trials');
legend({'image type 1','image type 2','image type 3'}); 
legend boxoff;

subplot(3,1,2);

XX=linspace(-pi,pi,30);

ix=(D.image_type==1)&(D.rev==1);
[N1,X1]=hist(atan2(D.x(ix),D.y(ix)),XX); 
ix=(D.image_type==2)&(D.rev==1); 
[N2,X2]=hist(atan2(D.x(ix),D.y(ix)),XX); 
ix=(D.image_type==3)&(D.rev==1); 
[N3,X3]=hist(atan2(D.x(ix),D.y(ix)),XX); 

plot(X1,N1/sum(N1),'-r',X2,N2/sum(N2),'-b',X3,N3/sum(N3),'-g','LineWidth',2); 
set(gca,'FontSize',16,'box','off','LineWidth',2); 
xlabel('angle of first fixation from origin'); 
ylabel('#trials');

subplot(3,1,3);

XX=linspace(-pi,pi,30);
XXe=[-inf XX(1:end-1)+diff(XX)/2 inf];

ix=(D.image_type==1)&(D.rev==1);
[N1,~,B]=histcounts(atan2(D.x(ix),D.y(ix)),XXe); 
D1=sqrt(accumarray(B',D.x(ix)).^2+accumarray(B',D.y(ix)).^2);
ix=(D.image_type==2)&(D.rev==1); 
[N2,~,B]=histcounts(atan2(D.x(ix),D.y(ix)),XXe); 
D2=sqrt(accumarray(B',D.x(ix)).^2+accumarray(B',D.y(ix)).^2);
ix=(D.image_type==3)&(D.rev==1); 
[N3,~,B]=histcounts(atan2(D.x(ix),D.y(ix)),XXe); 
D3=sqrt(accumarray(B',D.x(ix)).^2+accumarray(B',D.y(ix)).^2);

plot(XX,D1./N1','-r',XX,D2./N2','-b',XX,D3./N3','-g','LineWidth',2); 
% plot(XX,D1','-r',XX,D2','-b',XX,D3','-g','LineWidth',2); 
set(gca,'FontSize',16,'box','off','LineWidth',2); 
xlabel('angle of first fixation from origin'); 
ylabel('average distance');
% ylabel('total distance');
