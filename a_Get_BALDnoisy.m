% replace by a_Get_BALDscoreProp
function Dout=a_Get_BALDnoisy(D,DIM,yrM,pv,pw)
%%
% simulation: one "trial" only, uses GetEntropyGM_bruteforce (slow with brute force)
% function eval: looping through all, uses a_Model (fast with Huber)

% % Sumulation setup 1 begin
% type=1;
% xsf=20;
% ysf=20;
% im=15;
% % noise=0.01;
% 
% % see a_GenGPimage for input meaning
% scrsz = get(0,'ScreenSize');
% hf=figure('Name', 'Simulation Plot Window', 'position',[scrsz(3)/10 scrsz(4)/4 scrsz(3)/1.2 scrsz(4)/2],'menubar','none');
% 
% npix=770;  n=77;  nr=10;  px=linspace(1,npix,n);
% pvp=[1,log(0.01),0,0,0,0]; w=pvp(1); % parameters vector for perception
% pvd=[1,0,0.5];  % parameters vector for decision
% bm=1;
% lkm=-10;
% nm=0;
% mode=2;
% 
% % pvp=[w,lsigo,phit,bs,bs0];  % parameters vector for perception
% % pvd=[bd,lapse,prpa];  % parameters vector for decision
% 
% % npix=770;  n=77;  nr=25;  px=linspace(1,npix,n);
% % w=1;
% % pin=[npix/20, npix/20, npix/6, npix/30, log(0.5)];
% % prpa=0.3;
% % b1=0.3;
% % lapse=0.03;
% % Simulation setup 1 end

ntrial = D.Trials;
nrmax  = max(D.MaxRevealingTrial);

% function evaluation setup 1 begin
% setup input parameters
npix=770;
w=pv(9);
lsigo=pv(1);
phit=pv(2);
phis=pv(3);
phis0=pv(4);
phis_max=pv(5);
pvp=[w,lsigo,phit,phis,phis0,phis_max];  % parameters vector for perception
bd=pv(6);
lapse=pv(7);
prpa=pv(8);
pvd=[bd,lapse,prpa];  % parameters vector for decision
n=77;
bm=pw(1);  % softmax beta for motor choice
lkm=pw(2); % lapse rate for motor choice
nm=pw(3);  % SD of motor noise in pixels

% setup outputs
Dout=D;
RevPosX=NaN(nrmax,ntrial);
RevPosY=NaN(nrmax,ntrial);
Dout.MaxRevealingTrial=zeros(ntrial,1);

Dout.BALDxs = NaN(ntrial,nrmax);
Dout.PTILExs = NaN(ntrial,nrmax);
Dout.mla = NaN(ntrial,nrmax);
% function evaluation setup 1 end
 
% setup probe points
xx = linspace(1,npix,n);
yy = linspace(1,npix,n);
[x1,x2]=meshgrid(xx,yy);
z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
dxx = xx(2)-xx(1);

% setup gaussian motor noise conv map
if(mode==1 && nm>0)
    mu = [npix, npix]./2;
    covar = [nm^2, 0; 0, nm^2];
    F = mvnpdf([x1(:) x2(:)],mu,covar); F = reshape(F,length(x2),length(x1));
    ng=ceil(nm/dxx*5);
    Gmap = F(round(n/2-ng):round(n/2+ng),round(n/2-ng):round(n/2+ng));
    Gmap = Gmap/sum(sum(Gmap)); %normalize
end

% function evaluation setup 2 begin
for trial=1:ntrial
    yimage = DIM(trial).yimage;
    nr=D.MaxRevealingTrial(trial);
    yrin=yrM(trial,:)';
% function evaluation setup 2 end
    
% % Sumulation setup 2 begin
%     TT = load(strcat('im/Type',num2str(type),'X',num2str(xsf),'Y',num2str(ysf),'Num', num2str(im),'raw.mat'));
%     yimage = TT.yimage;
%     cumB=0;
%     yrin=0*yrM(3,:)';
% % Sumulation setup 2 end

    BALDxD = [0,0] + npix/2 + 1;  % in pixel
    count=1;
        
    while ( count<nr )    
        count=count+1;    
        %[lmla, lmlb, H_BALD, z] = a_Get_BALDmap_window(BALDxD, yimage, pin, w, prpa, b1, lapse, n, [1,npix], [1,npix]);
        [H_BALD, ~, lw, mu, var] = a_Model(BALDxD, yimage, pvp, yrin, pvd, n);
        
% % Sumulation setup 3 begin
%         [maxBs,ind]=max(H_BALD);
%         %ind=ceil(rand()*n*n);
%         hBs=H_BALD(ind);
%         
%         v=reshape(var,n*n,3); m=reshape(mu,n*n,3);
%         vx=v(ind,:); mx=m(ind,:); w=exp(lw)/sum(exp(lw));  w2=w(2:3)/sum(w(2:3));
%         B=GetEntropyGM_bruteforce(w,mx,vx)-w(1)*GetEntropyGM_bruteforce(1,mx(1),vx(1))-sum(w(2:3))*GetEntropyGM_bruteforce(w2,mx(2:3),vx(2:3)); 
% 
%         cumB=cumB+B;
%         
%         mla=exp(lw(1))/sum(exp(lw));  mlb= sum(exp(lw(2:3)))/sum(exp(lw));
%         entrog= -log(0.5) +mla*log(mla)+mlb*log(mlb);        
%         disp(sprintf('cum Bs=%f; Entropy gain=%f' ,cumB, entrog));
% % Sumulation setup 3 end

        matP = H_BALD.^(bm);  % softmax
        matP = matP/sum(matP);
        matP = (1-exp(lkm))*matP + exp(lkm)/n/n;  % lapse rate
        if(nm>0)
            Bmap = reshape(conv2(reshape(matP,n,n),Gmap,'same'),n*n,1); % convolve with motor noise
        else
            Bmap = matP;
        end
        intBmap = cumsum(Bmap);
        rseed=rand()*intBmap(end);
        ind = find(intBmap>rseed,1,'first');

        %hardmax
        %[~,ind]=max(H_BALD);            
        
        
%         %add motor noise directly
%         if(nm>0)
%             noise = nm/dxx;
%             indx=round(randn()*noise)*n;
%             while( ind+indx<1 || ind+indx>n*n )
%                 indx=round(randn()*noise)*n;
%             end
%             ind=ind+indx;
%             indy=round(randn()*noise);
%             while( floor((ind+indy)/n) ~= floor(ind/n) || ind+indy<1 || ind+indy>n*n )
%                 indy=round(randn()*noise);
%             end
%             ind=ind+indy;
%         end
        
        BALDxnew = z(ind,:); %z is in (x,y)            
        BALDxD = [BALDxD; BALDxnew];

% function evaluation setup 3 begin
        v=reshape(var,n*n,3); m=reshape(mu,n*n,3);
        vx=v(ind,:); mx=m(ind,:); w=exp(lw)/sum(exp(lw));  w2=w(2:3)/sum(w(2:3));
        B=GetEntropyGM_bruteforce(w,mx,vx)-w(1)*GetEntropyGM_bruteforce(1,mx(1),vx(1))-sum(w(2:3))*GetEntropyGM_bruteforce(w2,mx(2:3),vx(2:3)); 
        B=max(0,B);
        
        Dout.mla(trial,count-1)=exp(lw(1))/sum(exp(lw));
        Bxs=H_BALD(ind); 
        Dout.BALDxs(trial,count) = B;
        sBmap = sort(H_BALD);
        PTxs = find(sBmap<=Bxs,1,'last')/n/n*100;
        Dout.PTILExs(trial,count) = PTxs;
% function evaluation setup 3 end

% % Sumulation setup 4 begin
%         % reshape is to break rows into n subrows and put side by side
%         % z: [hor row 1 of image; hor row 2 of image; ...] --> so is H_BALD
%         
%         % plots
%         colormap(gray);
% 
%          subplot(1,2,1); imagesc(px, px, reshape(H_BALD,n,n)); hold on;         
%             plot( BALDxD(:,1), BALDxD(:,2), 'r.'); 
%             plot( BALDxnew(1), BALDxnew(2), 'ro','MarkerFaceColor','g', 'MarkerSize',5); hold off;
% 
%          subplot(1,2,2); imagesc(yimage); hold on;         
%             plot( BALDxD(:,1), BALDxD(:,2), 'r.'); 
%             plot( BALDxnew(1), BALDxnew(2), 'ro','MarkerFaceColor','g', 'MarkerSize',5);  hold off;
% 
%         % allows pressing q to exit the loop
%         if strcmp(get(hf,'currentcharacter'),'q')
%             break
%         end
% 
%         drawnow
%         pause
% % Sumulation setup 4 end
        
    end
    
% function evaluation setup 4 begin    
    Dout.MaxRevealingTrial(trial)=nr;
    RevPosX(1:nr,trial) = BALDxD(1:nr,1);
    RevPosY(1:nr,trial) = BALDxD(1:nr,2);
    disp(sprintf('trial=%d;',trial));
end

% from pixel to cm for display in experiment
% row[1]=10 y; column[1]= -10 x;
ss=20;  % screen size = [-10,10] cm
Dout.RevealPosX = (RevPosX'-1)/(npix-1)*ss-ss/2;  %Dout.RevealPosX(Dout.RevealPosX>9.9)=9.9;  Dout.RevealPosX(Dout.RevealPosX<-9.9)=-9.9;
Dout.RevealPosY = (1-RevPosY')/(npix-1)*ss+ss/2;  %Dout.RevealPosY(Dout.RevealPosY>9.9)=9.9;  Dout.RevealPosY(Dout.RevealPosY<-9.9)=-9.9;
% function evaluation setup 4 end    

end
