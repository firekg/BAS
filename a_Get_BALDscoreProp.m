function [varargout] = a_Get_BALDscoreProp(D,DIM,pv,pw,yrM,Gmap,z,mode,PbiasM,biasw)
%% still messy: different modes use different motor noise and lapse rate
% mode 1: likelihood and gradients (gradient not coded yet) (variable needed up to mode)
% mode 2: record BAS stat of D (variable needed up to mode)
%   Note! Moments.lpx = lpx, used for pattern likelihood in a_Model
% mode 2.5: record BAS maps of participants AL (yet to convolved with Motor noise)
% mode 2.6: record Max Ent maps of participants AL
% mode 3: simulate BAS with DIM, optimal (variable needed up to Gmap)
% mode 4: simulate BAS with Motor noise and bias (-Gmap)(use PbiasM=[]; baisw=zeros(25,1) for without bias)
% mode 5: simulate Max Ent with DIM, optimal (variable needed up to mode -Gmap)
% mode 6: heuristic: BAS+ Motor noise + bias + softmaxed log weights lw (see mode 4 for specifications)
% mode 7: posterior-dep order-indep fixation (variable needed up to Gmap)
% mode 7.5: posterior-dep order-dep fixation
% mode 8: same as mode 7, but use prob matching
% mode 9: type-sensitive heursitics with participants' saccade maps (like mode 7)
% mode 10: same as mode 9, but use prob matching
% mode 11: heursitic: gungho
% mode 12: same as mode 11, but use prob matching

%% set up parameters
npix=770;  ImageSize=20;  n=sqrt(size(z,1));  ss=20;  % screen size = [-10,10] cm
w=pv(9);
lsigo=pv(1);  sigo=exp(lsigo);
phit=pv(2);
lphis=pv(3);  phis=exp(lphis);
phis0=pv(4);
lbx=pv(5);  bx=exp(lbx);
shift = pv(12);
pvp=[w,sigo,phit,phis,phis0,bx,shift];  % parameters vector for perception
lbd=pv(6);  bd=exp(lbd);
lapse=pv(7);
prpa=pv(8);
prxpa = pv(10);
prCommon = pv(11);
pvd=[bd,lapse,prpa, prxpa, prCommon];  % parameters vector for decision

% initialize quantities
lbm=pw(1);  bm=exp(lbm); % softmax beta for motor choice
km=pw(2);  % lapse rate for motor choice
nm=pw(3);  % SD of motor noise in pixels
bias=pw(4);  % proportion of Pbias used
ntrial = D.Trials;
nrmax  = max(D.MaxRevealingTrial);

nz=size(z,1);

if(mode==1)
    E.Pxs = NaN(ntrial,nrmax);    
end

if(mode==2)
    E.RevealPosX = zeros(ntrial,nrmax);
    E.RevealPosY = zeros(ntrial,nrmax);
    E.mla = NaN(ntrial,nrmax);
    E.ml = NaN(ntrial,nrmax,3); %mlPA, mlSH, mlSV
    E.softmla = NaN(ntrial,nrmax);
%     E.BALDmap = NaN(ntrial,n*n,nrmax);  %this takes up a lot of space
%     E.BALDxs = NaN(ntrial,nrmax);
%     E.PTILExs = NaN(ntrial,nrmax);
end

if(mode==2.5 || mode==2.6)
    E.RevealPosX = zeros(ntrial,nrmax);
    E.RevealPosY = zeros(ntrial,nrmax);
    E.ml = NaN(ntrial,nrmax,3); %mlPA, mlSH, mlSV
    E.BALDxs = NaN(ntrial,nrmax);
    E.PTILExs = NaN(ntrial,nrmax);
    E.BALDxs = NaN(ntrial,nrmax);
    E.opt_xs = NaN(ntrial,nrmax); %euclidean dist from optimal point to xs
    %these take lots of memory
    %E.BFmap = NaN(ntrial,nz,nrmax);  %brute force
    E.BALDmap = NaN(ntrial,nz,nrmax);
    %need these two for BASmap-2015-09-09 or 14
    E.pred_mean = NaN(nz,3,ntrial, nrmax);
    E.pred_var = NaN(nz,3,ntrial, nrmax);
end

if(mode==5)
    ntrial = D.Trials;
    nrmax  = max(D.MaxRevealingTrial);    
    Dout=D;
    RevPosX=NaN(nrmax,ntrial);
    RevPosY=NaN(nrmax,ntrial);
    Dout.MaxRevealingTrial=zeros(ntrial,1);
    %Dout.map = NaN(ntrial,nz,nrmax);
    Dout.ml = NaN(ntrial,nrmax,3); %mlPA, mlSH, mlSV
end

if(mode==3 || mode==4 || mode==6)
    ntrial = D.Trials;
    nrmax  = max(D.MaxRevealingTrial);    
    Dout=D;
    RevPosX=NaN(nrmax,ntrial);
    RevPosY=NaN(nrmax,ntrial);
    Dout.MaxRevealingTrial=zeros(ntrial,1);
    % Dout.mla = NaN(ntrial,nrmax); % delete this?
    Dout.ml = NaN(ntrial,nrmax,3); %mlPA, mlSH, mlSV
    Dout.winner = NaN(ntrial,nrmax); %for mode==7
%     Dout.BALDxs = NaN(ntrial,nrmax);
%     Dout.PTILExs = NaN(ntrial,nrmax);
end

if(mode>=7)
    ntrial = D.Trials;
    nrmax  = max(D.MaxRevealingTrial);    
    Dout = D;
    RevPosX = NaN(nrmax,ntrial);
    RevPosY = NaN(nrmax,ntrial);
    Dout.MaxRevealingTrial = zeros(ntrial,1);
    Dout.ml = NaN(ntrial,nrmax,3); %mlPA, mlSH, mlSV
    Dout.winner = NaN(ntrial,nrmax);
end

if(mode==9 || mode==10)
   D = appendSaccadeMap(D);
end

% % setup plot
% scrsz = get(0,'ScreenSize');
% hf=figure('Name', 'Simulation Plot Window', 'position',[scrsz(3)/10 scrsz(4)/4 scrsz(3)/1.2 scrsz(4)/2],'menubar','none');

%% loop over trials
for trial=1:D.Trials

nr = D.MaxRevealingTrial(trial);
if (mode == 2 || mode == 2.5 || mode == 2.6)
    yD = D.RevealZ(trial,:);
else
    yimage = DIM(trial).yimage;
end

if(mode<3)        
xDx = 1+(npix-1)*(D.RevealPosX(trial,1:nr)/ImageSize+0.5);  % cm to pixel: from inverting last two lines of a_GET_BALD
xDy = 1-(npix-1)*(D.RevealPosY(trial,1:nr)/ImageSize-0.5);
xD = [xDx',xDy'];  
mla=NaN(nr,1);
softmla = NaN(nr,1);
matH_BALD=zeros(nz,nr);  xsH_BALD=zeros(1,nr);  matH_BF=zeros(nz,nr);
mlPA=NaN(nr,1);  mlSH=NaN(nr,1);  mlSV=NaN(nr,1);
for i=1:nr
    subxD = xD(1:i,:);
    subyD = yD(1:i)';
    yrin = yrM(trial,1:i)';  %simulating perceived noise
    if (mode==2)
        Moments.lpx = Gmap(trial,i,:);
        [softmla_t, ~, lw] = a_Model(subxD, subyD, pvp, Moments, pvd, yrin);
        softmla(i) = softmla_t;
    else    
        if(i==nr)
            zxs=[z;zeros(1,2)];  %add to z [0,0] as dummy    
        else
            zxs=[z;xD(i+1,:)];  %add to z the next revealing chosen
        end
        %[lw,matmu,matvar] = a_Model(subxD, yimage, pvp, [], pvd, yrin, zxs); % old version                
        [lw,matmu,matvar] = a_Model(subxD, subyD, pvp, [], pvd, yrin, zxs);
        %record to understand why MaxEnt isn't highest at image corners
        %2015-09-09
        E.pred_mean(:,:,trial,i) = matmu(1:nz,:);
        E.pred_var(:,:,trial,i) = matvar(1:nz,:);
        %BALD score calculation
        lwb=lwbound(lw); %to bound lw
        mu=reshape(matmu,numel(matmu),1);
        var=reshape(matvar,numel(matvar),1);
        if (mode==2.5)
            H_BALD = a_GetBALDGM_mex(lwb, mu, var);
            %bruteforce BAS score as a check
%             BF=zeros(nz,1);
%             w=exp(lwb)/sum(exp(lwb));  w2=w(2:3)/sum(w(2:3));
%             for indbf=1:nz+1
%                 mx=matmu(indbf,:);
%                 vx=matvar(indbf,:);
%                 BF(indbf) = GetEntropyGM_bruteforce(w,mx,vx)...
%                     -w(1)*GetEntropyGM_bruteforce(1,mx(1),vx(1))...
%                     -sum(w(2:3))*GetEntropyGM_bruteforce(w2,mx(2:3),vx(2:3));
%             end
        elseif (mode==2.6)
            H_BALD = a_GetMaxEnt_mex(lwb, mu, var);
        end
        H_BALD = max(0, H_BALD);
        matH_BALD(:,i)=H_BALD(1:end-1);
        %matH_BF(:,i)=BF(1:end-1);
        xsH_BALD(1,i)=H_BALD(end);
    end
    mla(i)=exp(lw(1))/sum(exp(lw));
    mlPA(i)=exp(lw(1))/sum(exp(lw));
    mlSH(i)=exp(lw(2))/sum(exp(lw));
    mlSV(i)=exp(lw(3))/sum(exp(lw));
end
if (mode==1)
    matP = matH_BALD.^(bm);  % softmax
    for i=1:nr-1    
        Pmap = matP(:,i)/sum(matP(:,i));
        Pmap = reshape(conv2(reshape(Pmap(:,i),n,n),Gmap,'same'),n*n,1); % convolve with motor noise
        Pmap = (1-km)*Pmap + km/n/n;  %lapse rate
        indx = find(z(:,1)>xD(i+1,1),1,'first');
        indx = min(indx,n*n-n);
        indy = find(z(indx:indx+n,2)>xD(i+1,2),1,'first');
        Pxs = Pmap(indx+indy)/sum(Pmap);       
        E.Pxs(trial,i+1)=Pxs;
    end
end        
if(mode==2)
    E.RevealPosX(trial,1:nr) = xD(:,1)';
    E.RevealPosY(trial,1:nr) = xD(:,2)';
    E.mla(trial,1:nr) = mla;
    E.ml(trial,1:nr,1) = mlPA;
    E.ml(trial,1:nr,2) = mlSH;
    E.ml(trial,1:nr,3) = mlSV;
    E.softmla(trial,1:nr) = softmla;
%     E.BALDmap(trial,:,1:nr)=matH_BALD;
%     E.BALDxs(trial,1:nr)=xsH_BALD;
%     Not needed for DataPrep
%     for i=1:nr-1
%         Bmap = matH_BALD(:,i); %raw BALD score
%         Bxs = xsH_BALD(i);
%         sBmap = sort(Bmap);
%         PTxs = find(sBmap<=Bxs,1,'last')/n/n*100;  %caused a problem for sig0=0.1, prob because nan sort of stuff
%         E.PTILExs(trial,i+1) = PTxs;            
%     end
end
if(mode==2.5 || mode==2.6)
    for i=1:nr-1
        %Convolve with motor noise--MC approx
        Bmap = matH_BALD(:,i); %raw BALD score
        Cmap = Bmap; %MotorNoiseConv(Bmap,xDp,z,npix,ss,nz);
        %optimal spot
        [~,ind]=max(Cmap);
        dvec = z(ind,:)-xD(i+1,:);
        E.opt_xs(trial,i+1) = norm(dvec); %is it +1?
        %PTILExs: may be slow because of the sort
        % indxs %find closet ind of xs in z now that Cmap is convolved
        Cxs = xsH_BALD(i); %Cmap(ind);
        sCmap = sort(Cmap);
        PTxs = find(sCmap<=Cxs,1,'last')/n/n*100;  %caused a problem for sig0=0.1, prob because nan sort of stuff
        if (isempty(PTxs))
            E.PTILExs(trial,i+1) = 100;
        else        
            E.PTILExs(trial,i+1) = PTxs;
        end
        %record
        E.RevealPosX(trial,1:nr) = xD(:,1)';
        E.RevealPosY(trial,1:nr) = xD(:,2)';
        E.ml(trial,1:nr,1) = mlPA;
        E.ml(trial,1:nr,2) = mlSH;
        E.ml(trial,1:nr,3) = mlSV;
        E.BALDmap(trial,:,i+1)=Cmap;
        E.BALDxs(trial,1:i+1)=Cxs;
    end
    %E.BFmap(trial,:,:)=matH_BF;
end
end
    
% simulation mode starts
if(mode==3 || mode==4 || mode==5 || mode==6)
BALDxD = [0,0] + npix/2 + 1 + randn(1,2)*1e-2;  % (despite) in pixel, 1st column still x, 2ns still y
count=1;        
    while ( count<=nr )    
        count=count+1;
        yrin = yrM(trial, 1:count-1)';  %simulating perceived noise
        [lw,matmu,matvar] = a_Model(BALDxD, yimage, pvp, [], pvd, yrin, z);
        mu=reshape(matmu,numel(matmu),1);
        var=reshape(matvar,numel(matvar),1);
        if(mode==3 || mode==4)
            % Huber approximated BAS score; still can return to visited place;
            % likely related to perception noise calculation: see a_Model
            lwb = lwbound(lw); %to bound lw
            H_BALD = a_GetBALDGM_mex(lwb, mu, var);
            H_BALD = max(0, H_BALD);
        elseif(mode==5)
            %MaxEnt calculation
            lwb = lwbound(lw); %to bound lw
            H_BALD = a_GetMaxEnt_mex(lwb, mu, var);
            H_BALD = max(0, H_BALD);
        elseif(mode==6)
            exponent = 5; % just a random choice
            slw = sofmaxLogWeights(lw, exponent);
            slwb = lwbound(slw); %to bound lw
            H_BALD = a_GetBALDGM_mex(slwb, mu, var);
            H_BALD = max(0, H_BALD);            
        end        
        if(isinf(bm))  %hardmax; changed but not tested!
            %[~,ind]=max(H_BALD);  %the original
            %matP=zeros(n); matP(ind)=1;  %the original: is it fluke that this worked?
            maxH = max(H_BALD);
            matP = zeros(nz,1);
            matP(H_BALD == maxH) = 1;
            matP = matP/sum(matP);
        else  %softmax
            matP = H_BALD.^(bm);
            matP = matP/sum(matP);
        end
        intBmap = cumsum(matP);
        rseed=rand()*intBmap(end);
        ind = find(intBmap>rseed,1,'first');
        BALDxnew = z(ind,:) + randn(1,2)*1e-2; %z is in (x,y)
        if (mode==4 || mode==5 || mode==6)
            if(rand()>biasw(count)) %sample from BAS with motor noise
                BALDxn = BALDxnew;
                BALDxnew = MotorNoise(BALDxD(count-1,:),BALDxn,npix,ss); %add realistic motor noise
            else %sample from bias map
                Pbias=PbiasM(:,count);
                intPbias = cumsum(Pbias);
                rseed=rand()*intPbias(end);
                ind = find(intPbias>rseed,1,'first');
                BALDxnew = z(ind,:)+ randn(1,2)*1e-2;                
            end
        end
        BALDxD = [BALDxD; BALDxnew];

%         %bruteforce BAS score as a check
%         vx=matvar(ind,:); mx=matmu(ind,:); w=exp(lw)/sum(exp(lw));  w2=w(2:3)/sum(w(2:3));
%         B=GetEntropyGM_bruteforce(w,mx,vx)-w(1)*GetEntropyGM_bruteforce(1,mx(1),vx(1))-sum(w(2:3))*GetEntropyGM_bruteforce(w2,mx(2:3),vx(2:3)); 
%         B=max(0,B);
        
        %Dout.mla(trial,count-1)=exp(lw(1))/sum(exp(lw)); %delete this?
        Dout.ml(trial,count-1,1)=exp(lw(1))/sum(exp(lw));
        Dout.ml(trial,count-1,2)=exp(lw(2))/sum(exp(lw));
        Dout.ml(trial,count-1,3)=exp(lw(3))/sum(exp(lw));

%         Not needed for DataPrep
%         if (mode==3) %need to update ind after MotorNoise to be used in mode==4
%             Bxs=H_BALD(ind);
%             Dout.BALDxs(trial,count) = Bxs;  
%             sBmap = sort(H_BALD);
%             PTxs = find(sBmap<=Bxs,1,'last')/n/n*100;
%             Dout.PTILExs(trial,count) = PTxs;
%         end

%         if (mode==5)
%             Dout.map(trial,:,count-1) = H_BALD;
%         end
        
%         % plots
%         colormap(gray);
%          px=linspace(1,npix,n);
%          subplot(1,3,1); imagesc(px, px, reshape(H_BALD,n,n)); hold on;         
%             plot( BALDxD(:,1), BALDxD(:,2), 'r.');
%             plot( BALDxnew(1), BALDxnew(2), 'ro','MarkerFaceColor','g', 'MarkerSize',5); hold off;
%             axis square;
%          
%          subplot(1,3,2); imagesc(px, px, reshape(matP,n,n)); hold on;         
%             plot( BALDxD(:,1), BALDxD(:,2), 'r.'); 
%             plot( [BALDxn(1),BALDxnew(1)], [BALDxn(2),BALDxnew(2)], 'y.-');
%             %plot( BALDxnew(1), BALDxnew(2), 'ro','MarkerFaceColor','g', 'MarkerSize',5); 
%             hold off;
%             axis square;
% 
%          subplot(1,3,3); imagesc(yimage); hold on;         
%             plot( BALDxD(:,1), BALDxD(:,2), 'r.'); 
%             plot( BALDxnew(1), BALDxnew(2), 'ro','MarkerFaceColor','g', 'MarkerSize',5);  hold off;
%             axis square;
%             
%         % allows pressing q to exit the loop
%         if strcmp(get(hf,'currentcharacter'),'q')
%             break
%         end
% 
%         drawnow
%         pause
    end
    Dout.MaxRevealingTrial(trial) = nr;
    RevPosX(1:nr,trial) = BALDxD(1:nr,1);
    RevPosY(1:nr,trial) = BALDxD(1:nr,2);
end
% simulation mode ends

% active sensing heuristics begins
if (mode>=7)
    moveType(:,1) = D.ImageID(:,2)==20; %pa
    moveType(:,2) = D.ImageID(:,2)==6;  %sh
    moveType(:,3) = D.ImageID(:,2)==30; %sv
    all = ones(600, 1);
    if(mode==7 || mode==7.5 || mode==9 || mode==11)
        softmaxRate = inf; % inf: hardmax, 1:prob matching
    elseif(mode==8 || mode==10 || mode==12)
        softmaxRate = 1; % inf: hardmax, 1:prob matching
    end
    count = 1;
    if (mode==7 || mode==7.5 || mode==8)
        BALDxD = sampleRevPos(D, all);
    else
        BALDxD = [0,0] + npix/2 + 1 + randn(1,2)*1e-2;
    end
    while ( count<nr )
        count=count+1;
        yrin = yrM(trial, 1:count-1)';  %simulating perceived noise
        [~, ~, lw] = a_Model(BALDxD, yimage, pvp, [], pvd, yrin);
        w = exp(lw);
        winner = winnerIs(w, softmaxRate);
        if(mode==7 || mode==8)
            BALDxnew = sampleRevPos(D, moveType(:,winner));
        elseif(mode==7.5)
            BALDxnew = sampleOrderedRevPos(D, count, moveType(:,winner));
        elseif(mode==9 || mode==10)
            posCur = BALDxD(count-1, :);
            BALDxnew = sampleSaccade(posCur, D, moveType(:,winner));
        elseif(mode==11 || mode==12)
            posCur = BALDxD(count-1, :);                
            BALDxnew = gungHo(winner, posCur);
            BALDxn = BALDxnew;
            npix = 770;
            ss = 20;
            BALDxnew = MotorNoise(posCur, BALDxn, npix, ss);
        end
        BALDxD = [BALDxD; BALDxnew];
        Dout.winner(trial,count-1) = winner;
    end
    Dout.MaxRevealingTrial(trial) = nr;
    RevPosX(1:nr,trial) = BALDxD(1:nr,1);
    RevPosY(1:nr,trial) = BALDxD(1:nr,2);    
end
% active sensing heuristics ends

% disp(sprintf('trial=%d;',trial));
end

% outputs
if(mode==1)
    varargout={-nansum(nansum(nanmax(-1000,log(E.Pxs))))};  % neg log liklihood
elseif(mode==2 || mode==2.5 || mode==2.6)
    varargout={E};
elseif(mode==3 || mode==4 || mode==5 || mode==6 || mode>=7)
    Dout.RevealPosX = (RevPosX'-1)/(npix-1)*ss-ss/2;  %Dout.RevealPosX(Dout.RevealPosX>9.9)=9.9;  Dout.RevealPosX(Dout.RevealPosX<-9.9)=-9.9;
    Dout.RevealPosY = (1-RevPosY')/(npix-1)*ss+ss/2;  %Dout.RevealPosY(Dout.RevealPosY>9.9)=9.9;  Dout.RevealPosY(Dout.RevealPosY<-9.9)=-9.9;
    varargout={Dout};
elseif(mode==5)
    Dout.RevealPosX = RevPosX';
    Dout.RevealPosY = RevPosY';
    varargout={Dout};
end
    
end

%% ########################################################################

%% for mode==8 setup 
function Dout = appendSaccadeMap(Din)
    x = Din.RevealPosX;
    y = Din.RevealPosY;
    x(x==0) = nan;
    y(y==0) = nan;
    dx = diff(x,1,2);
    dy = diff(y,1,2);
    Dout = Din;
    Dout.SaccadeMapX = dx;
    Dout.SaccadeMapY = dy;
end

%% bound lw to avoid numerical error in a_GetBALDGM_mex(lwb, mu, var)
function lwn=lwbound(lw)
    c=min(lw);
    if(c<-30)
        lwn=max(lw,-30);
        lnor=log(sum(exp(lwn)));
        lwn=lwn-lnor;
    else
        lwn=lw;
    end
end

%% signal-dependent endpoint variability
% didn't include saccade angle dependence
function [sig_R, sig_a]=sig_endpoint(R)
mode=2;  
if (mode==1) % vanBeer07
    m_R=0.025;  a_R=0.15;   b_R=0.11; %endpoint radius parameters    
    m_a=0.015;  a_a=0.044;  b_a=0.12; %endpoint angle parameters
    sig_R = R*(m_R+a_R*exp(-b_R*R));
    sig_a = R*(m_a+a_a*exp(-b_a*R));    
elseif (mode==2)
%     % Krappmann98 + vanBeers07 
%     % for amp: 0.03*R + 0.23 + 1.51
%     sig_R = (20/26.8)*(0.03*(R*26.8/20) + 0.23 + 1.51);  
%     % for dir: 0.014*R + 0.085 + 0.626
%     sig_a = (20/26.8)*(0.014*(R*26.8/20) + 0.085 + 0.626);
    % from own experimental measurement (both input and output is in cm)
    % old sdTan: intercept = 0.3345, slope = 0.1011 (all in cm)
    % old sdNor: intercept = 0.2664, slope = 0.004257 (all in cm)
    % sdTan: intercept = 0.2989, slope = 0.1407 (all in cm)
    % sdNor: intercept = 0.2823, slope = 0.01341 (all in cm)
    % sig_R = 0.14*R + 0.30;  
    % sig_a = 0.013*R + 0.28;
    % madsdTan: intercept = 0.294, slope = 0.1301 (all in cm) 2015-05-29
    % madsdNor: intercept = 0.2724, slope = 0.01098 (all in cm)
    sig_R = 0.13*R + 0.29;  
    sig_a = 0.011*R + 0.27;
end
end

%% caculate motor noise: 
% an approximation that uses an rotated oval covariance 
% to describe the projected mag and direction variance
function xDnew=MotorNoise(xDp,xDn,npix,ss)
    %step1: from pixel to cm
    d=42; %distance from eye to fixation cross in cm
    xDp(1) = (xDp(1)-1)/(npix-1)*ss-ss/2;  xDp(2) = (1-xDp(2))/(npix-1)*ss+ss/2;
    xDn(1) = (xDn(1)-1)/(npix-1)*ss-ss/2;  xDn(2) = (1-xDn(2))/(npix-1)*ss+ss/2;
    %step2: define fixation vectors
    fp=[xDp(1),xDp(2),d]; %fixation previous in x,y,z, cm
    fn=[xDn(1),xDn(2),d]; %fixation now in x,y,x, cm
    %step3: find angle between fixation vectors and motor noise
    % R=acos( (fp*fn')/sqrt(fp*fp')/sqrt(fn*fn') )*180/pi;  %R in degree
    Rcm = norm(xDn-xDp);  % R in cm
    [sig_R, sig_a]=sig_endpoint(Rcm);
    %step4: find difference vector and its perpendicular vector
    dv=fn-fp;  dv(3)=[];  dv=dv/sqrt(dv*dv'); %unit difference vector
    du(1)=-dv(2);  du(2)=dv(1); %unit perpendicular vector
    %step5: compute covariance matrix in x,y
    U=[dv',du']; L=[sig_R,0;0,sig_a].^2;
    CM=U*L*U';
    %step6: sample noise from covariance
    mn = mvnrnd([0,0],CM);
    %step7: add bias
    % errTan: intercept = 0.2643, slope = -0.2284 (all in cm)
    bias = -0.23*Rcm + 0.26;
    bn = bias*dv;
    %step8: from cm to pixel
    xDnew=xDn+mn+bn;
    xDnew=max(-10,xDnew);  xDnew=min(10,xDnew); %bounds
    xDnew(1) = 1+(npix-1)*(xDnew(1)/ss+0.5);
    xDnew(2) = 1-(npix-1)*(xDnew(2)/ss-0.5);
    %mn = mn*npix/ss;
end

%% to convolve a map with the saccadic noise:
% function Cmap=MotorNoiseConv(Bmap,xDp,z,npix,ss,nz)
% 
% nsamp=100;
% Cmap = zeros(nz,1); %size=Bmap
% 
% for ii=1:nz
%     xsamp=MotorNoise(xDp,z(ii,:),npix,ss,nsamp); %add realistic motor noise
%     ind=;%get xsamp ind to get their Bmap values
%     %what happens when sample out of bound?
%     Bsamp = Bmap(ind,:);
%     Cmap(ii) = mean(Bsamp);
% end
% 
% end

%% for mode==6 (not used anymore)
function slw = sofmaxLogWeights(lw, exponent)
    slw = exponent*lw;
    lognor = log(sum(exp(slw)));
    slw = slw - lognor;
end

%% choose winner
function winner = winnerIs(w, softmaxRate)
    if(isinf(softmaxRate))
        [~, winner] = max(w);
    elseif(softmaxRate == 1)
        wCum = cumsum(w);
        rseed = rand();
        ind = find(wCum > rseed, 1, 'first');
        winner = ind;
    end
end

%% order-independent sample from revMap
function revPos = sampleRevPos(D, moveType)
    % remove irrelevant move types and nan elements in x
    posX = D.RevealPosX;
    posX(moveType,:) = [];
    posX = posX(:);
    posX(posX==0) = [];
    % remove irrelevant move types and nan elements in y
    posY = D.RevealPosY;
    posY(moveType,:) = [];
    posY = posY(:);
    posY(posY==0) = [];
    % sample from valid moves
    ind = ceil(rand()*numel(posX));    
    [posXnew, posYnew] = cm2Pixel(posX(ind), posY(ind));
    revPos = [posXnew, posYnew];    
end

%% order-dependent sample revealings from participant's data (variants of mode==7)
function revPos = sampleOrderedRevPos(D, irev, moveType)
    indPool = find(D.RevealPosX(:,irev)~=0 & moveType); %extract relevant revealed ind
    numInPool = numel(indPool);
    randInd = ceil(rand()*numInPool);
    ind = indPool(randInd); %use only revealed ind
    posX = D.RevealPosX(ind,irev);
    posY = D.RevealPosY(ind,irev);
    [posXnew, posYnew] = cm2Pixel(posX, posY);
    revPos = [posXnew, posYnew];    
end

%% order-independent sample from saccadeMap
function posNew = sampleSaccade(posCur, D, moveType) %xCur in pixel
    % remove irrelevant move types and nan elements in x
    sacX = D.SaccadeMapX;
    sacX(moveType,:) = [];
    sacX = sacX(:);
    sacX(isnan(sacX)) = [];
    % remove irrelevant move types and nan elements in y
    sacY = D.SaccadeMapY;
    sacY(moveType,:) = [];
    sacY = sacY(:);
    sacY(isnan(sacY)) = [];
    % remove saccades that go outside the boundaries
    [xcm, ycm] =  pixel2Cm(posCur(1), posCur(2)); 
    xN = xcm + sacX;
    yN = ycm + sacY;
    outside = (xN<-10 | xN>10 | yN<-10 | yN>10);
    xN(outside) = [];
    yN(outside) = [];
    % sample from valid moves
    ind = ceil(rand()*numel(xN));
    [xNew, yNew] = cm2Pixel(xN(ind), yN(ind));
    posNew = [xNew, yNew];    
end
        
%% heursitic gungho
function revPos = gungHo(winner, curPos) % curPos in pixel
    stepCm = 2;
    npix = 770;
    ImageSize = 20;
    step = round((npix-1)*(stepCm/ImageSize));
    margin = 0; % non-zero margin does not work with motor noise

    if (winner == 1)  %PA
        moves(1,:) = [step, 0];     %go right
        moves(2,:) = [-step, 0];    %go left
        moves(3,:) = [0, step];     %go up
        moves(4,:) = [0, -step];    %go down
    elseif (winner == 2)  %SH
        moves(1,:) = [step, 0];
        moves(2,:) = [-step, 0];
    elseif (winner == 3)  %SV
        moves(1,:) = [0, step];
        moves(2,:) = [0, -step];
    end
    tryPos = repmat(curPos, size(moves,1), 1) + moves;
    tryPos(tryPos(:,1) < 1+margin, :) = [];
    tryPos(tryPos(:,2) < 1+margin, :) = [];
    tryPos(tryPos(:,1) > npix-margin, :) = [];
    tryPos(tryPos(:,2) > npix-margin, :) = [];
    ind = ceil(rand()*size(tryPos,1));
    revPos = tryPos(ind,:);
end

%% 
function [xcm, ycm] = pixel2Cm(xp, yp)
    npix = 770;
    ss = 20;
    xcm = (xp-1)/(npix-1)*ss - ss/2;
    ycm = (1-yp)/(npix-1)*ss + ss/2;
end

%% 
function [xp, yp] = cm2Pixel(xcm, ycm)
    npix = 770;
    ss = 20;
    xp = 1 + (npix-1)*(xcm/ss + 0.5);
    yp = 1 - (npix-1)*(ycm/ss - 0.5);
end


