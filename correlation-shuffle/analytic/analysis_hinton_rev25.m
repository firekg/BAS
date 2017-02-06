%4 - image type patch, sh sv, average
%participant 1 2 3 av
%condition 1 data 2 bas
% ssh -Y -4 -o 'GSSAPIKeyExchange no' hinton
% /misc/apps/matlab/matlabR2015a/bin/matlab -nodesktop -nodisplay -r "subj=1;" < analysis_hinton.m  &> hinton1.out &

% uncomment for hinton
% addpath(genpath('~/Dropbox/matlab/'))
% addpath(genpath('/users/lightspeed'));
% addpath(genpath('/users/fastfit'));
% startup

for subj = 1:3

BAS_FLAG = 1;
load REV_pr_Maxent
nrep=1000; % repetitions of splitting i.e. bootstrap samples
nrev=25; % maximum number of revealings
nimage=3;
sd=1.5*20/27.8;
% sd=3*20/27.8;

%%

count=0;
tcount=0;

%sets up data an sizes for input
for j=1:3  %patchy  stripy-h stripy-v
    j
    if BAS_FLAG
        data1=REV(j,subj,1);  %pull out data  subjects
        data2=REV(j,subj,2);  %pull out data   BAS
        
        data.MaxRevealingTrial=[data1.MaxRevealingTrial; data2.MaxRevealingTrial];
        data.RevealPosX =      [data1.RevealPosX      ;  data2.RevealPosX];
        data.RevealPosY=       [data1.RevealPosY      ;  data2.RevealPosY];
    else
        data=REV(j,subj,1);  %pull out data
    end
        
    %save all fixations in one big array
    for k=1:length(data.MaxRevealingTrial) %trials
        tcount=tcount+1;
        trial.image_type(tcount)=j;  %trial data
        trial.maxrev(tcount)=data.MaxRevealingTrial(k);
        trial.weight(tcount)=ceil(data.MaxRevealingTrial(k)/5); % weighting for sampling

        for rev=1:data.MaxRevealingTrial(k) % loop over revealings for the trial
            count=count+1;
            fix.x(count)=data.RevealPosX(k,rev);
            fix.y(count)=data.RevealPosY(k,rev);
            fix.maxrev(count)=data.MaxRevealingTrial(k); % revealings in each trial
            fix.rev(count)=rev; % revealing number
            fix.trial(count)=tcount; % trial number
            fix.weight(count)=ceil(rev/5); % weighting for sampling
        end
    end
end


%%

S=bindata([trial.image_type ;trial.maxrev]',ones(cols(trial.image_type),1));
trialmin=floor(min(S.n(:))/2); % number of trials we will sample into each half

% hard-wired from incorrect
% trialmin = 15;

nselect = repmat(trialmin,1,5)*[5:5:25]'; %number of fixations to sample in each half for correlations
ntotal=nimage*nselect;


%%
reshuffle=0
for rep=1:nrep  %splittings
    fprintf('parti%d sampling rep %d;\n', subj, rep);
    % split into two groups of trials randomly
    for j=1:3
        % extract randomly selected  trials but equal with regard
        
        n{j,1}=[];
        n{j,2}=[];
        for k=[5:5:25]  %make this even for each
            s=find(trial.image_type==j & trial.maxrev==k);
            
            if BAS_FLAG
                nperm= randsample(s(1:length(s)/2),trialmin);
                n{j,1}=[n{j,1}; nperm'];  %trials in first split
                nperm= randsample(s(1+length(s)/2:end),trialmin);
                n{j,2}=[n{j,2}; nperm']; %trials in second split
            else
                nperm= randsample(s,2*trialmin);
                n{j,1}=[n{j,1}; nperm(1:trialmin)'];  %trials in first split
                n{j,2}=[n{j,2}; nperm(trialmin+1:end)']; %trials in second split
            end
        end
        
        for k=1:2
            s = find(ismember(fix.trial,n{j,k})) ; %find the selected trials for this half
            indp = randperm(length(s));
            indp(nselect+1:end) = [];
            indtr = s(indp);
            indrev = fix.rev(s(indp));
            [B,I] = sort(indrev);
            mnRev{rep,j,k} = B;
            mn{rep,j,k} = indtr(I);
        end
                
    end
end
       

%% presave stuff

fprintf('subj%d trailmin=%d\n', subj, trialmin); 

% [xa,xb]=meshgrid(fix.x,fix.x);
% d1 = normpdf(xa-xb,0,sqrt(2)*sd);
% clearvars xa xb;
% [ya,yb]=meshgrid(fix.y,fix.y);
% d2 = normpdf(ya-yb,0,sqrt(2)*sd);
% clearvars ya yb;
% dist = d1.*d2;
% clearvars d1 d2;

L = length(fix.x);
s = 1000;
nB = ceil(L/s)
dist = zeros(L);
size(dist)

for iB = 1:nB
    for jB = 1:nB
        if iB<nB
            blockindi = (iB-1)*s+1:iB*s;
        else
            blockindi = (iB-1)*s+1:L;
        end
        if jB<nB
            blockindj = (jB-1)*s+1:jB*s;
        else
            blockindj = (jB-1)*s+1:L;
        end
        [xa,xb] = meshgrid(fix.x(blockindj), fix.x(blockindi));
        [ya,yb] = meshgrid(fix.y(blockindj), fix.y(blockindi));
        dist(blockindi, blockindj) = normpdf(xa-xb,0,sqrt(2)*sd).*normpdf(ya-yb,0,sqrt(2)*sd);
    end
    iB
end
        
% mnRev is the same everytime, good
r = mnRev{1,1,1}';
weight = repmat( 1./(6 - ceil(r/5)), 3, 3);

mask=repmat(-1/ntotal,[ntotal nimage]);
for i=1:nimage
    mask((i-1)*nselect+1:i*nselect,i) = ((nimage-1)/ntotal);
end

mask = mask.*weight;

maskprod=zeros(ntotal,ntotal,nimage,nimage);
for i=1:nimage
    for j=1:nimage
        maskprod(:,:,i,j)=mask(:,i)*mask(:,j)';
        w=maskprod(:,:,i,j);
        newmask{i,j}=w(:)';
    end
end

%%
for rep=1:nrep  %splittings
    fprintf('parti%d correlating rep %d;\n', subj, rep);
    s1=[mn{rep,1,1} ; mn{rep,2,1}; mn{rep,3,1}]';
    s2=[mn{rep,1,2} ; mn{rep,2,2}; mn{rep,3,2}]';
    c(:,:,rep)= new_mapcorr(s1,s2,dist,newmask);
end

%%
%combine  within and between
within=zeros(nrep,1);
between=zeros(nrep,1);

for j=1:3
    for k=1:3
        if j==k
            within = within + squeeze(c(j,k,:));
        else
            between = between + squeeze(c(j,k,:));
        end
    end
end

%normalise
within=within/3;
between=between/6;

%uncomment/m for hinton
str = sprintf('save subj%iprMaxEntRev25 within between',subj);
eval(str);

clear;
end


