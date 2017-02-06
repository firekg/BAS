%4 - image type patch, sh sv, average
%participant 1 2 3 av
%condition 1 data 2 bas
% ssh -Y -4 -o 'GSSAPIKeyExchange no' hinton
% /misc/apps/matlab/matlabR2015a/bin/matlab -nodesktop -nodisplay -r "subj=1;" < analysis_hinton.m  &> hinton1.out &


% addpath(genpath('~/Dropbox/matlab/'))
% addpath(genpath('/users/lightspeed'));
% addpath(genpath('/users/fastfit'));
% startup

for subj = 1:3

BAS_FLAG = 0;
load REV_pr
% load REV
nrep=1000; % repetitions of splitting i.e. bootstrap samples
nrev=25; % maximum number of revealings
%minrev=5; % we can go from here as is limited by number of 25th revealings
minrev=1;
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

%fprintf('Got to data and fix.\n');

%%
S=bindata([trial.image_type ;trial.maxrev]',ones(cols(trial.image_type),1));
trialmin=floor(min(S.n(:))/2); % number of trials we will sample into each half
%nselect=trialmin*21; %number of fixations to sample in each half for correlations
nselect=trialmin*5;
ntotal=nimage*nselect;

%fprintf('Got to set up parameters.\n');

%% presave stuff

[xa,xb]=meshgrid(fix.x,fix.x);
[ya,yb]=meshgrid(fix.y,fix.y);
dist=normpdf(xa-xb,0,sqrt(2)*sd).*normpdf(ya-yb,0,sqrt(2)*sd);

mask=repmat(-1/ntotal,[ntotal nimage]);
for i=1:nimage
    mask((i-1)*nselect+1:i*nselect,i)=((nimage-1)/ntotal);
end
maskprod=zeros(ntotal,ntotal,nimage,nimage);
for i=1:nimage
    for j=1:nimage
        maskprod(:,:,i,j)=mask(:,i)*mask(:,j)';
        w=maskprod(:,:,i,j);
        newmask{i,j}=w(:)';
    end
end


%fprintf('Finish presave stuff.\n');

%%
% profile clear
% profile on
reshuffle=0
for rep=1:nrep  %splittings
    rep
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
            %s = ismember(fix.trial,n{j,k}) ; %Scott trying
            
            for rev=minrev:nrev
                if reshuffle  %this is for reshuffling
                    p=fix.weight .* s; %choose all the fixations from these trials with weights
                    
                    ps=p; %extract simply by weights with no regard to order
                    e=datasample(1:length(p),nselect,'Replace',false,'Weights',ps+eps);
                else
                    e=[];
                    for i=1:rev
                        ps=s(fix.rev(s)==i); %choose all the fixations from these trials with weights
                        e=[e randsample(ps,ceil(nselect/rev))];
                    end
                    e=randsample(e,nselect);
                end
                mn{rep,j,k,rev}=e;
            end
        end
    end
end

%fprintf('Made mn.\n');


%%
for rep=1:nrep  %splittings
    rep
    for rev=minrev:nrev
        s1=[mn{rep,1,1,rev} ; mn{rep,2,1,rev}; mn{rep,3,1,rev}]';
        s2=[mn{rep,1,2,rev} ; mn{rep,2,2,rev}; mn{rep,3,2,rev}]';
        c(:,:,rev,rep)= new_mapcorr(s1,s2,dist,newmask);
    end
end

%%
%combine  within and between
within=zeros(nrev,nrep);
between=zeros(nrev,nrep);

for j=1:3
    for k=1:3
        for rev=minrev:nrev
            if j==k
                within(rev,:) =within(rev,:)+squeeze(c(j,k,rev,:))';
            else
                between(rev,:)=between(rev,:)+squeeze(c(j,k,rev,:))';
            end
        end
    end
end

%normalise
within=within/3;
between=between/6;

str=sprintf('save subj%iprBas1to25 minrev within between',subj);
eval(str);

clear;
end
