%4 - image type patch, sh sv, average
%participant 1 2 3 av
%condition 1 data 2 bas
clear all
clf
BAS_FLAG=0;
load REV
nrep=10; % repetitions of splitting i.e. bootstrap samples
nrev=25; % maximum number of revealings
minrev=5; % we can go from here as is limited by number of 25th revealings
nimage=3;
sd=1.5*20/27.8;

profile off
%%
for subj=1:4
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
        for k=1:data.Trials %trials
            tcount=tcount+1;
            trial.image(tcount)=j;
            trial.num(tcount)=tcount;
            tria.maxrev(tcount)=data.MaxRevealingTrial(k);
            trial.weight(tcount)=ceil(data.MaxRevealingTrial(k)/5); % weighting for sampling
            
            for rev=1:data.MaxRevealingTrial(k) % loop over revealings for the trial
                count=count+1;
                
                fix.x(count)=data.RevealPosX(k,rev);
                fix.y(count)=data.RevealPosY(k,rev);
                fix.maxrev(count)=data.MaxRevealingTrial(k); % revealings in each trial
                fix.im_type(count)=j;
                fix.rev(count)=rev; % revealing number
                fix.trial(count)=tcount; % trial number
                fix.weight(count)=ceil(rev/5); % weighting for sampling
            end
        end
    end
    
    
    
    %%
    S=bindata([trial.image ;tria.maxrev]',ones(cols(trial.image),1));
    trialmin=floor(min(S.n(:))/2);
    nselect=trialmin*25;
    
    ntotal=nimage*nselect;
    
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
        end
    end
    
    %%
    reshuffle=0
    for rep=1:nrep  %splittings
        rep
        % split into two groups of trials randomly
        for j=1:3
            % extract randomly selected  trials but equal with regard
            % to number of revealings in each half
            if BAS_FLAG
                ns=length(trial.image==j);
                n1=trialmin;
                
                n{j,1}=1:n1;  %trials in first split
                n{j,2}=n1+1:2*trialmin; %trials in 2nd split
            else
                n{j,1}=[];
                n{j,2}=[];
                for k=[5:5:25]  %make this even for each
                    s=find(trial.image==j & tria.maxrev==k);
                    nperm= datasample(s,2*trialmin,'Replace',false);
                    n{j,1}=[n{j,1}; nperm(1:trialmin)'];  %trials in first split
                    n{j,2}=[n{j,2}; nperm(trialmin+1:end)']; %trials in second split
                end
            end
            
            for k=1:2
                s = ismember(fix.trial,n{j,k}) ; %find the selected trials for this half
                p=fix.weight .* s; %choose all the fixations from these trials with weights
                
                for rev=minrev:nrev
                    if reshuffle  %this is for reshuffling
                        ps=p; %extract simply by weights with no regard to order
                        e=datasample(1:length(p),nselect,'Replace',false,'Weights',ps+eps);
                    else
                        e=[];
                        for i=1:rev
                            ps=p.* (fix.rev==i); %choose all the fixations from these trials with weights
                            e=[e datasample(1:length(p),ceil(nselect/rev),'Replace',false,'Weights',ps+eps)];
                        end
                        e=datasample(e,nselect,'Replace',false);
                    end
                    
                    mn{rep,j,k,rev}=e;
                end
            end
        end
    end
    
    
    
    
    %%
    for rep=1:nrep  %splittings
        for rev=minrev:nrev
            s1=[mn{rep,1,1,rev} ; mn{rep,2,1,rev}; mn{rep,3,1,rev}]';
            s2=[mn{rep,1,2,rev} ; mn{rep,2,2,rev}; mn{rep,3,2,rev}]';
            c(:,:,rev,rep)= mapcorr(s1,s2,dist,maskprod);
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
    
    % plot
    mysubplot(4,2,subj,reshuffle+1)
    hold on
    
    [Hp,Hs]=shadeplot(minrev:25,mean(within(minrev:end,:)'),std(within(minrev:end,:)'),'g','g')
    set(Hs,'FaceAlpha',0.2)
    
    [Hp,Hs]=shadeplot(minrev:25,mean(between(minrev:end,:)'),std(between(minrev:end,:)'),'r','r')
    set(Hs,'FaceAlpha',0.2)
    pause(0.1)
end



%%
for k=1:8
    subplot(4,2,k)
    axis([minrev 25 -0.4 0.4])
end
profile report

%export_fig shuffle.pdf

