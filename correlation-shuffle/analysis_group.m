%4 - image type patch, sh sv, average
%partcipant 1 2 3 av
%condition 1 data 2 bas

clear all
clf
load REV
nrep=10; % repetitions of splitting i.e. bootstrap samples
nrev=25; % maximum number of revealings
ndim=77; %size of grid to bin saccades before smoothing
xx=linspace(-15,15,ndim);  %grid
 minrev=1; % we can go from here as is limited by number of 25th revealings
 nmin=150;
%minrev=5; % we can go from here as is limited by number of 25th revealings
%nmin=375;
sd=1;
%%
for subj=1
    clear sel fix ntrial data
    %sets up data an sizes for input
    for j=1:3  %patchy  stripy-h stripy-v
        j
        data=REV(j,subj,1);  %pull out data
        ntrial(j)=length(data.MaxRevealingTrial); %trials in dataset
        revealings{j}=data.MaxRevealingTrial; % revealings in each trial
    
        clf
plot(data.RevealPosX(:,1),data.RevealPosY(:,1),'o')
axis equal
        pause
        %save all fixations in one big array
        count=0;
        nfix=sum(data.MaxRevealingTrial);
        fix{j}.sac=zeros(nfix,ndim*ndim); %this will store all the smoothed saccades as a vector
        
        for k=1:ntrial(j) %trials
            for rev=1:data.MaxRevealingTrial(k) % loop over revealings for the trial
                count=count+1;
                
                x1=data.RevealPosX(k,rev);
                y1=data.RevealPosY(k,rev);
                tmp=ksdensity2d([x1 y1],xx,xx,[sd sd]); %smooth each fixation
                
                fix{j}.sac(count,:)=tmp(:); % save individual fixations as a vector
                fix{j}.rev(count)=rev; % revealing number
                fix{j}.trial(count)=k; % trial number
                fix{j}.weight(count)=ceil(rev/5); % weighting for sampling
            end
        end
        clear data
        
    end
    
    
    %%
    for reshuffle=0
        for rep=1:nrep  %splittings
            rep
            % split into two groups of trials randomly
            for j=1:3
                % extract randomly selected  trials but equal with regard
                % to number of revealings in each half
                n{j,1}=[];
                n{j,2}=[];
                flag=1; %this us used to altternate which half gets the odd number of trials
                for k=[5:5:25]  %make this even for each
                    s=find(revealings{j}==k);
                    
                    ns=length(s);
                    if j==1
                        ns=round(ns/2);
                    end
                    if flag
                        n1=floor(ns/2);
                    else
                        n1=ceil(ns/2);
                    end
                    if rem(ns,2)==1  %odd so switch flag
                        flag=1-flag;
                    end
                    
                    nperm= datasample(s,ns,'Replace',false);
                    n{j,1}=[n{j,1}; nperm(1:n1)];  %trials in first split
                    n{j,2}=[n{j,2}; nperm(n1+1:end)]; %trials in second split
                end
                
         
                for k=1:2
                    s = ismember(fix{j}.trial,n{j,k}); %find the selected trials for this half
                    p=fix{j}.weight .* s; %choose all the fixations from these trials with weights
                    p=fix{j}.weight; %choose all the fixations from these trials with weights

                    for rev=minrev:nrev
                        if reshuffle  %this is for reshuffling
                            ps=p; %extract simply by weights with no regard to order
                            e=datasample(1:length(p),nmin,'Replace',false,'Weights',ps+eps);
                        else
                            e=[];
                            for i=1:rev
                                ps=p.* (fix{j}.rev==i); %choose all the fixations from these trials with weights
                                e=[e datasample(1:length(p),floor(nmin/rev),'Replace',false,'Weights',ps+eps)];
                            end
                        end
                        mn{rep,j,k,rev}=mean(fix{j}.sac(e,:));
                    end
                end
            end
            
            
            
            if 0
                for j=1:3
                    mysubplot(4,3,1,j)
                    imagesc(reshape(mn{j,1,5},ndim,ndim))
                    mysubplot(4,3,2,j)
                    imagesc(reshape(mn{j,2,5},ndim,ndim))
                    mysubplot(4,3,3,j)
                    imagesc(reshape(mn_corr{j,1,5},ndim,ndim))
                    mysubplot(4,3,4,j)
                    imagesc(reshape(mn_corr{j,2,5},ndim,ndim))
                    %colorbar
                end
                pause
            end
        end
        
        gsum=zeros(nrev,ndim^2);
        
        %calculate mean across the three image types for both halves
        for rev=minrev:nrev
            for rep=1:nrep  %splitting
                for j=1:3
                    for k=1:2
                        gsum(rev,:)= gsum(rev,:)+mn{rep,j,k,rev};
                    end
                end
            end
        end
        gmean=gsum/(6*nrep);
                
        
        %calculate correlations
        c=zeros(3,2,nrev,nrep);  % to store correlations

        for rep=1:nrep  %splittings
            for j=1:3
                for k=1:3
                    for rev=minrev:nrev
                        if 0
                            c(j,k,rev,rep)=dsim(mn{rep,j,1,rev},mn{rep,k,2,rev},0.5);
                        else
                            tmp=corrcoef(mn{rep,j,1,rev}-gmean(rev,:),mn{rep,k,2,rev}-gmean(rev,:));
                            c(j,k,rev,rep)=tmp(1,2);
                        end
                    end
                end
            end
        end
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
        
        [Hp,Hs]=shadeplot(minrev:25,mean(within(minrev:end,:),2)',std(within(minrev:end,:)'),'g','g')
        set(Hs,'FaceAlpha',0.2)
        
        [Hp,Hs]=shadeplot(minrev:25,mean(between(minrev:end,:),2)',std(between(minrev:end,:)'),'r','r')
        set(Hs,'FaceAlpha',0.2)
        pause(0.1)
    end
    
    
end


%%
for k=1:8
    subplot(4,2,k)
    axis([minrev 25 -0.6 0.6])
end

%export_fig shuffle.pdf

%%
