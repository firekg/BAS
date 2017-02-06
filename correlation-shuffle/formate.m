%partcipant 1 2 3 av
%condition 1 data 2 bas

clear all
clf
load REV
count=0;

subj=1
clear sel fix ntrial data
%sets up data an sizes for input
for j=1:3  %patchy  stripy-h stripy-v
    j
    data=REV(j,subj,1);  %pull out data
    ntrial=length(data.MaxRevealingTrial); %trials in dataset
    
    %save all fixations in one big array
    
    for k=1:ntrial %trials
        for rev=1:5% loop over revealings for the trial
            count=count+1;
            D.maxrev(count)=data.MaxRevealingTrial(k);
            D.x(count)=data.RevealPosX(k,rev);
            D.y(count)=data.RevealPosY(k,rev);
            D.rev(count)=rev;
            D.trial(count)=k;
            D.image_type(count)=j;
        end
    end
   
end

save for_mate D
