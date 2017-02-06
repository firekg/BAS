%4 - image type patch, sh sv, average
%participant 1 2 3 av
%condition 1 data 2 bas

for subj=1:3
    
    
    str=sprintf('load subj%i ',subj);
    eval(str);
    
    
    
    % plot
    mysubplot(4,2,subj,1)
    hold on
    
    [Hp,Hs]=shadeplot(minrev:25,mean(within(minrev:end,:)'),std(within(minrev:end,:)'),'g','g')
    set(Hs,'FaceAlpha',0.2)
    
    [Hp,Hs]=shadeplot(minrev:25,mean(between(minrev:end,:)'),std(between(minrev:end,:)'),'r','r')
    set(Hs,'FaceAlpha',0.2)
    pause(0.1)
end


for subj=1:4
    
    
    str=sprintf('load subj%ibas ',subj);
    eval(str);
    
    
    
    % plot
    mysubplot(4,2,subj,2)
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
    axis([minrev 25 -0.6 0.6])
end


%export_fig corr.pdf

