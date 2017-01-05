%% diagnostic plot for exp performance
P = PERF;
Px = Pmix;
x = 5:5:25;
clf;
% plot PERF
for parti = 1:4
    subplot(4,4, 4*(parti-1)+1);
    hold on;
    errorbar(x,P(:,parti,2,1),P(:,parti,2,2),'-bs','MarkerFaceColor',[0,0,1],'MarkerSize',5);
    errorbar(x,P(:,parti,1,1),P(:,parti,1,2),'-ro','MarkerFaceColor',[1,0,0],'MarkerSize',5);
    %errorbar(x,P(:,parti,3,1),P(:,parti,3,2),'-kd','MarkerFaceColor',[0,0,0],'MarkerSize',5);
    %errorbar(x,P(:,parti,4,1),P(:,parti,4,2),'-d','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5],'MarkerSize',5);    
    errorbar(x,P(:,parti,5,1),P(:,parti,5,2),'-kd','MarkerFaceColor',[0,0,0],'MarkerSize',5);
    axis([0 27 0 1]);
    title(['P', num2str(parti)]);
    ylabel('P(correct)');
    xlabel('Revealing number');
end

% plot PERF_ImgVsRev
for parti = 1:4
    for img = 1:3
        for rev = 1:3
            subplot(4,4, 4*(parti-1)+1+img);
            hold on;
            cond = rev + 3*(img-1);
            % marker for rev Type
            if (rev==1) sym = '-o'; end
            if (rev==2) sym = '->'; end
            if (rev==3) sym = '-^'; end
            if (img == rev)
                cl = [0,0,0];
            else
                cl = 0.5*[1,1,1];
            end
            errorbar(x,Px(:,parti,cond,1),Px(:,parti,cond,2), sym, 'color', cl);   
            axis([0 27 0 1])
            ylabel('P(correct)');
            xlabel('Revealing number');
        end
    end
end
subplot(4,4,2);  title('Img = PA');
subplot(4,4,3);  title('Img = SH');
subplot(4,4,4);  title('Img = SV');

set(gcf,'Color',[1,1,1]);

%% diagnostic plot for choice biases
B = Bias;
Bx = BiasMix;
x = 5:5:25;
clf;
% plot PERF
for parti = 1:4
    subplot(4,4, 4*(parti-1)+1);
    hold on;
    plot(x,B(:,parti,2,1),'-bs','MarkerFaceColor',[0,0,1],'MarkerSize',5);
    plot(x,B(:,parti,1,1),'-ro','MarkerFaceColor',[1,0,0],'MarkerSize',5);
    plot(x,B(:,parti,3,1),'-kd','MarkerFaceColor',[0,0,0],'MarkerSize',5);
    plot(x,B(:,parti,4,1),'-d','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5],'MarkerSize',5);    
%     plot(x,B(:,parti,2,2),'--bs','MarkerFaceColor',[0,0,1],'MarkerSize',5);
%     plot(x,B(:,parti,1,2),'--ro','MarkerFaceColor',[1,0,0],'MarkerSize',5);
%     plot(x,B(:,parti,3,2),'--kd','MarkerFaceColor',[0,0,0],'MarkerSize',5);
%     plot(x,B(:,parti,4,2),'--d','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5],'MarkerSize',5);    
    axis([0 27 0 1]);
    title(['P', num2str(parti)]);
    ylabel('P(choice = PA)');
    xlabel('Revealing number');
end

% plot PERF_ImgVsRev
for parti = 1:4
    for img = 1:3
        for rev = 1:3
            subplot(4,4, 4*(parti-1)+1+img);
            hold on;
            cond = rev + 3*(img-1);
            % marker for rev Type
            if (rev==1) sym = 'o'; end
            if (rev==2) sym = '>'; end
            if (rev==3) sym = '^'; end
            if (img == rev)
                cl = [0,0,0];
            else
                cl = 0.5*[1,1,1];
            end
            plot(x,Bx(:,parti,cond,1), [sym, '-'], 'color', cl);   
            %plot(x,Bx(:,parti,cond,2), [sym, '--'], 'color', cl);   
            axis([0 27 0 1])
            ylabel('P(choice = PA)');
            xlabel('Revealing number');
        end
    end
end
subplot(4,4,2);  title('Img = PA');
subplot(4,4,3);  title('Img = SH');
subplot(4,4,4);  title('Img = SV');

set(gcf,'Color',[1,1,1]);

%% diagnostic plot for model performance
P = mPERF;
Px = mPmix;
% P = lpxPERF;
% Px = lpxPmix;
x = 1:25;
clf;
% plot PERF
for parti = 1:4
    subplot(4,4, 4*(parti-1)+1);
    hold on;
    errorbar(x,P(:,parti,2,1),P(:,parti,2,2),'-b','MarkerFaceColor',[0,0,1]);
    errorbar(x,P(:,parti,1,1),P(:,parti,1,2),'-r','MarkerFaceColor',[1,0,0]);
    errorbar(x,P(:,parti,3,1),P(:,parti,3,2),'-k','MarkerFaceColor',[0,0,0],'LineWidth',2);
    errorbar(x,P(:,parti,4,1),P(:,parti,4,2),'-','Color',[0.5,0.5,0.5],'LineWidth',2);    
    %errorbar(x,P(:,parti,5,1),P(:,parti,5,2),'-k','MarkerFaceColor',[0,0,0],'LineWidth',2);
    axis([0 27 0 1]);
    title(['P', num2str(parti)]);
    ylabel('P(correct)');
    xlabel('Revealing number');
end

% plot PERF_ImgVsRev
for parti = 1:4
    for img = 1:3
        for rev = 1:3
            subplot(4,4, 4*(parti-1)+1+img);
            hold on;
            cond = rev + 3*(img-1);
            % marker for rev Type
            if (rev==1) sym = '-o'; end
            if (rev==2) sym = '->'; end
            if (rev==3) sym = '-^'; end
            if (img == rev)
                cl = [0,0,0];
            else
                cl = 0.5*[1,1,1];
            end
            errorbar(x,Px(:,parti,cond,1),Px(:,parti,cond,2), sym, 'color', cl);   
            axis([0 27 0 1])
            ylabel('P(correct)');
            xlabel('Revealing number');
        end
    end
end
subplot(4,4,2);  title('Img = PA');
subplot(4,4,3);  title('Img = SH');
subplot(4,4,4);  title('Img = SV');

set(gcf,'Color',[1,1,1]);

%% plot only average perf both human and model -- this is very much like plot_mperf

P = PERF;
x = 5:5:25;

mP = mPERF;
mx = 1:25;

clf;
% plot PERF
for parti = 1:4
    subplot(1, 4, parti);
    hold on;
    errorbar(x,P(:,parti,2,1),P(:,parti,2,2),'bs','MarkerFaceColor',[0,0,1]);
    errorbar(x,P(:,parti,1,1),P(:,parti,1,2),'ro','MarkerFaceColor',[1,0,0]);
    errorbar(x,P(:,parti,3,1),P(:,parti,3,2),'kd','MarkerFaceColor',[0,0,0]);
    errorbar(x,P(:,parti,4,1),P(:,parti,4,2),'^','Color',0.5*[1,1,1],'MarkerFaceColor',0.5*[1,1,1]);    
    plot(mx,mP(:,parti,2,1),'b-');
    plot(mx,mP(:,parti,1,1),'r-');
    plot(mx,mP(:,parti,3,1),'k-');
    plot(mx,mP(:,parti,4,1),'-','Color',0.5*[1,1,1]);    
    axis([0 27 0.4 1]);
    title(['P', num2str(parti)]);
    ylabel('P(correct)');
    xlabel('Revealing number');
    axis square
end

set(gcf,'Color',[1,1,1]);
