load('expDataWPix.mat');

%%
opengl software

for cond = 1:3
    if(cond==1)
        condName = 'AL';
    elseif(cond==2)  
        condName = 'PLR';
    elseif(cond==3)  
        condName = 'PLB';
    end
    for parti=1:4
        if(parti==1)
            partiName = 'SC';
        elseif(parti==2)  
            partiName = 'EP';
        elseif(parti==3)  
            partiName = 'BZ';
        elseif(parti==4)
            partiName = 'AVG';
        end

        if(parti<4)
            D = eval([partiName,condName,'wPix']);
        else
            D = a_Combine_Data( eval(['SC',condName,'wPix']),...
                                eval(['EP',condName,'wPix']),...
                                eval(['BZ',condName,'wPix']));
        end
        
        revMax = D.MaxRevealingTrial;
        correct = D.AnswerChoice == D.AnswerReal;
        decTime = D.StateAnswerTime/1000;
        
        fprintf(['parti %d ', condName, ', avg time: %f\n'], parti, nanmean(decTime));

        % Correct/Incorrect times as a function of revealing number
%         correctTime = zeros(5,2); %mean and std
%         mistakeTime = zeros(5,2);
%         nrev = [5,10,15,20,25];
%         for ind = 1:5
%             rev = nrev(ind);
%             times = decTime( correct & revMax==rev );
%             correctTime(ind,1) = nanmean(times);
%             correctTime(ind,2) = sqrt(nanvar(times)/numel(times));
%             times = decTime( ~correct & revMax==rev );
%             mistakeTime(ind,1) = nanmean(times);
%             mistakeTime(ind,2) = sqrt(nanvar(times)/numel(times));
%         end
%         
%         subplot(3,4, (cond-1)*4 + parti);
%         hold on;
%         cl=[1.0,0.0,0.0];  a_ErrorShade(nrev, correctTime(:,1)-correctTime(:,2),correctTime(:,1)+correctTime(:,2),cl,0.5);    
%         cl=[0.0,0.0,1.0];  a_ErrorShade(nrev, mistakeTime(:,1)-mistakeTime(:,2),mistakeTime(:,1)+mistakeTime(:,2),cl,0.5);    
%         plot(nrev, correctTime(:,1), 'r-o');
%         plot(nrev, mistakeTime(:,1), 'b-s');
% 
%         axis([5,25,0,25]);
%         xlabel('Revealing number');
%         ylabel('Decision time (s)');
%         title([partiName,'-',condName]);
% 
%         if (parti==4)
%            legend('correct','incorrect');
%         end

        % Answer correctness as a function of decTime
        [time, order] = sort(decTime);
        nTr = numel(decTime);
        correct = correct(order);
        nbin = 10;
        tbin = nan(nbin,1); %mean, SE
        right = nan(nbin,3); %mean, low, upp
        for ibin = 1:nbin
            ind1 = (ibin-1)*floor(nTr/nbin) + 1;
            if (ibin == nbin)
                ind2 = nTr;
            else
                ind2 =  ibin*floor(nTr/nbin);
            end
            tbin(ibin) = mean(time(ind1:ind2));
            [phat,pci] = binofit(sum(correct(ind1:ind2)), ind2-ind1+1, 0.05);
            right(ibin,1) = phat;
            right(ibin,2) = abs(pci(1)-phat);
            right(ibin,3) = abs(pci(2)-phat);
        end
        
        p = polyfit(time, correct, 1);
        x1 = linspace(0, tbin(end), 100);
        y1 = polyval(p,x1);
        
        subplot(3,4, (cond-1)*4 + parti);
        hold on;
        errorbar(tbin, right(:,1), right(:,2), right(:,3), 'o');
        plot(x1, y1, 'k-');
        axis([0,inf,0,1]);
        ylabel('Fraction correct');
        xlabel('Decision time (s)');
        title([partiName,'-',condName]);
    end
end

set(gcf,'Color',[1,1,1]);
