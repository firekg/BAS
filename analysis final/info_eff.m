function wbl_ab = info_eff(INFOC)
% %% fit the curves (won't work with INFOC)
% xx = 1:25;
% for parti = 1:4
%     %AL
%     info = INFO(:,parti,1,1);
%     f = @(par) sum(( wblcdf(xx',par(1),par(2)) - info(:) ).^2);
%     optpar = fminsearch(f, [15,1]);
%     wbl_a(1,parti) = optpar(1);
%     wbl_b(1,parti) = optpar(2);
% 
%     %Random
%     info = INFO(:,parti,2,1);
%     f = @(par) sum(( wblcdf(xx',par(1),par(2)) - info(:) ).^2);
%     optpar = fminsearch(f, [15,1]);
%     wbl_a(2,parti) = optpar(1);
%     wbl_b(2,parti) = optpar(2);
% 
%     %BAS
%     info = INFO(:,parti,3,1); 
%     f = @(par) sum(( wblcdf(xx',par(1),par(2)) - info(:) ).^2);
%     optpar = fminsearch(f, [15,1]);
%     wbl_a(3,parti) = optpar(1);
%     wbl_b(3,parti) = optpar(2);
%     
%     %noisy BAS
%     info = INFO(:,parti,5,1); 
%     f = @(par) sum(( wblcdf(xx',par(1),par(2)) - info(:) ).^2);
%     optpar = fminsearch(f, [15,1]);
%     wbl_a(4,parti) = optpar(1);
%     wbl_b(4,parti) = optpar(2);    
% end

%% joint fits
nboot = 200;
ncond = 7;
nparti = 4;
wbl_ab = zeros(ncond+1, nparti, nboot);

xx = 1:25;
for iboot = 1:nboot
    for parti = 1:nparti
        %AL
        infoAL = INFOC(:,parti,1,1,iboot);
        infoRan = INFOC(:,parti,2,1,iboot);
        infoBAS = INFOC(:,parti,3,1,iboot); 
        infoNBAS = INFOC(:,parti,5,1,iboot);
        infoHeu1 = INFOC(:,parti,7,1,iboot);
        infoHeu2 = INFOC(:,parti,8,1,iboot);
        infoHeu3 = INFOC(:,parti,9,1,iboot);
        f = @(par) sum( ( wblcdf(xx',par(2),par(1)) - infoAL(:) ).^2 ...
                    + ( wblcdf(xx',par(3),par(1)) - infoRan(:) ).^2 ...
                    + ( wblcdf(xx',par(4),par(1)) - infoBAS(:) ).^2 ...
                    + ( wblcdf(xx',par(5),par(1)) - infoNBAS(:) ).^2 ...
                    + ( wblcdf(xx',par(6),par(1)) - infoHeu1(:) ).^2 ...
                    + ( wblcdf(xx',par(7),par(1)) - infoHeu2(:) ).^2 ...
                    + ( wblcdf(xx',par(8),par(1)) - infoHeu3(:) ).^2 );
        options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000);        
        optpar = fminsearch(f, [1, 15,15,15,15, 15,15,15], options);
        wbl_ab(:,parti,iboot) = optpar;
    end
    fprintf('boot %d done.\n', iboot);
end

%% report efficiencies

% % BAS wrt AL
% effv = wbl_ab(1+1,:,:)./wbl_ab(3+1,:,:);
% effv = sort(effv,3);
% eff_avg = mean(effv,3);
% eff_low = effv(:,:,round(nboot*0.05));
% eff_high = effv(:,:,round(nboot*0.95));
% for parti=1:3
%    fprintf('parti %d BAS/AL: %f(%f,%f)\n', parti, eff_avg(parti), eff_low(parti), eff_high(parti));  
% end
% fprintf('parti avg BAS/AL: %f(%f,%f)\n', mean(eff_avg), mean(eff_low), mean(eff_high));  
% 
% 
% % nBAS wrt AL
% effv = wbl_ab(1+1,:,:)./wbl_ab(4+1,:,:);
% effv = sort(effv,3);
% eff_avg = mean(effv,3);
% eff_low = effv(:,:,round(nboot*0.05));
% eff_high = effv(:,:,round(nboot*0.95));
% for parti=1:4
%    fprintf('parti %d nBAS/AL: %f(%f,%f)\n', parti, eff_avg(parti), eff_low(parti), eff_high(parti));  
% end
% fprintf('parti avg nBAS/AL: %f(%f,%f)\n', mean(eff_avg), mean(eff_low), mean(eff_high));  

% AL wrt random
effv = wbl_ab(2+1,:,:)./wbl_ab(1+1,:,:);
effv = sort(effv,3);
eff_avg = mean(effv,3);
eff_low = effv(:,:,round(nboot*0.05));
eff_high = effv(:,:,round(nboot*0.95));
for parti=1:3
   fprintf('parti %d AL/random: %f(%f,%f)\n', parti, eff_avg(parti), eff_low(parti), eff_high(parti));  
end
fprintf('parti avg AL/random: %f(%f,%f)\n', mean(eff_avg), mean(eff_low), mean(eff_high));  

% % AL wrt Heu1
% effv = wbl_ab(5+1,:,:)./wbl_ab(1+1,:,:);
% effv = sort(effv,3);
% eff_avg = mean(effv,3);
% eff_low = effv(:,:,round(nboot*0.05));
% eff_high = effv(:,:,round(nboot*0.95));
% for parti=1:3
%    fprintf('parti %d AL/Heu1: %f(%f,%f)\n', parti, eff_avg(parti), eff_low(parti), eff_high(parti));  
% end
% fprintf('parti avg AL/Heu1: %f(%f,%f)\n', mean(eff_avg), mean(eff_low), mean(eff_high));  
% 
% % AL wrt Heu2
% effv = wbl_ab(6+1,:,:)./wbl_ab(1+1,:,:);
% effv = sort(effv,3);
% eff_avg = mean(effv,3);
% eff_low = effv(:,:,round(nboot*0.05));
% eff_high = effv(:,:,round(nboot*0.95));
% for parti=1:3
%    fprintf('parti %d AL/Heu2: %f(%f,%f)\n', parti, eff_avg(parti), eff_low(parti), eff_high(parti));  
% end
% fprintf('parti avg AL/Heu2: %f(%f,%f)\n', mean(eff_avg), mean(eff_low), mean(eff_high));  
% 
% % AL wrt Heu3
% effv = wbl_ab(7+1,:,:)./wbl_ab(1+1,:,:);
% effv = sort(effv,3);
% eff_avg = mean(effv,3);
% eff_low = effv(:,:,round(nboot*0.05));
% eff_high = effv(:,:,round(nboot*0.95));
% for parti=1:3
%    fprintf('parti %d AL/Heu3: %f(%f,%f)\n', parti, eff_avg(parti), eff_low(parti), eff_high(parti));  
% end
% fprintf('parti avg AL/Heu3: %f(%f,%f)\n', mean(eff_avg), mean(eff_low), mean(eff_high));  
% 
% % AL wrt nBAS
% effv = wbl_ab(4+1,:,:)./wbl_ab(1+1,:,:);
% effv = sort(effv,3);
% eff_avg = mean(effv,3);
% eff_low = effv(:,:,round(nboot*0.05));
% eff_high = effv(:,:,round(nboot*0.95));
% for parti=1:3
%    fprintf('parti %d AL/nBAS: %f(%f,%f)\n', parti, eff_avg(parti), eff_low(parti), eff_high(parti));  
% end
% fprintf('parti avg AL/nBAS: %f(%f,%f)\n', mean(eff_avg), mean(eff_low), mean(eff_high));  
% 
% % random wrt nBAS
% effv = wbl_ab(4+1,:,:)./wbl_ab(2+1,:,:); 
% effv = sort(effv,3);
% eff_avg = mean(effv,3);
% eff_low = effv(:,:,round(nboot*0.05));
% eff_high = effv(:,:,round(nboot*0.95));
% for parti=1:3
%    fprintf('parti %d Rand/nBAS: %f(%f,%f)\n', parti, eff_avg(parti), eff_low(parti), eff_high(parti));  
% end
% fprintf('parti avg Rand/nBAS: %f(%f,%f)\n', mean(eff_avg), mean(eff_low), mean(eff_high));  




end

%% plot fits
% rev = 1:0.01:150;
% xx = 1:25;
% nparti = 4;
% figure(1)
% for parti = 1:4
%     subplot(1,nparti,parti);
%     plot(xx, INFO(:,parti,1,1), 'ro'); hold on;
%     plot(xx, INFO(:,parti,2,1), 'bs');
%     plot(xx, INFO(:,parti,3,1), 'kd');
%     plot(xx, INFO(:,parti,5,1), 'go');
%     % fit with no shared parameters
%     plot(rev, wblcdf(rev, wbl_a(1,parti), wbl_b(1,parti)), 'r-');
%     plot(rev, wblcdf(rev, wbl_a(2,parti), wbl_b(2,parti)), 'b-');
%     plot(rev, wblcdf(rev, wbl_a(3,parti), wbl_b(3,parti)), 'k-');
%     plot(rev, wblcdf(rev, wbl_a(4,parti), wbl_b(4,parti)), 'g-');
%     % fit with shared parameters
%     plot(rev, wblcdf(rev, wbl_ab(1,parti), wbl_ab(5,parti)), 'r:');
%     plot(rev, wblcdf(rev, wbl_ab(2,parti), wbl_ab(5,parti)), 'b:');
%     plot(rev, wblcdf(rev, wbl_ab(3,parti), wbl_ab(5,parti)), 'k:');
%     plot(rev, wblcdf(rev, wbl_ab(4,parti), wbl_ab(5,parti)), 'g:');
% end

%% find how many revealings to get to some bits
% npt = 100;
% rcond = 4;
% prob = linspace(0.6, 0.95, npt);
% eff = zeros(4,4); %in terms of revealing number
% effRatio = zeros(npt,nparti,rcond);
% eff_jf = zeros(4,4);
% effRatio_jf = zeros(npt,nparti,rcond);
% 
% % fit with no shared parameters
% for pt = 1:npt
%     target = 1+prob(pt).*log2(prob(pt))+(1-prob(pt)).*log2(1-prob(pt));
%     for parti = 1:4
%         for cond = 1:4
%             info_wbl = wblcdf(rev, wbl_a(cond,parti), wbl_b(cond,parti));
%             [val, ind] = min(abs(info_wbl - target));
%             eff(cond, parti) = rev(ind);
%         end
%     end
%     effRatio(pt,:,1) = eff(2,:)./eff(1,:); %AL/rand
%     effRatio(pt,:,2) = eff(3,:)./eff(1,:); %AL/BAS
%     effRatio(pt,:,3) = eff(4,:)./eff(1,:); %AL/nBAS    
%     effRatio(pt,:,4) = eff(4,:)./eff(2,:); %rand/nBAS    
% end
% 
% % joint fit
% for pt = 1:npt
%     target = 1+prob(pt).*log2(prob(pt))+(1-prob(pt)).*log2(1-prob(pt));
%     for parti = 1:4
%         for cond = 1:4
%             info_wbl = wblcdf(rev, wbl_ab(cond,parti), wbl_ab(5,parti));
%             [val, ind] = min(abs(info_wbl - target));
%             eff_jf(cond, parti) = rev(ind);
%         end
%     end
%     effRatio_jf(pt,:,1) = eff_jf(2,:)./eff_jf(1,:); %AL/rand
%     effRatio_jf(pt,:,2) = eff_jf(3,:)./eff_jf(1,:); %AL/BAS
%     effRatio_jf(pt,:,3) = eff_jf(4,:)./eff_jf(1,:); %AL/nBAS    
%     effRatio_jf(pt,:,4) = eff_jf(4,:)./eff_jf(2,:); %rand/nBAS    
% end
% 
% figure(2);
% for parti = 1:nparti
%     subplot(1,nparti,parti);
%     % indep fit
%     plot(prob,effRatio(:,parti,1), 'b-', 'LineWidth', 2); %AL/rand
%     hold on;
%     plot(prob,effRatio(:,parti,2), 'k-', 'LineWidth', 2); %AL/BAS
%     plot(prob,effRatio(:,parti,3), 'r-', 'LineWidth', 2); %AL/nBAS
%     plot(prob,effRatio(:,parti,4), 'c-', 'LineWidth', 2); %rand/nBAS
%     % joint fit results
%     plot(prob,effRatio_jf(:,parti,1), 'b:', 'LineWidth', 2); %AL/rand
%     plot(prob,effRatio_jf(:,parti,2), 'k:', 'LineWidth', 2); %AL/BAS
%     plot(prob,effRatio_jf(:,parti,3), 'r:', 'LineWidth', 2); %AL/nBAS
%     plot(prob,effRatio_jf(:,parti,4), 'c:', 'LineWidth', 2); %rand/nBAS
%     axis([0.6, 1, 0 ,4]);
% end
% legend('rand/AL','BAS/AL','nBAS/AL', 'nBAS/rand');

%% get a sense of entropy
% prob = 0.75;
% 1+prob.*log2(prob)+(1-prob).*log2(1-prob)
% --> p=0.75: bits=0.1887
% plot(p,y)

%% get a sense of Weibull
% yy = wblcdf(xx,10,5);
% plot(xx,yy)

%% 2015-06-25: an proposal for computing efficiency, not used
% function INFOR = infoRatio(INFOC)
% 
% [nrev, nparti, ncond, ~, nboot] = size(INFOC);
% INFOR = zeros(nrev, nparti, ncond, ncond, nboot);
% for rowcond = 1:ncond
%     for colcond= 1:ncond
%         info1 = INFOC(:,:,rowcond,1,:);
%         info2 = INFOC(:,:,colcond,1,:);
%         nor = INFOC(:,:,3,1,:);
%         base = INFOC(:,:,2,1,:);
%         INFOR(:,:,rowcond,colcond,:) = ((info1-base) - (info2-base))./(nor-base);
%         %INFOR(:,:,rowcond,colcond,:) = (info1-base)./(info2-base);
%         %INFOR(:,:,rowcond,colcond,:) = info1./info2;
%     end
% end
% 
% % conditions: AL, PLR, PLB, PLaB, nBAS, pBAS, heu1, heu2, heu3;
% 
% % ratios has 3 layers: nrev; cond; mean (lowCI + highCI)
% ncondr = 7;
% ratios = zeros(nrev, ncondr, 3);
% 
% %cond: 'AL/PLR', 'BAS/AL', 'nBAS/AL', 'AL/pi-od', 'AL/pd-oi', 'AL/pd-od'
% condr1 = [1,3,5,1,1,1,3];
% condr2 = [2,1,1,7,8,9,2];
% 
% for icondr = 1:ncondr
%      temp = INFOR(:,4,condr1(icondr),condr2(icondr),:);
%      temp = sort(temp,5);
%      ratios(:,icondr,1) = mean(temp,5);
%      ratios(:,icondr,2) = temp(:,1,1,1,round(nboot*0.05));
%      ratios(:,icondr,3) = temp(:,1,1,1,round(nboot*0.95));     
% end
% 
% specr(1,:) = 'r-o';
% specr(2,:) = 'k-o';
% specr(3,:) = 'g-o';
% specr(4,:) = 'r-s';
% specr(5,:) = 'r-^';
% specr(6,:) = 'r->';
% specr(7,:) = 'k-s';
% revX = 1:25;
% 
% for icondr = 1:ncondr
%     plot(revX, ratios(:,icondr,1), specr(icondr,:)); hold on;
%     %plot(revX, ratios(:,icondr,2), 'c:');
%     %plot(revX, ratios(:,icondr,3), 'm:');
% end
% axis([0,25,0,1]);
% %legend('AL/PLR', 'BAS/AL', 'nBAS/AL', 'AL/pi-od', 'AL/pd-oi', 'AL/pd-od');
% 
% cutoff = 15;
% fprintf('AL/PLR: %f\n', mean(ratios(cutoff:nrev,1)));
% fprintf('BAS/AL: %f\n', mean(ratios(cutoff:nrev,2)));
% fprintf('nBAS/AL: %f\n', mean(ratios(cutoff:nrev,3)));
% fprintf('AL/pi-od: %f\n', mean(ratios(cutoff:nrev,4)));
% fprintf('AL/pd-oi: %f\n', mean(ratios(cutoff:nrev,5)));
% fprintf('AL/pd-od: %f\n', mean(ratios(cutoff:nrev,6)));
% 
% end







