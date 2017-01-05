function [CORR] = corr_comb_SMboot(RdrevMap, SMboot)
%% Combine SIM and PROP data from fear
% initialize

if (nargin==0)
    %load('SMboot_data/SMboot1.mat');
    load('newNoisePROP_data/newNoiseSMboot1.mat');
    SMbootC = SMboot; % combined SMboot
    % start combining data
    for task=2:20
        %load(['SMboot_data/SMboot',num2str(task),'.mat']);
        load(['newNoisePROP_data/newNoiseSMboot',num2str(task),'.mat']);
        SMbootC = [SMbootC, SMboot];
    end
    %save('SMbootC.mat','SMbootC');
elseif (nargin==2)
    %SMbootC = SMboot;
    % 2015-09-07    
    SMbootC = SMboot(1).SMboot;
    for task = 2:10
        SMbootC = [SMbootC, SMboot(task).SMboot];
    end
end

% % Compute Correlation curves
% CORR has 5 layers:
% revealing number 25x: 1:25
% moments x3: (mean, upper 95%, lower 95%)
% pattern (mis)match 2x: same, different
% conditions 3x: self-self, self-BAS, self-nBAS
% participants 4x: SC, EP, BZ, AVG
s =size(SMbootC);
nboot = s(2);  %hard-coded
temp_same = mean(SMbootC(:,:,1:3,:,:),3);  temp_same = sort(temp_same,2);
temp_diff = mean(SMbootC(:,:,4:9,:,:),3);  temp_diff = sort(temp_diff,2);
CORR(:,1,1,:,:) = mean(temp_same,2);
CORR(:,2,1,:,:) = temp_same(:,round(nboot*0.95),1,:,:);
CORR(:,3,1,:,:) = temp_same(:,round(nboot*0.05),1,:,:);
CORR(:,1,2,:,:) = mean(temp_diff,2);
CORR(:,2,2,:,:) = temp_diff(:,round(nboot*0.95),1,:,:);
CORR(:,3,2,:,:) = temp_diff(:,round(nboot*0.05),1,:,:);

% % do direct map correlation for self-nBAS
% %RdrevMap has 6 layers: 770; 770; pattern; parti; cond; rev number
% rM = zeros(2,2,9);
% rLOM = rM;
% rUPM = rM;
% pM = rM;
% for rev = 1:25
%     for parti = 1:4
%         for zz = 1:9
%             if(zz==1)
%                 xx = RdrevMap(:,:,1,parti,1,rev);
%                 yy = RdrevMap(:,:,1,parti,3,rev);                
%             elseif(zz==2)
%                 xx = RdrevMap(:,:,2,parti,1,rev);
%                 yy = RdrevMap(:,:,2,parti,3,rev);                
%             elseif(zz==3)
%                 xx = RdrevMap(:,:,3,parti,1,rev);
%                 yy = RdrevMap(:,:,3,parti,3,rev);                
%             elseif(zz==4)
%                 xx = RdrevMap(:,:,1,parti,1,rev);
%                 yy = RdrevMap(:,:,2,parti,3,rev);                
%             elseif(zz==5)
%                 xx = RdrevMap(:,:,1,parti,1,rev);
%                 yy = RdrevMap(:,:,3,parti,3,rev);                
%             elseif(zz==6)
%                 xx = RdrevMap(:,:,2,parti,1,rev);
%                 yy = RdrevMap(:,:,1,parti,3,rev);                
%             elseif(zz==7)
%                 xx = RdrevMap(:,:,2,parti,1,rev);
%                 yy = RdrevMap(:,:,3,parti,3,rev);                
%             elseif(zz==8)
%                 xx = RdrevMap(:,:,3,parti,1,rev);
%                 yy = RdrevMap(:,:,1,parti,3,rev);                
%             elseif(zz==9)
%                 xx = RdrevMap(:,:,3,parti,1,rev);
%                 yy = RdrevMap(:,:,2,parti,3,rev);                
%             end
%             [r,p,rLO,rUP] = corrcoef(xx(:),yy(:));
%             rM(:,:,zz) = r;
%             pM(:,:,zz) = p;
%             rLOM(:,:,zz) = rLO;
%             rUPM(:,:,zz) = rUP;
%         end
%         % within-type
%         rtemp = mean(rM(:,:,1:3),3);
%         ptemp = mean(pM(:,:,1:3),3);
%         rLOtemp = mean(rLOM(:,:,1:3),3);
%         rUPtemp = mean(rUPM(:,:,1:3),3);
%         CORR(rev,1,1,3,parti) = rtemp(1,2); %mean
%         CORR(rev,2,1,3,parti) = rUPtemp(1,2); %upper
%         CORR(rev,3,1,3,parti) = rLOtemp(1,2); %lower
%         % across-type
%         rtemp = mean(rM(:,:,4:9),3);
%         ptemp = mean(pM(:,:,4:9),3);
%         rLOtemp = mean(rLOM(:,:,4:9),3);
%         rUPtemp = mean(rUPM(:,:,4:9),3);
%         CORR(rev,1,2,3,parti) = rtemp(1,2); %mean
%         CORR(rev,2,2,3,parti) = rUPtemp(1,2); %upper
%         CORR(rev,3,2,3,parti) = rLOtemp(1,2); %lower
%     end
% end

% SMboot has 5 layers:
% revealing number 25x: 1:25
% bootstrap samples bootx: 1:nboot
% pattern (mis)match 9x: PA-PA, SH-SH, SV-SV; PA-SH, PA-SV, SH-PA, SH-SV, SV-PA, SV-SH
% conditions 3x: self-self, self-BAS, self-nBAS
% participants 4x: SC, EP, BZ, AVG
% test zero-cross for CORR(25,2,2,1,1)
for parti=1:4
    test = temp_same(25,:,1,2,parti)<0;
    temp = sum(test);
    disp(sprintf('Parti %d; Same: %f.', parti, temp/1000));
    test = temp_diff(25,:,1,2,parti)>0;
    temp = sum(test);
    disp(sprintf('Parti %d; Diff: %f.', parti, temp/1000));
end

%save('CORR.mat','CORR');

end

%% 2015-09-22
% P(corr_correct > corr_incorrect)
% 
% SMbootC = CombSMcorrect(1).SMboot;
% for task = 2:10
%     SMbootC = [SMbootC, CombSMcorrect(task).SMboot];
% end
% correct_same = mean(SMbootC(25,:,1:3,2,4),3);
% 
% SMbootC = CombSMincorrect(1).SMboot;
% for task = 2:10
%     SMbootC = [SMbootC, CombSMincorrect(task).SMboot];
% end
% incorrect_same = mean(SMbootC(25,:,1:3,2,4),3);
% 
% p = sum(correct_same > incorrect_same)/1000;
% 
% fprintf('P(rho_correct > rho_incorrect) = %f\n', p);




