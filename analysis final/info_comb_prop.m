function [INFO, INFOC] = info_comb_prop(PROP)
% INFO has 5 layers: 
% revealing number x25 (1:25); 
% parti x4; 
% condiitons x5 (AL,PLR,BAS,mBAS,softmaxBAS); 
% moments x2 (mean,SD)
% number of prop

nprop = 20;
nsimu = 10;
% ncond = 12; %the 12th for newNoisePROP
% ncond = 7; %6 for heuristicPROP
ncond = 4; % 4 for PROP

if (nargin==0)
    INFOC = nan(25, 4, ncond, 2, nprop*nsimu);
    for isimu = 1:nsimu
%         load(['appendedPROP',num2str(isimu),'.mat']); %after performing appendToProp on newNoisePROP
%         load(['newNoisePROP_data/newNoisePROP',num2str(isimu),'.mat']);
%         load(['heuristics_data/heuristicPROP',num2str(isimu),'.mat']);
%         load(['PROP_data/PROP',num2str(isimu),'.mat']);
        load(['expPROP_data/expPROP',num2str(isimu),'.mat']);
        for iprop = 1:nprop
            [INFO] = info_prop(PROP(:,:,iprop));        
            INFOC(:,:,:,:, (isimu-1)*nprop+iprop) = INFO;
        end
    end
elseif (nargin==1)
    INFOC = nan(25, 4, ncond, 2, nprop);
    for iprop = 1:nprop
        [INFO] = info_prop(PROP(:,:,iprop));        
        INFOC(:,:,:,:,iprop) = INFO;
    end    
end

INFO = nanmean(INFOC,5);

end

%% Compute information gain
function [INFO] = info_prop(PROP)
% PROP has 2 layers: participants x3 (SC,EP,BZ); conditions x5 (ALprop, PLRprop, BASprop, mBASprop, sBASprop)
% INFO has 4 layers: revealing number x25 (1:25); parti x4; condiitons x5 (AL,PLR,BAS,mBAS,sBAS); stat x2 (mean,upper,lower)

% ncond = 12; %the 12th for newNoisePROP
% ncond = 7; %6 for heuristicPROP
ncond = 4; %4 for PROP
INFO = nan(25,4,ncond,2);
for parti = 1:3
    for cond = 1:ncond
        m = PROP(parti, cond).mla;  
        M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  
        M= -M/log(0.5);
        INFO(:, parti, cond, 1) = nanmean(M);  
        INFO(:, parti, cond, 2) = sqrt(nanvar(M)./sum(~isnan(M)));        
    end
   %disp(sprintf('parti=%d;',parti));
end
for parti=4
    for cond = 1:ncond
        m=[PROP(1, cond).mla; PROP(2, cond).mla; PROP(3, cond).mla];     
        M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  
        M= -M/log(0.5);
        INFO(:, parti, cond, 1) = nanmean(M);  
        INFO(:, parti, cond, 2) = sqrt(nanvar(M)./sum(~isnan(M)));
    end
    %disp(sprintf('parti=%d;',parti));
end
%save('INFOf.mat','INFO');

end

%% To get the slope of the information curves, I used for example:

% for linear slope
% % (1:25)'\INFO(:,4,3,1) or regress(INFO(:,4,3,1),(1:25)')
% %  note: intercept set to 0.
% slopes = zeros(3,3); % SC,EP,BZ; AL,PLR,PLB
% for ii=1:3
%     for jj=1:3
%         slopes(ii,jj)=regress(INFO(:,ii,jj,1),(1:25)');
%     end
% end

