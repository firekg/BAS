function [CorrBoot, CORR] = corr_DataModel_bootstrap(REV)
%% modified from SM_bootstrap
% correlate between data and model within and across parti
% SM_bootstrap includes data-data and data-model for within parti

nrev = 25;
nboot = 1000; %number of bootstrap samples
nparti = 3; % no across parti for AVG, so only 3
ncond = 9;

dataId = 1;
modelId = 2;

nmap = 77; %P and Q map size
GaussFilter = a_GaussFilter(nmap);

CorrBoot = zeros(nparti, nparti, nrev, nboot, ncond);
% CorrBoot has 5 layers:
% data parti: P1, P2, P3, AVG
% model parti: P1, P2, P3, AVG
% revealing number 25x: 1:25
% bootstrap samples 50x: 1:nboot
% condition: PA-PA, PA-SH, PA-SV, SH-PA, SH-SH, SH-SV, SV-PA, SV-SH, SV-SV

for modelParti = 1:nparti
    for dataParti = 1:nparti
        rev1PA = get_rev_position(REV(1, modelParti, modelId));
        rev1SH = get_rev_position(REV(2, modelParti, modelId));
        rev1SV = get_rev_position(REV(3, modelParti, modelId));
        rev2PA = get_rev_position(REV(1, dataParti, dataId));
        rev2SH = get_rev_position(REV(2, dataParti, dataId));
        rev2SV = get_rev_position(REV(3, dataParti, dataId));
        for irev = 1:nrev %1:nrev
            rev1PAsub = get_subrev(rev1PA, irev);
            rev1SHsub = get_subrev(rev1SH, irev);
            rev1SVsub = get_subrev(rev1SV, irev);
            rev2PAsub = get_subrev(rev2PA, irev);
            rev2SHsub = get_subrev(rev2SH, irev);
            rev2SVsub = get_subrev(rev2SV, irev);
            for iboot = 1:nboot %(can use parfor here)
                rev1PAsub_half = get_half(rev1PAsub);
                rev1SHsub_half = get_half(rev1SHsub);
                rev1SVsub_half = get_half(rev1SVsub);
                rev2PAsub_half = get_half(rev2PAsub);
                rev2SHsub_half = get_half(rev2SHsub);
                rev2SVsub_half = get_half(rev2SVsub);
                % patt maps
                map1PA = get_rev_map_mod3(rev1PAsub_half, GaussFilter, nmap);
                map1SH = get_rev_map_mod3(rev1SHsub_half, GaussFilter, nmap);
                map1SV = get_rev_map_mod3(rev1SVsub_half, GaussFilter, nmap);
                map2PA = get_rev_map_mod3(rev2PAsub_half, GaussFilter, nmap);
                map2SH = get_rev_map_mod3(rev2SHsub_half, GaussFilter, nmap);
                map2SV = get_rev_map_mod3(rev2SVsub_half, GaussFilter, nmap);
                % mean maps
                map1AVG = (map1PA + map1SH +map1SV)/3;
                map2AVG = (map2PA + map2SH +map2SV)/3;
                % mean-corrected maps
                map1PAmc = map1PA - map1AVG;
                map1SHmc = map1SH - map1AVG;
                map1SVmc = map1SV - map1AVG;
                map2PAmc = map2PA - map2AVG;
                map2SHmc = map2SH - map2AVG;
                map2SVmc = map2SV - map2AVG;
                % testing
                %subplot(1,2,1);
                %imagesc(map1SHmc);
                %subplot(1,2,2);
                %imagesc(map2SHmc);
                %corrcoef(map1SHmc(:), map2SHmc(:))
                %pause
                for icond = 1:ncond
                    if(icond==1)
                        corr1 = map1PAmc;
                        corr2 = map2PAmc;
                    elseif(icond==2) 
                        corr1 = map1PAmc;
                        corr2 = map2SHmc;
                    elseif(icond==3)
                        corr1 = map1PAmc;
                        corr2 = map2SVmc;
                    elseif(icond==4)
                        corr1 = map1SHmc;
                        corr2 = map2PAmc;
                    elseif(icond==5)
                        corr1 = map1SHmc;
                        corr2 = map2SHmc;
                    elseif(icond==6)
                        corr1 = map1SHmc;
                        corr2 = map2SVmc;
                    elseif(icond==7)
                        corr1 = map1SVmc;
                        corr2 = map2PAmc;
                    elseif(icond==8)
                        corr1 = map1SVmc;
                        corr2 = map2SHmc;
                    elseif(icond==9)
                        corr1 = map1SVmc;
                        corr2 = map2SVmc;
                    end
                    c = corrcoef(corr1(:), corr2(:)); %2-by-2
                    CorrBoot(modelParti, dataParti, irev, iboot, icond) = c(1,2);
                end
            end
        end
        fprintf('ModelParti %d, DataParti %d; \n', modelParti, dataParti);    
    end
end

% % Compute Correlation curves
% CORR has 5 layers:
CORR = nan(3,nrev,2,2,nparti);
% moments x3: (mean, upper 95%, lower 95%)
% revealing number 25x: 1:25
% parti-cond x2: within parti, across parti
% patt-cond x2: same patt, diff parti
% participants x3: SC, EP, BZ

samePattId = [1,5,9];
diffPattId = [2,3,4,6,7,8];
for parti = 1:nparti
    within_same = mean(CorrBoot(parti, parti, :, :, samePattId), 5);
    within_diff = mean(CorrBoot(parti, parti, :, :, diffPattId), 5);    
    acrossInd = [1,2,3];  acrossInd(parti) = [];
    across_same = mean(mean(CorrBoot(parti, acrossInd, :, :, samePattId), 5), 2);
    across_diff = mean(mean(CorrBoot(parti, acrossInd, :, :, diffPattId), 5), 2);
    
    within_same = sort(within_same, 4);
    within_diff = sort(within_diff, 4);
    across_same = sort(across_same, 4);
    across_diff = sort(across_diff, 4);

    CORR(1,:,1,1,parti) = mean(within_same, 4);
    CORR(2,:,1,1,parti) = within_same(1,1,:,round(nboot*0.95),1);
    CORR(3,:,1,1,parti) = within_same(1,1,:,round(nboot*0.05),1);

    CORR(1,:,1,2,parti) = mean(within_diff, 4);
    CORR(2,:,1,2,parti) = within_diff(1,1,:,round(nboot*0.95),1);
    CORR(3,:,1,2,parti) = within_diff(1,1,:,round(nboot*0.05),1);

    CORR(1,:,2,1,parti) = mean(across_same, 4);
    CORR(2,:,2,1,parti) = across_same(1,1,:,round(nboot*0.95),1);
    CORR(3,:,2,1,parti) = across_same(1,1,:,round(nboot*0.05),1);

    CORR(1,:,2,2,parti) = mean(across_diff, 4);
    CORR(2,:,2,2,parti) = across_diff(1,1,:,round(nboot*0.95),1);
    CORR(3,:,2,2,parti) = across_diff(1,1,:,round(nboot*0.05),1);

%     test = within_same(25,:,1,3,parti)<0;
%     temp = sum(test);
%     fprintf('Parti %d; Same: %f.\n', parti, temp/1000);
%     test = temp_diff(25,:,1,3,parti)>0;
%     temp = sum(test);
%     fprintf('Parti %d; Diff: %f.\n', parti, temp/1000);
end


end

%% get a random half of the locations
function [rev_half]=  get_half(rev)
    npt = numel(rev.nth);
    perm = randperm(npt);
    npt_half = round(npt/2);
    rev_half.x = rev.x(perm(1:npt_half));
    rev_half.y = rev.y(perm(1:npt_half));
    rev_half.nth = rev.nth(perm(1:npt_half));
end

%% get a subset of locations bounded by revealing number
function [revsub] = get_subrev(rev, maxRev)
    ind = rev.nth<=maxRev;
    revsub.x = rev.x(ind);
    revsub.y = rev.y(ind);
    revsub.nth = rev.nth(ind);
end

%% get revealing positinos
function [revPos] = get_rev_position(D)

ntr = D.Trials;

revx = D.RevealPosX;
revy = D.RevealPosY;
revn = repmat(1:25, ntr, 1);

revx(isnan(revx)) = 0;
revy(isnan(revy)) = 0;

revPos.x = revx(:);
revPos.y = revy(:);
revPos.nth = revn(:);

ind = revPos.x==0; %indices of 0
revPos.x(ind) =[];
revPos.y(ind) =[];
revPos.nth(ind) =[];

end

%% modified version of a_Get_RevealMap_mod2
function revMap = get_rev_map_mod3(rev, GaussFilter, nmap)

npix = 770;
ImageSize = 20;
nbin = npix/nmap;
revMap = zeros(nmap,nmap);

npos = numel(rev.x);  %bootstrap size

%compute weights
wr = (6 - ceil(rev.nth/5))/5;  %hard-wired weights

for i=1:npos    
    xDx = 1+(npix-1)*(rev.x(i)/ImageSize+0.5);
    xDy = 1-(npix-1)*(rev.y(i)/ImageSize-0.5);
    xDx = max(1,round(xDx/nbin));
    xDy = max(1,round(xDy/nbin));
    revMap(xDy,xDx) = revMap(xDy,xDx) + 1*wr(i); %x-->col, y-->row
end

revMap = conv2(revMap,GaussFilter,'same');
revMap = revMap/sum(revMap(:)); %normalize

end
