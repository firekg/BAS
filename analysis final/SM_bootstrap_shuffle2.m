function [SMboot] = SM_bootstrap_shuffle2(REV)
%% modified from SM_bootstrap_shuffle - 2015-09-24
% REV has 3 layers: 
% image type (PA, SH, SV, AVG)
% parti (P1, P2, P3, AVG)
% condition (data, BAS)

nrev = 25;
nboot = 10; %number of bootstrap samples
npatt = 9;
ncond = 2; %2 if no noisy BAS
nparti = 4;

SMboot = zeros(nrev, nboot, npatt, ncond, nparti);
% SMboot has 5 layers:
% revealing number 25x: 1:25
% bootstrap samples bootx: 1:nboot
% pattern (mis)match 9x: PA-PA, SH-SH, SV-SV; PA-SH, PA-SV, SH-PA, SH-SV, SV-PA, SV-SH
% conditions 3x: self-self, self-BAS, self-nBAS
% participants 4x: SC, EP, BZ, AVG

nmap = 77; %P and Q map size
GaussFilter = a_GaussFilter(nmap);

for parti = 1:nparti
    %avoid correlation at 1st revealing
    P1=REV(1,parti,1);  P2=REV(2,parti,1);  P3=REV(3,parti,1);
    [P1x, P1y, P1n] = get_rev_position(P1);
    [P2x, P2y, P2n] = get_rev_position(P2);
    [P3x, P3y, P3n] = get_rev_position(P3);
    [P1x, P1y, P1n] = equilizeSampNum(P1x, P1y, P1n);
    [P2x, P2y, P2n] = equilizeSampNum(P2x, P2y, P2n);
    [P3x, P3y, P3n] = equilizeSampNum(P3x, P3y, P3n);
    for zz = 1:npatt
        for iboot = 1:nboot            
            for irev = 1:nrev
                
                % permute trials, divide in two halves:
                % pattern 1
                nnr = numel(P1n);  perm = randperm(nnr);  nnrh = round(nnr/2);
                P1nh1 = P1n(perm(1:nnrh));  P1nh2 = P1n(perm(nnrh+1:nnr));
                P1xh1 = P1x(perm(1:nnrh));  P1xh2 = P1x(perm(nnrh+1:nnr));
                P1yh1 = P1y(perm(1:nnrh));  P1yh2 = P1y(perm(nnrh+1:nnr));
                % pattern 2
                nnr = numel(P2n);  perm = randperm(nnr);  nnrh = round(nnr/2);
                P2nh1 = P2n(perm(1:nnrh));  P2nh2 = P2n(perm(nnrh+1:nnr));
                P2xh1 = P2x(perm(1:nnrh));  P2xh2 = P2x(perm(nnrh+1:nnr));
                P2yh1 = P2y(perm(1:nnrh));  P2yh2 = P2y(perm(nnrh+1:nnr));
                % pattern 3
                nnr = numel(P3n);  perm = randperm(nnr);  nnrh = round(nnr/2);
                P3nh1 = P3n(perm(1:nnrh));  P3nh2 = P3n(perm(nnrh+1:nnr));
                P3xh1 = P3x(perm(1:nnrh));  P3xh2 = P3x(perm(nnrh+1:nnr));
                P3yh1 = P3y(perm(1:nnrh));  P3yh2 = P3y(perm(nnrh+1:nnr));
                
                % permute then restrict ind on first half (shuffling step)
                nnr = numel(P1nh1);  perm = randperm(nnr);  P1nh1 = P1nh1(perm);
                nnr = numel(P2nh1);  perm = randperm(nnr);  P2nh1 = P2nh1(perm);
                nnr = numel(P3nh1);  perm = randperm(nnr);  P3nh1 = P3nh1(perm);
                ind = P1nh1<=irev;  
                xQ1 = P1xh1(ind);  yQ1 = P1yh1(ind);  nQ1 = P1nh1(ind);
                ind = P2nh1<=irev;  
                xQ2 = P2xh1(ind);  yQ2 = P2yh1(ind);  nQ2 = P2nh1(ind);
                ind = P3nh1<=irev;  
                xQ3 = P3xh1(ind);  yQ3 = P3yh1(ind);  nQ3 = P3nh1(ind);

                % permute then restrict ind on second half (shuffling step)
                nnr = numel(P1nh2);  perm = randperm(nnr);  P1nh2 = P1nh2(perm);
                nnr = numel(P2nh2);  perm = randperm(nnr);  P2nh2 = P2nh2(perm);
                nnr = numel(P3nh2);  perm = randperm(nnr);  P3nh2 = P3nh2(perm);
                ind = P1nh2<=irev;  
                xP1 = P1xh2(ind);  yP1 = P1yh2(ind);  nP1 = P1nh2(ind);
                ind = P2nh2<=irev;  
                xP2 = P2xh2(ind);  yP2 = P2yh2(ind);  nP2 = P2nh2(ind);
                ind = P3nh2<=irev;  
                xP3 = P3xh2(ind);  yP3 = P3yh2(ind);  nP3 = P3nh2(ind);
                     
                P1map = get_rev_map_mod3(xP1, yP1, nP1, GaussFilter, nmap);
                P2map = get_rev_map_mod3(xP2, yP2, nP2, GaussFilter, nmap);
                P3map = get_rev_map_mod3(xP3, yP3, nP3, GaussFilter, nmap);
                Pavmap = (P1map+P2map+P3map)/3;
                Q1map = get_rev_map_mod3(xQ1, yQ1, nQ1, GaussFilter, nmap);                
                Q2map = get_rev_map_mod3(xQ2, yQ2, nQ2, GaussFilter, nmap);                
                Q3map = get_rev_map_mod3(xQ3, yQ3, nQ3, GaussFilter, nmap);
                Qavmap = (Q1map+Q2map+Q3map)/3;
                % mean-correct
                P1map=P1map-Pavmap;  P2map=P2map-Pavmap;  P3map=P3map-Pavmap;
                Q1map=Q1map-Qavmap;  Q2map=Q2map-Qavmap;  Q3map=Q3map-Qavmap;
                % P-Q 6 cases: PA-PA, SH-SH, SV-SV
                %              PA-SH, PA-SV, SH-PA, SH-SV, SV-PA, SV-SH
                if(zz==1)      P=P1map;  Q=Q1map;
                elseif(zz==2)  P=P2map;  Q=Q2map;
                elseif(zz==3)  P=P3map;  Q=Q3map;
                elseif(zz==4)  P=P1map;  Q=Q2map;
                elseif(zz==5)  P=P1map;  Q=Q3map;
                elseif(zz==6)  P=P2map;  Q=Q1map;
                elseif(zz==7)  P=P2map;  Q=Q3map;
                elseif(zz==8)  P=P3map;  Q=Q1map;
                elseif(zz==9)  P=P3map;  Q=Q2map;
                end
                % correlation
                c = corrcoef(P(:),Q(:)); %2-by-2
                % record SMboot
                SMboot(irev,iboot,zz,1,parti) = c(1,2);
                %fprintf('c=%d.\n',c(1,2));    
            end
        end
    fprintf('zz %d/%d.\n',zz, npatt);    
    end
    fprintf('parti=%d\n',parti);    
end

end

%%
function [revPosX, revPosY, revPosN] = get_rev_position(D)

ntr = D.Trials;

revx = D.RevealPosX;
revy = D.RevealPosY;
revn = repmat(1:25, ntr, 1);

revx(isnan(revx)) = 0;
revy(isnan(revy)) = 0;

revPosX = revx(:);
revPosY = revy(:);
revPosN = revn(:);

ind = revPosX==0; %indices of 0
revPosX(ind) =[];
revPosY(ind) =[];
revPosN(ind) =[];

end

%%
function [revxO, revyO, revnO] = equilizeSampNum(revx, revy, revn)

revxO=[];
revyO=[];
revnO=[];

for ind =1:numel(revn)
    if revn(ind)<=5
        f = 12;
    elseif revn(ind)>5 && revn(ind)<=10
        f = 15;
    elseif revn(ind)>10 && revn(ind)<=15
        f = 20;
    elseif revn(ind)>15 && revn(ind)<=20
        f = 30;
    elseif revn(ind)>20 && revn(ind)<=25
        f = 60;
    end
    tempx = repmat(revx(ind),f,1);
    tempy = repmat(revy(ind),f,1);
    tempn = repmat(revn(ind),f,1);
    revxO = [revxO; tempx];
    revyO = [revyO; tempy];
    revnO = [revnO; tempn];
end

end

%% modified version of a_Get_RevealMap_mod2
function filterMap = get_rev_map_mod3(revx, revy, revn, GaussFilter, nmap)

npix = 770;
ImageSize = 20;
nbin = npix/nmap;
revMap = zeros(nmap,nmap);

npos = numel(revx);  %bootstrap size

%compute weights
% wr = (6 - ceil(revn/5))/5;  %hard-wired weights

xDx = 1+(npix-1)*(revx/ImageSize+0.5);
xDy = 1-(npix-1)*(revy/ImageSize-0.5);
xDx = max(1,round(xDx/nbin));
xDy = max(1,round(xDy/nbin));
for i=1:npos    
    %revMap(xDy(i),xDx(i)) = revMap(xDy(i),xDx(i)) + 1*wr(i); %x-->col, y-->row
    revMap(xDy(i),xDx(i)) = revMap(xDy(i),xDx(i)) + 1; %x-->col, y-->row
end

% for i=1:npos
%     xDx = 1+(npix-1)*(revx(i)/ImageSize+0.5);
%     xDy = 1-(npix-1)*(revy(i)/ImageSize-0.5);
%     xDx = max(1,round(xDx/nbin));
%     xDy = max(1,round(xDy/nbin));
%     revMap(xDy,xDx) = revMap(xDy,xDx) + 1*wr(i); %x-->col, y-->row
% end

filterMap = conv2(revMap,GaussFilter,'same');
filterMap = filterMap/sum(filterMap(:)); %normalize

end

        
    
    