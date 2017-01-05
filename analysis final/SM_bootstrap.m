function [SMboot] = SM_bootstrap(REV)
%% modified from DataPrep.m -- calculate similarity measure (LONG!!)

nrev = 25;
nboot = 100; %number of bootstrap samples
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
%nmap = 770; %P and Q map size
GaussFilter = a_GaussFilter(nmap);

for parti = 1:nparti
for pq = 2 %1:2 %1:ncond %cases: human-human, human-optimal, human-noisy_BAS    
    if(pq==1)
        % Q1=P1;  Q2=P2;  Q3=P3;  %Qav=Pav;
        %avoid correlation at 1st revealing
        P1=REV(1,parti,1);  P2=REV(2,parti,1);  P3=REV(3,parti,1);  %Pav=REV(4,parti,1);
        [P1x, P1y, P1n] = get_rev_position(P1);
        [P2x, P2y, P2n] = get_rev_position(P2);
        [P3x, P3y, P3n] = get_rev_position(P3);        
    elseif(pq==2)
        P1=REV(1,parti,1);  P2=REV(2,parti,1);  P3=REV(3,parti,1);  %Pav=REV(4,parti,1);
        Q1=REV(1,parti,2);  Q2=REV(2,parti,2);  Q3=REV(3,parti,2);  %Qav=REV(4,parti,2);  %noiseless BAS
        [P1x, P1y, P1n] = get_rev_position(P1);
        [P2x, P2y, P2n] = get_rev_position(P2);
        [P3x, P3y, P3n] = get_rev_position(P3);        
        [Q1x, Q1y, Q1n] = get_rev_position(Q1);
        [Q2x, Q2y, Q2n] = get_rev_position(Q2);
        [Q3x, Q3y, Q3n] = get_rev_position(Q3);        
    elseif(pq==3)
        P1=REV(1,parti,1);  P2=REV(2,parti,1);  P3=REV(3,parti,1);  %Pav=REV(4,parti,1);
        Q1=REV(1,parti,3);  Q2=REV(2,parti,3);  Q3=REV(3,parti,3);  %Qav=REV(4,parti,3);  %noisy BAS   
        [P1x, P1y, P1n] = get_rev_position(P1);
        [P2x, P2y, P2n] = get_rev_position(P2);
        [P3x, P3y, P3n] = get_rev_position(P3);        
        [Q1x, Q1y, Q1n] = get_rev_position(Q1);
        [Q2x, Q2y, Q2n] = get_rev_position(Q2);
        [Q3x, Q3y, Q3n] = get_rev_position(Q3);        
    end
    for zz = 1:npatt
        for irev = 1:nrev %1:nrev number of revealings
            ind = P1n<=irev;  xP1 = P1x(ind);  yP1 = P1y(ind);  nP1 = P1n(ind);
            ind = P2n<=irev;  xP2 = P2x(ind);  yP2 = P2y(ind);  nP2 = P2n(ind);
            ind = P3n<=irev;  xP3 = P3x(ind);  yP3 = P3y(ind);  nP3 = P3n(ind);
            if(pq>1)
                ind = Q1n<=irev;  xQ1 = Q1x(ind);  yQ1 = Q1y(ind);  nQ1 = Q1n(ind);
                ind = Q2n<=irev;  xQ2 = Q2x(ind);  yQ2 = Q2y(ind);  nQ2 = Q2n(ind);
                ind = Q3n<=irev;  xQ3 = Q3x(ind);  yQ3 = Q3y(ind);  nQ3 = Q3n(ind);                    
            end    
            for iboot = 1:nboot %number of bootstraps
                % bootstrapped sample                
                if(pq==1)
                    % 1: PA
                    nnr = numel(nP1);
                    perm = randperm(nnr);
                    nnrh = round(nnr/2);
                    sxP1 = xP1(perm(1:nnrh));
                    syP1 = yP1(perm(1:nnrh));
                    snP1 = nP1(perm(1:nnrh));
                    sxQ1 = xP1(perm(nnrh+1:nnr));
                    syQ1 = yP1(perm(nnrh+1:nnr));
                    snQ1 = nP1(perm(nnrh+1:nnr));
                    % 2: SH
                    nnr = numel(nP2);
                    perm = randperm(nnr);
                    nnrh = round(nnr/2);
                    sxP2 = xP2(perm(1:nnrh));
                    syP2 = yP2(perm(1:nnrh));
                    snP2 = nP2(perm(1:nnrh));
                    sxQ2 = xP2(perm(nnrh+1:nnr));
                    syQ2 = yP2(perm(nnrh+1:nnr));
                    snQ2 = nP2(perm(nnrh+1:nnr));
                    % 3: SV
                    nnr = numel(nP3);
                    perm = randperm(nnr);
                    nnrh = round(nnr/2);
                    sxP3 = xP3(perm(1:nnrh));
                    syP3 = yP3(perm(1:nnrh));
                    snP3 = nP3(perm(1:nnrh));
                    sxQ3 = xP3(perm(nnrh+1:nnr));
                    syQ3 = yP3(perm(nnrh+1:nnr));
                    snQ3 = nP3(perm(nnrh+1:nnr));
                elseif (pq>1)
                    % 1: PA
                    nnr = numel(nP1);
                    perm = randperm(nnr);  % note: classical bootstrap with replacement doesn't change thigns much
                    nnrh = round(nnr/2);
                    sxP1 = xP1(perm(1:nnrh));
                    syP1 = yP1(perm(1:nnrh));
                    snP1 = nP1(perm(1:nnrh));
                    nnr = numel(nQ1);
                    perm = randperm(nnr);
                    nnrh = round(nnr/2);
                    sxQ1 = xQ1(perm(1:nnrh));
                    syQ1 = yQ1(perm(1:nnrh));
                    snQ1 = nQ1(perm(1:nnrh));
                    % 2: SH
                    nnr = numel(nP2);
                    perm = randperm(nnr);
                    nnrh = round(nnr/2);
                    sxP2 = xP2(perm(1:nnrh));
                    syP2 = yP2(perm(1:nnrh));
                    snP2 = nP2(perm(1:nnrh));
                    nnr = numel(nQ2);
                    perm = randperm(nnr);
                    nnrh = round(nnr/2);
                    sxQ2 = xQ2(perm(1:nnrh));
                    syQ2 = yQ2(perm(1:nnrh));
                    snQ2 = nQ2(perm(1:nnrh));
                    % 3: SV
                    nnr = numel(nP3);
                    perm = randperm(nnr);
                    nnrh = round(nnr/2);
                    sxP3 = xP3(perm(1:nnrh));
                    syP3 = yP3(perm(1:nnrh));
                    snP3 = nP3(perm(1:nnrh));
                    nnr = numel(nQ3);
                    perm = randperm(nnr);
                    nnrh = round(nnr/2);
                    sxQ3 = xQ3(perm(1:nnrh));
                    syQ3 = yQ3(perm(1:nnrh));
                    snQ3 = nQ3(perm(1:nnrh));
                end
                P1map = get_rev_map_mod3(sxP1, syP1, snP1, GaussFilter, nmap);
                P2map = get_rev_map_mod3(sxP2, syP2, snP2, GaussFilter, nmap);
                P3map = get_rev_map_mod3(sxP3, syP3, snP3, GaussFilter, nmap);
                Pavmap = (P1map+P2map+P3map)/3;
                Q1map = get_rev_map_mod3(sxQ1, syQ1, snQ1, GaussFilter, nmap);                
                Q2map = get_rev_map_mod3(sxQ2, syQ2, snQ2, GaussFilter, nmap);                
                Q3map = get_rev_map_mod3(sxQ3, syQ3, snQ3, GaussFilter, nmap);
                Qavmap = (Q1map+Q2map+Q3map)/3;
                % mean-correct
                P1map=P1map-Pavmap;  P2map=P2map-Pavmap;  P3map=P3map-Pavmap;
                Q1map=Q1map-Qavmap;  Q2map=Q2map-Qavmap;  Q3map=Q3map-Qavmap;
                %testing
                %subplot(1,2,1);
                %imagesc(P2map);
                %subplot(1,2,2);
                %imagesc(Q2map);
                %pause
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
%                 if (0)
%                     subplot(1,2,1); imagesc(P); axis square
%                     subplot(1,2,2); imagesc(Q); axis square
%                     pause
%                 end
                % correlation
                c = corrcoef(P(:),Q(:)); %2-by-2
                % record SMboot
                SMboot(irev,iboot,zz,pq,parti) = c(1,2);
                %fprintf('c=%d.\n',c(1,2));    
            end
        end
    fprintf('zz %d/%d.\n',zz, npatt);    
    end
    fprintf('parti=%d; pq=%d.\n',parti,pq);    
end
end

%save('SMbootC.mat','SMboot');


end

%% get revealing positinos
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

%% modified version of a_Get_RevealMap_mod2
function filterMap = get_rev_map_mod3(revx, revy, revn, GaussFilter, nmap)

npix = 770;
ImageSize = 20;
nbin = npix/nmap;
revMap = zeros(nmap,nmap);

npos = numel(revx);  %bootstrap size

%compute weights
% wr = (6 - ceil(revn/5))/5;  %hard-wired weights
wr = 1./(6 - ceil(revn/5));  %2015-09-24

% for i=1:npos    
%     xDx = 1+(npix-1)*(revx(i)/ImageSize+0.5);
%     xDy = 1-(npix-1)*(revy(i)/ImageSize-0.5);
%     xDx = max(1,round(xDx/nbin));
%     xDy = max(1,round(xDy/nbin));
%     revMap(xDy,xDx) = revMap(xDy,xDx) + 1*wr(i); %x-->col, y-->row
% end

xDx = 1+(npix-1)*(revx/ImageSize+0.5);
xDy = 1-(npix-1)*(revy/ImageSize-0.5);
xDx = max(1,round(xDx/nbin));
xDy = max(1,round(xDy/nbin));
for i=1:npos    
    revMap(xDy(i),xDx(i)) = revMap(xDy(i),xDx(i)) + 1*wr(i); %x-->col, y-->row
end

filterMap = conv2(revMap,GaussFilter,'same');
filterMap = filterMap/sum(filterMap(:)); %normalize

end
