function checkCorr_sy(REV)
% 2015-09-30: checking correlations

same = zeros(25,3,4); %rev, type, parti
diff = zeros(25,3,4);

for irev = 5:25
    [w, a] = getCorr(REV, irev);
    same(irev,:,:) = w;
    diff(irev,:,:) = a;
    fprintf('rev%d done\n', irev);
end

clf;
for parti = 1:4
    subplot(1,4,parti);
    hold on;
    plot(1:25, same(:,1,parti,:), 'r-');
    plot(1:25, same(:,2,parti), 'r:');
    plot(1:25, same(:,3,parti), 'r:');
    plot(1:25, diff(:,1,parti), 'b-');
    plot(1:25, diff(:,2,parti), 'b:');
    plot(1:25, diff(:,3,parti), 'b:');
    title(sprintf('parti%d', parti));
    axis([0, 25, -1, 1]);
end

end

%%
function [within, across] = getCorr(REV, nrev)
% within, across 2 layers: mean ,lowCI, highCI; parti

nsamp = 10; %bootstrap

within = zeros(3,4); %type, parti
across = zeros(3,4);

c = zeros(3,3,4,nsamp); %type, type, parti, nsamp
parfor isamp = 1:nsamp
    [revx1, revy1, revx2, revy2] = getRevXY(REV,nrev);
    %fprintf('getRevXY samp%d done\n', isamp);
    c(:,:,:,isamp) = corrMap_anal(revx1, revy1, revx2, revy2);
    fprintf('samp%d done\n', isamp);
end

within_s = (c(1,1,:,:) + c(2,2,:,:) + c(3,3,:,:))/3;
across_s = (c(1,2,:,:) + c(1,3,:,:) +...
          c(2,1,:,:) + c(2,3,:,:) +...
          c(3,1,:,:) + c(3,2,:,:))/6;
      
for parti = 1:4
    temp = within_s(1,1,parti,:);
    x = sort(temp(:));
    within(1,parti) = mean(x); %mean
    within(2,parti) = x(round(nsamp*0.05)); %lower
    within(3,parti) = x(round(nsamp*0.95)); %upper
    
    temp = across_s(1,1,parti,:);
    x = sort(temp(:));
    across(1,parti) = mean(x); %mean
    across(2,parti) = x(round(nsamp*0.05)); %lower
    across(3,parti) = x(round(nsamp*0.95)); %upper
end

end

%%
function [revx1, revy1, revx2, revy2] = getRevXY(REV,nrev)

nrev = 5;
nsamp = floor(375/2); %for constructing image

revx1 = zeros(nsamp,4,3); %nsamp, parti, type
revy1 = zeros(nsamp,4,3);
revx2 = zeros(nsamp,4,3);
revy2 = zeros(nsamp,4,3);

for parti = 1:4
    for type = 1:3
        % remove trials that are too short
        revx = REV(type,parti,1).RevealPosX(:,1:nrev);
        revy = REV(type,parti,1).RevealPosY(:,1:nrev);
        yes0 = revx(:,nrev) == 0;
        revx(yes0,:) = [];
        revy(yes0,:) = [];
        % trial preserving permutation
        [ntr, ~] = size(revx);
        perm = randperm(ntr);
        permx = revx(perm,:);
        permy = revy(perm,:);
        % sample from each revealing number
        n = floor(nsamp/nrev);
        perm1 = perm(1:n);
        perm2 = perm(n+1:2*n);
        x1 = permx(perm1,:);
        y1 = permy(perm1,:);
        x2 = permx(perm2,:);
        y2 = permy(perm2,:);
        x1 = x1(:);
        x2 = x2(:);
        y1 = y1(:);
        y2 = y2(:);
        % removed used trials
        permx(1:2*n,:) = [];
        permy(1:2*n,:) = [];
        % sample remainder points - not trial preserving
        n = rem(nsamp,nrev);
        % linearize and permute again
        permx = permx(:);
        permy = permy(:);
        perm = randperm(numel(permx));
        permx = permx(perm);
        permy = permy(perm);
        % sample from pool
        xr1 = permx(1:n);
        yr1 = permy(1:n);
        xr2 = permx(n+1:2*n);
        yr2 = permy(n+1:2*n);
        % put in
        revx1(:,parti,type) = [x1;xr1];
        revx2(:,parti,type) = [x2;xr2];
        revy1(:,parti,type) = [y1;yr1];
        revy2(:,parti,type) = [y2;yr2];
    end
end

end

%%
function corr = corrMap_anal(revx1, revy1, revx2, revy2)

corr = zeros(3,3,4); %type, type, parti
sig = 1*20/27.8;
var = sig.^2;

for parti = 1:4
    x1 = revx1(:,parti, :);
    y1 = revy1(:,parti, :);
    x2 = revx2(:,parti, :);
    y2 = revy2(:,parti, :);
    % lienarize
    x1 = x1(:);
    y1 = y1(:);
    x2 = x2(:);
    y2 = y2(:);
    % difference weights
    ns = numel(x1);
    % 3 col for 3 types
    w(:,1) = [2./3./ns.*ones(ns,1); -1./3./ns.*ones(ns,1); -1./3./ns.*ones(ns,1)]; 
    w(:,2) = [-1./3./ns.*ones(ns,1); 2./3./ns.*ones(ns,1); -1./3./ns.*ones(ns,1)]; 
    w(:,3) = [-1./3./ns.*ones(ns,1); -1./3./ns.*ones(ns,1); 2./3./ns.*ones(ns,1)];
    % compute matrix
    for it = 1:3
    for jt = 1:3
        % between 1 and 2
        dx = bsxfun(@minus,x1,x2');
        dy = bsxfun(@minus,y1,y2');
        temp = exp(-0.5*( dx.^2./var + dy.^2./var ) );
        temp = repmat(temp,3,3);
        wm = w(:,it)*w(:,jt)';
        num = wm.*temp;
        % between 1 and 1
        dx = bsxfun(@minus,x1,x1');
        dy = bsxfun(@minus,y1,y1');
        temp = exp(-0.5*( dx.^2./var + dy.^2./var ) );
        temp = repmat(temp,3,3);
        wm = w(:,it)*w(:,it)';
        den1 = wm.*temp;
        % between 2 and 2
        dx = bsxfun(@minus,x2,x2');
        dy = bsxfun(@minus,y2,y2');
        temp = exp(-0.5*( dx.^2./var + dy.^2./var ) );
        temp = repmat(temp,3,3);
        wm = w(:,jt)*w(:,jt)';
        den2 = wm.*temp;
        corr(it,jt,parti) = sum(num(:)) / sqrt(sum(den1(:))*sum(den2(:)));
        %this is still wrong
    end
    end    
end

end

%%
% function corr = corrMap(revx1, revy1, revx2, revy2)
% 
% corr = zeros(3,3,3); %type, type, parti
% 
% nmap = 77;
% gFilter = GaussFilter(nmap);
% filterMap1 = zeros(nmap, nmap, 3);
% filterMap2 = zeros(nmap, nmap, 3);
% diffMap1 = zeros(nmap, nmap, 3);
% diffMap2 = zeros(nmap, nmap, 3);
% 
% for parti = 1:3
%     for type = 1:3
%         x1 = revx1(:,parti, type);
%         x2 = revx2(:,parti, type);
%         y1 = revy1(:,parti, type);
%         y2 = revy2(:,parti, type);
%         filterMap1(:,:,type) = rev_map(x1, y1, gFilter, nmap);
%         filterMap2(:,:,type) = rev_map(x2, y2, gFilter, nmap);
%     end
%     meanMap1 = mean(filterMap1,3);    
%     meanMap2 = mean(filterMap2,3);    
%     meanMap = (meanMap1 + meanMap2)/2;    
%     for type = 1:3
%         diffMap1(:,:,type) = filterMap1(:,:,type) - meanMap;
%         diffMap2(:,:,type) = filterMap2(:,:,type) - meanMap;
%     end
%     for itype = 1:3
%         for jtype = 1:3
%             imap = diffMap1(:,:,itype);
%             jmap = diffMap2(:,:,jtype);
%             c = corrcoef(imap(:),jmap(:)); %2-by-2
%             corr(itype, jtype, parti) = c(1,2);
%         end
%     end
% end
% 
% end

%%
% function gFilter = GaussFilter(mapsize)
% 
% % mapsize: number of rows of map to be filtered
% npix=770;  bin=npix/mapsize;
% var=(20/bin)^2; %variance is fixed at 20 for npix=770
% n=ceil(sqrt(var)*6);
% [x1,x2]=meshgrid(1:n,1:n);
% z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % in pixel coordinate
% gFilter = zeros(n*n,1);
% for i=1:n*n
%     gFilter(i) = mvnpdf(z(i,:),[1,1]*n/2,[1,0;0,1]*var);
% end
% gFilter = reshape(gFilter,n,n);
% 
% end

%%
% function filterMap = rev_map(revx, revy, GaussFilter, nmap)
% 
% npix = 770;
% ImageSize = 20;
% nbin = npix/nmap;
% revMap = zeros(nmap,nmap);
% npos = numel(revx);
% 
% %compute weights
% % wr = (6 - ceil(revn/5))/5;  %hard-wired weights
% 
% xDx = 1+(npix-1)*(revx/ImageSize+0.5);
% xDy = 1-(npix-1)*(revy/ImageSize-0.5);
% xDx = max(1,round(xDx/nbin));
% xDy = max(1,round(xDy/nbin));
% for i=1:npos    
%     %revMap(xDy(i),xDx(i)) = revMap(xDy(i),xDx(i)) + 1*wr(i); %x-->col, y-->row
%     revMap(xDy(i),xDx(i)) = revMap(xDy(i),xDx(i)) + 1; %x-->col, y-->row
% end
% 
% filterMap = conv2(revMap,GaussFilter,'same');
% filterMap = filterMap/sum(filterMap(:)); %normalize
% 
% end

%% mainly for testing
% function corr = corrMap_old(REV)
% 
% corr = zeros(3,3,3); %type, type, parti
% 
% nmap = 77;
% gFilter = GaussFilter(nmap);
% filterMap1 = zeros(nmap, nmap, 3);
% filterMap2 = zeros(nmap, nmap, 3);
% diffMap1 = zeros(nmap, nmap, 3);
% diffMap2 = zeros(nmap, nmap, 3);
% 
% for parti = 1:3
%     for type = 1:3
% 
% %         [ntr, ~] = size(REV(type,parti,1).RevealPosX);
% %         perm = randperm(ntr);        
% %         perm1 = perm(1:75);
% %         perm2 = perm(76:150);
% %         revx1 = REV(type,parti,1).RevealPosX(perm1,1);
% %         revy1 = REV(type,parti,1).RevealPosY(perm1,1);
% %         revx2 = REV(type,parti,1).RevealPosX(perm2,1);
% %         revy2 = REV(type,parti,1).RevealPosY(perm2,1);
%         
%         [revx1, revy1, revx2, revy2] = getRevXY(REV,1);        
%         filterMap1(:,:,type) = rev_map(revx1(:,parti,type), revy1(:,parti,type), gFilter, nmap);
%         filterMap2(:,:,type) = rev_map(revx2(:,parti,type), revy2(:,parti,type), gFilter, nmap);
%     end
%     meanMap1 = mean(filterMap1,3);    
%     meanMap2 = mean(filterMap2,3);    
%     meanMap = (meanMap1 + meanMap2)/2;    
%     for type = 1:3
%         diffMap1(:,:,type) = filterMap1(:,:,type) - meanMap;
%         diffMap2(:,:,type) = filterMap2(:,:,type) - meanMap;
% %         subplot(3,3, (type-1)*3 + parti);
% %         colormap('gray');
% %         %plot(revx,revy,'.');
% %         %axis([-10, 10, -10, 10]);
% %         imagesc(diffMap(:,:,type));
% %         axis square
% %         title(sprintf('parti%d type%d',parti, type));        
%     end
%     for itype = 1:3
%         for jtype = 1:3
%             imap = diffMap1(:,:,itype);
%             jmap = diffMap2(:,:,jtype);
%             c = corrcoef(imap(:),jmap(:)); %2-by-2
%             corr(itype, jtype, parti) = c(1,2);
%         end
%     end
% end
% 
% end