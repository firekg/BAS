%% a variant of a_Get_RevealMap with a single bootstrap
function revMap = a_Get_RevealMap_mod(D,rev,tr,GaussFilter, mode)
% rev can be a vector -- specify the revealing indices
% tr can be a vector -- specify the trial indices
% mode 1 = bootstrap; mode 2 = extract from all

npix=770; ImageSize=20; nbin=10;
n=npix/nbin;
revMap=zeros(n,n);

% create selection indices
nrev=numel(rev);
ntr=numel(tr);
revxv = reshape(D.RevealPosX(tr,rev),ntr*nrev,1); %revealing x vector
revyv = reshape(D.RevealPosY(tr,rev),ntr*nrev,1); %revealing y vector
nboot_ind = find(revxv); %indices of usable trials
nboot = numel(nboot_ind);  %bootstrap size

%disp(sprintf('nboot=%d;',nboot)); 

for i=1:nboot    
    if (mode==1)  
        ind = nboot_ind(ceil(rand(1)*nboot));
    elseif (mode==2) 
        ind=nboot_ind(i);
    end
    %xDx = 1+(npix-1)*(D.RevealPosX(trial,ivspec)/ImageSize+0.5);  % cm to pixel: from inverting last two lines of a_GET_BALD
    %xDy = 1-(npix-1)*(D.RevealPosY(trial,ivspec)/ImageSize-0.5);
    xDx = 1+(npix-1)*(revxv(ind)/ImageSize+0.5);
    xDy = 1-(npix-1)*(revyv(ind)/ImageSize-0.5);
    xDx = max(1,round(xDx/nbin));  % can somehow be 0
    xDy = max(1,round(xDy/nbin));
    revMap(xDy,xDx) = revMap(xDy,xDx) + 1; %x-->col, y-->row
end

revMap = conv2(revMap,GaussFilter,'same');
revMap = revMap/sum(sum(revMap)); %normalize
