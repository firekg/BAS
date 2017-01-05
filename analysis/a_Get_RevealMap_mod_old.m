%% a variant of a_Get_RevealMap with a single bootstrap
function revMap = a_Get_RevealMap_mod2(D,rev,tr,GaussFilter, n, mode)
% rev can be a vector -- specify the revealing indices
% tr can be a vector -- specify the trial indices
% mode 1 = bootstrap; mode 2 = extract from all

npix=770; ImageSize=20;
nbin=npix/n;
revMap=zeros(n,n);

% create selection indices
nrev=numel(rev);
ntr=numel(tr);
revxv = reshape(D.RevealPosX(tr,rev),ntr*nrev,1); %revealing x vector
revyv = reshape(D.RevealPosY(tr,rev),ntr*nrev,1); %revealing y vector
revxv(isnan(revxv))=0;  revyv(isnan(revyv))=0;
nboot_ind = find(revxv); %indices of usable trials
nboot = numel(nboot_ind);  %bootstrap size
%disp(sprintf('nboot=%d;',nboot));

%compute weights
revnum = reshape(repmat(rev,ntr,1),ntr*nrev,1); %revealing number vector
wr = ceil(revnum/5); %hard-wired weight on revealing: only good if number of revealings balanced across trials
revImID = reshape(repmat(D.ImageID(tr,1),1,nrev),ntr*nrev,1); %reshape(tr-by-nrev)
wz = revImID; %hard-wired weight on pattern z: only good if PA is twice SH,SV; then weight = ImageID

if (mode==1)  
    ind = nboot_ind(ceil(rand(nboot,1)*nboot));
elseif (mode==2)
    ind=nboot_ind;
end

for i=1:numel(ind)    
    %xDx = 1+(npix-1)*(D.RevealPosX(trial,ivspec)/ImageSize+0.5);  % cm to pixel: from inverting last two lines of a_GET_BALD
    %xDy = 1-(npix-1)*(D.RevealPosY(trial,ivspec)/ImageSize-0.5);
    xDx = 1+(npix-1)*(revxv(ind)/ImageSize+0.5);
    xDy = 1-(npix-1)*(revyv(ind)/ImageSize-0.5);
    xDx = max(1,round(xDx/nbin));
    xDy = max(1,round(xDy/nbin));
    w=wr(ind)*wz(ind);
    revMap(xDy,xDx) = revMap(xDy,xDx) + 1*w; %x-->col, y-->row
end

revMap = conv2(revMap,GaussFilter,'same');
%revMap=revMap+1;  %for bias map
revMap = revMap/sum(revMap(:)); %normalize
