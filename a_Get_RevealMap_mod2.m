%% a variant of a_Get_RevealMap with a single bootstrap
function revMap = a_Get_RevealMap_mod2(D,rev,tr,GaussFilter, n, mode)
% rev can be a vector -- specify the revealing indices
% tr can be a vector -- specify the trial indices
% mode 1 = bootstrap
% mode 2 = extract from all
% mode 3 = extract from all with center removed

npix=770; ImageSize=20;
nbin=npix/n;
revMap=zeros(n,n);

% create selection indices
nrev=numel(rev);
ntr=numel(tr);
revxv = reshape(D.RevealPosX(tr,rev),ntr*nrev,1); %revealing x vector
revyv = reshape(D.RevealPosY(tr,rev),ntr*nrev,1); %revealing y vector
revxv(isnan(revxv))=0;  revyv(isnan(revyv))=0;

if (mode==3) % remove starting at center bias
    center_ind = (sqrt(revxv.^2 + revyv.^2 ) < 0.37*0.5);
    revxv(center_ind)=0;  revyv(center_ind)=0;
end

nboot_ind = find(revxv); %indices of usable trials
nboot = numel(nboot_ind);  %bootstrap size
%disp(sprintf('nboot=%d;',nboot));

%compute weights
revnum = reshape(repmat(rev,ntr,1),ntr*nrev,1); %revealing number vector
% wr = (6 - ceil(revnum/5))/5; %hard-wired weight on revealing: only good if number of revealings balanced across trials
wr = 1./(6 - ceil(revnum/5)); % 2015-09-23: or is it this?
% wr = ceil(revnum/5); % used to use this wrong weight
revImID = reshape(repmat(D.ImageID(tr,1),1,nrev),ntr*nrev,1); %reshape(tr-by-nrev)
wz = revImID; %hard-wired weight on pattern z: only good if PA is twice SH,SV; then weight = ImageID

for i=1:nboot    
    if (mode==1)  
        ind = nboot_ind(ceil(rand(1)*nboot));
    elseif (mode>=2)
        ind=nboot_ind(i);
    end
    %xDx = 1+(npix-1)*(D.RevealPosX(trial,ivspec)/ImageSize+0.5);  % cm to pixel: from inverting last two lines of a_GET_BALD
    %xDy = 1-(npix-1)*(D.RevealPosY(trial,ivspec)/ImageSize-0.5);
    xDx = 1+(npix-1)*(revxv(ind)/ImageSize+0.5);
    xDy = 1-(npix-1)*(revyv(ind)/ImageSize-0.5);
    xDx = max(1,round(xDx/nbin));
    xDy = max(1,round(xDy/nbin));
    xDx = min(770,round(xDx/nbin));
    xDy = min(770,round(xDy/nbin));
    w=wr(ind)*wz(ind);
    revMap(xDy,xDx) = revMap(xDy,xDx) + 1*w; %x-->col, y-->row
end

revMap = conv2(revMap,GaussFilter,'same');
%revMap=revMap+1;  %for bias map
revMap = revMap/sum(revMap(:)); %normalize
