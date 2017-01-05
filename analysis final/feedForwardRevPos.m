function [ffAL] = feedForwardRevPos(AL)

ffAL = AL;
nrev = max(AL.MaxRevealingTrial); %=25

for irev = 1:nrev
    pooled(irev).revInd = find(AL.RevealPosX(:,irev)~=0); %extract revealed ind
    pooled(irev).maxN = numel(pooled(irev).revInd);
end

ntrial = AL.Trials; %=600
for itrial = 1:ntrial
    nrev = AL.MaxRevealingTrial(itrial);
    for irev = 1:nrev
        randInd = ceil(rand()*pooled(irev).maxN);
        ind = pooled(irev).revInd(randInd); %use only revealed ind
        ffAL.RevealPosX(itrial, irev) = AL.RevealPosX(ind,irev);
        ffAL.RevealPosY(itrial, irev) = AL.RevealPosY(ind,irev);
    end
end

end
