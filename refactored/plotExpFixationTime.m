function plotExpFixationTime(D, lineSpec)
    interSacTime = D.RevealTime;
    interSacTime(interSacTime == 0) = nan;
    for ind = 1:600
        interSacTime(ind,:) = interSacTime(ind,:) - interSacTime(ind,1);
    end
    mat = diff(interSacTime, 1, 2);
    mat(:,1) = nan;
    colN = size(mat, 2);
    avgVec = nanmean(mat, 1);
    varVec = nanvar(mat, 1);
    numVec = sum(1-isnan(mat), 1);
    semVec = sqrt(varVec./numVec);
    errorbar(1:colN, avgVec, semVec, lineSpec);
    hold on;
    xlabel('Revealing number');
    ylabel('Inter-revealing time (sec)');
end