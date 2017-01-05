%%
% toDo: clean structure duplication
function simuVis(simuObsMat, utilMeasMat)
    [noiseIndNum, costIndNum] = size(simuObsMat);
    for noiseInd = 1:noiseIndNum
        for costInd = 1:costIndNum
            simuObs = simuObsMat(noiseInd, costInd).simuObs;
            utilMeas = utilMeasMat(noiseInd, costInd).utilMeas;
            vis(noiseInd, costInd) = calVisQuant(simuObs, utilMeas);
        end
    end
    
    figure(1) % Time Vs Transition probability
    for noiseInd = 1:noiseIndNum
        subplot(2, noiseIndNum, noiseInd);        
            plotVisQuant(vis(noiseInd, 1).isStay, 'b-o');
            plotVisQuant(vis(noiseInd, 2).isStay, 'r-o');
        subplot(2, noiseIndNum, noiseIndNum+noiseInd);
            plotVisQuant(vis(noiseInd, 1).isReturn, 'b-o');
            plotVisQuant(vis(noiseInd, 2).isReturn, 'r-o');
    end
    set(gcf,'Color',[1,1,1]);
    % export_fig('TimeVsTransitionProb.pdf')
    
    figure(2) % Score gap property
    for noiseInd = 1:noiseIndNum
        subplot(2, noiseIndNum, noiseInd);
            plotVisQuant(vis(noiseInd, 1).logGapBas, 'b-o');
        subplot(2, noiseIndNum, noiseIndNum+noiseInd);
            histVisQuant(vis(noiseInd, 1).logGapBas);
    end
    set(gcf,'Color',[1,1,1]);
    % export_fig('ScoreGapProperty.pdf')

    figure(3) % Time vs Category certainty
    for noiseInd = 1:noiseIndNum
        subplot(1, noiseIndNum, noiseInd);
            plotVisQuant(vis(noiseInd, 1).certainty, 'b-o');
            plotVisQuant(vis(noiseInd, 2).certainty, 'r-o');
    end
    set(gcf,'Color',[1,1,1]);
    % export_fig('TimeVsCertainty.pdf')
    
    figure(4) % Fixation duration property
    for noiseInd = 1:noiseIndNum
        subplot(2, noiseIndNum, noiseInd);
%             plotPercentiles(vis(noiseInd, 1).fixateIndDuration.mat, 95, 'k-s');
            plotVisQuant(vis(noiseInd, 1).fixateIndDuration, 'b-o');
%             plotPercentiles(vis(noiseInd, 2).fixateIndDuration.mat, 95, 'k-s');
            plotVisQuant(vis(noiseInd, 2).fixateIndDuration, 'r-o');        
        subplot(2, noiseIndNum, noiseIndNum+noiseInd);
            plotVisQuant(vis(noiseInd, 1).fixateDuration, 'b-o');
            plotVisQuant(vis(noiseInd, 2).fixateDuration, 'r-o');        
    end
    set(gcf,'Color',[1,1,1]);
    %export_fig('FixationDurationProperty.pdf')

%     figure(5) % Fixation Trajectory
%     for noiseInd = 1:noiseIndNum
%         subplot(1, noiseIndNum, noiseInd);
%             plotPosIndVsTime(vis(noiseInd, 1).uniquePosInd.mat, [0,0,1]);
%             plotPosIndVsTime(vis(noiseInd, 2).uniquePosInd.mat, [1,0,0]);
%     end
%     set(gcf,'Color',[1,1,1]);
%     %export_fig('FixationTrajectory.pdf')
end

% calVisQuant group start
function vis = calVisQuant(simuObs, utilMeas)
    trialNum = numel(simuObs);
    maxObsTime = findMaxField(simuObs, 'time');
    mat = nan(trialNum, maxObsTime);
    vis.isStay.mat = mat;
    vis.isReturn.mat = mat;
    vis.logMaxBas.mat = mat;
    vis.logGapBas.mat = mat;
    vis.certainty.mat = mat;
    vis.uniquePosInd.mat = mat;
    vis.fixateIndDuration.mat = mat;
    vis.fixateDuration.mat = mat;
    for trialInd = 1:trialNum
        trialObs = simuObs(trialInd);
        endInd = trialObs.time;
        vis.isStay.mat(trialInd, 1:endInd) = calStayProb(trialObs);
        vis.isReturn.mat(trialInd, 1:endInd) = calReturnProb(trialObs);
        vis.logMaxBas.mat(trialInd, 1:endInd) = log(utilMeas(trialInd).maxBas(1:endInd));
        vis.logGapBas.mat(trialInd, 1:endInd) = log(utilMeas(trialInd).gapBas(1:endInd));
        vis.certainty.mat(trialInd, 1:endInd) = utilMeas(trialInd).certainty(1:endInd);
        vis.uniquePosInd.mat(trialInd, 1:endInd) = calUniquePos(trialObs);
        fixateInd = calFixateInd(trialObs);
        [fixateDuration, fixateNum] = calFixationDuration(fixateInd);
        vis.fixateIndDuration.mat(trialInd, 1:fixateNum) = fixateDuration;
        vis.fixateDuration.mat(trialInd, 1:endInd) = fixateDuration(fixateInd);
    end
    vis.isStay.yName = 'Stay probability';
    vis.isReturn.yName = 'Return probability';
    vis.logMaxBas.yName = 'log(Max BAS score)';
    vis.logGapBas.yName = 'log(Score gap)';
    vis.certainty.yName = 'Category certainty';
    vis.uniquePosInd.yName = 'Position index';
    vis.fixateIndDuration.yName = 'Average duration';
    vis.fixateDuration.yName = 'Average duration';
    
    vis.isStay.yLim = [0, 1];
    vis.isReturn.yLim = [0, 1];
    vis.logMaxBas.yLim = [-10, 0];
    vis.logGapBas.yLim = [-10, 0];
    vis.certainty.yLim = [0.5, 1];
    vis.uniquePosInd.yLim = calYLimits(vis.uniquePosInd.mat);
    vis.fixateIndDuration.yLim = calYLimits(vis.fixateIndDuration.mat);
    vis.fixateDuration.yLim = calYLimits(vis.fixateDuration.mat);
        
    vis.isStay.xName = 'Time';
    vis.isReturn.xName = 'Time';
    vis.logMaxBas.xName = 'Time';
    vis.logGapBas.xName = 'Time';
    vis.certainty.xName = 'Time';
    vis.uniquePosInd.xName = 'Time';
    vis.fixateIndDuration.xName = 'Fixation Index';
    vis.fixateDuration.xName = 'Time';      

    vis.isStay.xLim = calXLimits(vis.isStay.mat);
    vis.isReturn.xLim = calXLimits(vis.isReturn.mat);
    vis.logMaxBas.xLim = calXLimits(vis.logMaxBas.mat);
    vis.logGapBas.xLim = calXLimits(vis.logGapBas.mat);
    vis.certainty.xLim = calXLimits(vis.certainty.mat);
    vis.uniquePosInd.xLim = calXLimits(vis.uniquePosInd.mat);
    vis.fixateIndDuration.xLim = calXLimits(vis.fixateIndDuration.mat);
    vis.fixateDuration.xLim = calXLimits(vis.fixateDuration.mat);    
end

function maxTime = findMaxField(simuObs, field)
    trialNum = numel(simuObs);
    fieldVec = zeros(trialNum, 1);
    for trialInd = 1:trialNum
       fieldVec(trialInd) = simuObs(trialInd).(field); 
    end
    maxTime = max(fieldVec);
end

function isStay = calStayProb(trialObs)
    rawGridInd = trialObs.rawGridInd;
    maxObsTime = trialObs.time;
    isStay = zeros(maxObsTime, 1);
    for timeInd = 2:maxObsTime
       if ( rawGridInd(timeInd) == rawGridInd(timeInd-1))
           isStay(timeInd) = 1;
       else
           isStay(timeInd) = 0;
       end
    end
end

function isReturn = calReturnProb(trialObs)
    rawGridInd = trialObs.rawGridInd;
    maxObsTime = trialObs.time;
    isReturn = zeros(maxObsTime, 1);
    for timeInd = 2:maxObsTime
       oldIndNum = sum(rawGridInd(1:timeInd-1) == rawGridInd(timeInd));
       isSameInd = rawGridInd(timeInd) == rawGridInd(timeInd-1);
       if ( oldIndNum > 0 && isSameInd == 0 )
           isReturn(timeInd) = 1;
       else
           isReturn(timeInd) = 0;
       end
    end
end

function uniquePosInd = calUniquePos(trialObs)
    field = trialObs.rawPosXY;
    [~, ~, uniquePosInd] = unique(field, 'rows', 'stable');
end

function fixateInd = calFixateInd(trialObs)
    rawGridInd = trialObs.rawGridInd;
    maxObsTime = trialObs.time;
    fixateInd = ones(maxObsTime, 1);
    for timeInd = 1:maxObsTime-1
       if ( rawGridInd(timeInd) == rawGridInd(timeInd+1))
           fixateInd(timeInd+1) = fixateInd(timeInd);
       else
           fixateInd(timeInd+1) = fixateInd(timeInd) + 1;
       end
    end
end

function [fixateDuration, fixateNum] = calFixationDuration(fixateInd)
    fixateNum = max(fixateInd);
    fixateDuration = nan(fixateNum, 1);
    for ind = 1:fixateNum
        fixateDuration(ind) = sum(fixateInd == ind);
    end
    fixateDuration(fixateDuration>50) = nan; %bad: magic number 50 to kill outliers
end

function yLim = calYLimits(quantMat)
    yLim = [nanmin(quantMat(:)), nanmax(quantMat(:))/2]; % bad: magic number 2
    if ( yLim(2) == yLim(1) )
        yLim(2) = yLim(1) + 1;
    end
    if ( sum(1-isreal(quantMat(:))) )
        fprintf('Bug: imaginary number exists!\n')
    end
end

function xLim = calXLimits(quantMat)
    nMat = 1 - isnan(quantMat);
    nVec = sum(nMat, 1);
    thresh = 1;
    nOverThresh = (nVec > thresh);
    xLim = [1, find(nOverThresh, 1, 'last')];
end
% group end

function plotVisQuant(quant, lineSpec)
    plotAvgSem(quant.mat, lineSpec);
    xlabel(quant.xName);
    ylabel(quant.yName);
%     quant.xName
%     quant.yName
%     quant.xLim
%     quant.yLim
    xlim(quant.xLim);
    ylim(quant.yLim);
end

function plotAvgSem(mat, lineSpec)
    mat(isinf(mat)) = nan;
    %mat(~isreal(mat)) = nan;
    colN = size(mat, 2);
    avgVec = nanmean(mat, 1);
    varVec = nanvar(mat, 1);
    numVec = sum(1-isnan(mat), 1);
    semVec = sqrt(varVec./numVec);
    plot(1:colN, avgVec, lineSpec);
%     errorbar(1:colN, avgVec, semVec, lineSpec);
    hold on;
end

function plotPercentiles(mat, percentile, lineSpec)
    [~, colN] = size(mat);
    linePercentile = zeros(1, colN);
    for ind = 1:colN
        col = mat(:,ind);
        realN = max(1, sum(~isnan(col)));
        col = sort(col);
        linePercentile(ind) = col(ceil(realN*percentile/100));
    end
    plot(1:colN, linePercentile, lineSpec);
    hold on;
end

function histVisQuant(quant)
    data = quant.mat(:);
    data(isinf(data)) = nan;
    hist(data, 50, 'b');
    hold on;
    axis fill;
    xlabel(quant.yName);
    ylabel('Counts');
end

% plotPosIndVSTime group start
function plotPosIndVsTime(posIndMat, colorSpec)
    [trialNum, ~] = size(posIndMat);
    for trialInd = 1:trialNum
        uniqueInd = posIndMat(trialInd,:);
        time = 1:sum(~isnan(uniqueInd));
        plotTimeDuration(uniqueInd, time, trialInd/trialNum/2, colorSpec);
    end
    xlabel('Time');
    ylabel('Position Index');
end

function plotTimeDuration(aField, time, fieldShift, colorSpec)
    entryNum = numel(time);
    for ind = 1:entryNum
        if (ind < entryNum && aField(ind) == aField(ind+1))
            plot(time(ind:ind+1), aField(ind:ind+1) + fieldShift, '-', 'LineWidth', 1, 'Color', colorSpec);
        else
            plot(time(ind), aField(ind) + fieldShift, '.', 'MarkerSize', 5, 'Color', colorSpec);
        end
        hold on;
    end
end
% group end

