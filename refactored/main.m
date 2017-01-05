%%
% refactoring conventions
% function name stricture: verb-adj-noun-prep
% compute is a collection of cals

function main(imageData)    
    param.disp.dispNoiseSD = 1;
    param.percept.obsNoiseSD = 0.1;
    param.percept.forgettingRate = nan;
    param.percept.accuityFallOffRate =  nan;
    param.percept.accuityRadius = nan;
    param.decis.softmaxRate = 0;
    param.decis.lapseRate = 0;
    param.decis.priorTypeVector = [0.5, 0.25, 0.25];
    param.move.cost = 0;
    param.simu.trialNum = 20;
    param.simu.obsTime = 10;
    
%     simulateBAS(imageData, param);

    dispNoiseSD = [0.5, 1, 2];
    moveCost = [0, exp(-7)];
    trialNum = [2, 2, 2];
    obsTime = [30, 30, 30];
    noiseLevelNum = numel(dispNoiseSD);
    costLevelNum = numel(moveCost);
    simuObsWithVaryingNoise = struct([]);
    utilMeasWithVaryingNoise = struct([]);
    for parInd = 1:noiseLevelNum
        for costLevelInd = 1:costLevelNum
            param.disp.dispNoiseSD = dispNoiseSD(parInd);
            param.move.cost = moveCost(costLevelInd);
            param.simu.trialNum = trialNum(parInd);
            param.simu.obsTime = obsTime(parInd);
            [simuObs, utilMeas] = simulateBAS(imageData, param);
            simuObsWithVaryingNoise(parInd, costLevelInd).simuObs = simuObs;
            utilMeasWithVaryingNoise(parInd, costLevelInd).utilMeas = utilMeas;
        end
    end
    
    figure(1)
    for parInd = 1:noiseLevelNum
        % Position Index vs Time
        subplot(3, noiseLevelNum, parInd);
        plotPosIndVsTime(simuObsWithVaryingNoise(parInd, 1).simuObs, [0,0,1]);
        plotPosIndVsTime(simuObsWithVaryingNoise(parInd, 2).simuObs, [1,0,0]);
        axis([0 obsTime(parInd)+1 0 obsTime(parInd)+1]);
        axis fill;
        xlabel('Position Index');
        ylabel('Time');
        
        % Fixation Index vs Average duration
        subplot(3, noiseLevelNum, noiseLevelNum+parInd);
        plotDuration(simuObsWithVaryingNoise(parInd, 1).simuObs, 'b-o');
        plotDuration(simuObsWithVaryingNoise(parInd, 2).simuObs, 'r-o');
        axis([0 obsTime(parInd)+1 1 10]);
        axis fill;
        set(gca, 'XTick', 0:5:obsTime(parInd));
        set(gca, 'XTickLabel', 0:5:obsTime(parInd));
        xlabel('Fixation Index');
        ylabel('Average duration');
        
        % Time probe vs Average duration
        subplot(3, noiseLevelNum, 2*noiseLevelNum+parInd);
        probeInterval = 1;
        timeProbe = 1:probeInterval:obsTime(parInd);
        plotDuration(simuObsWithVaryingNoise(parInd, 1).simuObs, 'b-o', timeProbe);
        plotDuration(simuObsWithVaryingNoise(parInd, 2).simuObs, 'r-o', timeProbe);
        probeNum = numel(timeProbe); 
        axis([0 probeNum+3 1 10]); % bad: mageic number 3
        axis fill;
        set(gca, 'XTick', 0:5/probeInterval:probeNum+1);
        set(gca, 'XTickLabel', (0:5/probeInterval:probeNum+1)*probeInterval);
        xlabel('Time probe');
        ylabel('Average duration');        
    end
    set(gcf,'Color',[1,1,1]);
    % export_fig('FixationProperty.pdf')
    
    figure(2)
    for parInd = 1:noiseLevelNum
        subplot(3, noiseLevelNum, parInd);
        logIt = 1;
        utilMeas = utilMeasWithVaryingNoise(parInd, 1).utilMeas;
        plotTimeVsUtilField(utilMeas, 'maxBas', logIt, 'b-o');
        utilMeas = utilMeasWithVaryingNoise(parInd, 2).utilMeas;
        plotTimeVsUtilField(utilMeas, 'maxBas', logIt, 'r-o');
        axis([0 obsTime(parInd)+1 -15 0]);
        axis fill;
        set(gca, 'XTick', 0:5:obsTime(parInd));
        set(gca, 'XTickLabel', 0:5:obsTime(parInd));
        xlabel('Time');
        ylabel('log(Maximum BAS score)');

        subplot(3, noiseLevelNum, noiseLevelNum+parInd);
        logIt = 0;
        utilMeas = utilMeasWithVaryingNoise(parInd, 1).utilMeas;
        plotTimeVsUtilField(utilMeas, 'gapBas', logIt, 'b-o');
        utilMeas = utilMeasWithVaryingNoise(parInd, 2).utilMeas;
        plotTimeVsUtilField(utilMeas, 'gapBas', logIt, 'r-o');
        %axis([0 obsTime(parInd)+1 -15 0]);
        axis fill;
        set(gca, 'XTick', 0:5:obsTime(parInd));
        set(gca, 'XTickLabel', 0:5:obsTime(parInd));
        xlabel('Time');
        ylabel('log(Score gap)');
        
        subplot(3, noiseLevelNum, 2*noiseLevelNum+parInd);
        logIt = 0;
        utilMeas = utilMeasWithVaryingNoise(parInd, 1).utilMeas;
        plotTimeVsUtilField(utilMeas, 'certainty', logIt, 'b-o');
        utilMeas = utilMeasWithVaryingNoise(parInd, 2).utilMeas;
        plotTimeVsUtilField(utilMeas, 'certainty', logIt, 'r-o');
        axis([0 obsTime(parInd)+1 0.5 1]);
        axis fill;
        set(gca, 'XTick', 0:5:obsTime(parInd));
        set(gca, 'XTickLabel', 0:5:obsTime(parInd));
        xlabel('Time');
        ylabel('Category certainty');                
    end
    set(gcf,'Color',[1,1,1]);
    % export_fig('TimeVsUtilMeas.pdf')
end

function [simuObs, utilMeasure] = simulateBAS(imageData, param)
% ToDo: combine recordSimuBAS & watchSimuBAS
% ToDo: test results against old version
% ToAdd: introduce time into obs?
% ToFix: dispNoiseSD = 100, non-smooth BAS scores
% ToCheck: tiny dispNoiseSD --> no returns
% ToAdd: error detection of map calculation
% ToAdd: [paramGrad] = calLogLikelihoodGradients()
    if (nargout==0)
        mode = 'watch';
    elseif (nargout ==2)
        mode = 'record';
    end
    probeGridXY = prepProbeGridXY();
    trialNum = param.simu.trialNum;
    obsTime = param.simu.obsTime;
    if (strcmp(mode, 'watch'))
        figureHandle = setupSimuPlot();
    elseif(strcmp(mode, 'record'))
        % toDo: combine into initUtilMeas
        maxBas = zeros(obsTime, 1);
        gapBas = zeros(obsTime, 1);
        certainty = zeros(obsTime, 1);
    end
    for trialInd = 1:trialNum
        [newObs, actualImage] = initSimuTrial(imageData, probeGridXY, trialInd);
        obs = initObs();
        for obsInd = 1:obsTime
            obs = updateObs(obs, newObs, actualImage);
            gaussMix = computeExpectedPixDistr(obs, param, probeGridXY);
            [newObs, scoreMapLinearized] = findPosWithMaxScore(gaussMix, probeGridXY, obs, param);
            if (strcmp(mode, 'watch'))
                plotSimuStep([obs.posXY; newObs.posXY], actualImage, scoreMapLinearized);
                if strcmp(get(figureHandle,'currentcharacter'), 'q')
                    break
                end
            elseif (strcmp(mode, 'record'))
                % toDo: combine into calUtilMeas
                maxBas(obsInd) = max(scoreMapLinearized); 
                gapBas(obsInd) = max(scoreMapLinearized)-scoreMapLinearized(obs.rawGridInd(obsInd));
                [~, logPosteriorCategP] = calChoiceProb(gaussMix.logWeights, param.decis);
                certainty(obsInd) = max( 1-exp(logPosteriorCategP), exp(logPosteriorCategP) );
            end
        end
        if (strcmp(mode, 'watch'))
            if strcmp(get(figureHandle,'currentcharacter'), 'q')
                break
            end
        elseif (strcmp(mode, 'record'))            
            simuObs(trialInd) = obs;
            % toDo: combine into recordUtilMeas
            utilMeasure(trialInd).maxBas = maxBas;
            utilMeasure(trialInd).gapBas = gapBas;
            utilMeasure(trialInd).certainty = certainty;
        end
    end
    if (strcmp(mode, 'watch'))
        close(figureHandle);
    end
end

function probeGridXY = prepProbeGridXY()
    expConst = defineExperimentConstants();
    probeRowNum = 110;
    pointsAlongDim = linspace(-expConst.imageLengthCm/2, expConst.imageLengthCm/2, probeRowNum);
    [gridX, gridY] = meshgrid(pointsAlongDim,pointsAlongDim);
    probeGridXY = [reshape(gridX, probeRowNum^2, 1), reshape(gridY, probeRowNum^2, 1)];
end

function expConst = defineExperimentConstants()
    expConst.pixelsPerDim = 770;
    expConst.pixelNum = 770^2;
    expConst.imageLengthCm = 20;  % screen size = [-10,10] cm
end

% initSimuTrial group start
function [newObs, actualImage] = initSimuTrial(imageData, gridXY, trialInd)
    pointXY = [0, 0];
    newObs.gridInd = findClosestInd(gridXY, pointXY);
    newObs.posXY = gridXY(newObs.gridInd, :);
    actualImage = imageData(trialInd).yimage;
end

function ind = findClosestInd(gridXY, pointXY)
    [~, ind] = min( (gridXY(:,1)-pointXY(1,1)).^2 + (gridXY(:,2)-pointXY(1,2)).^2 );
end
% group end

function obs = initObs()
    obs.posXY = []; % these 6 are for simulation efficiency
    obs.gridInd = [];
    obs.stayTime = [];
    obs.actualPix = [];
    obs.dispRandn = [];
    obs.num = 0;
    obs.rawPosXY = []; % these 3 are for plotting
    obs.rawGridInd = [];
    obs.time = 0;
end

% appendToObs group start
function obsOut = updateObs(obsIn, newObs, actualImage)
    if (obsIn.num == 0)
        ind = 0;
    else
        ind = obsIn.gridInd == newObs.gridInd;
    end
    if (sum(ind)==0)
        mode = 'newObs';
    elseif (sum(ind)==1)
        mode = 'returnObs';
    end
    if (strcmp(mode, 'newObs'))
        obsOut.posXY = [obsIn.posXY; newObs.posXY];
        obsOut.gridInd = [obsIn.gridInd; newObs.gridInd];
        obsOut.stayTime = [obsIn.stayTime; 1];
        newActualPix = getActualPix(newObs.posXY(1,1), newObs.posXY(1,2), actualImage);
        obsOut.actualPix = [obsIn.actualPix; newActualPix];
        obsOut.dispRandn = [obsIn.dispRandn; randn()];
        obsOut.num = obsIn.num + 1;
    elseif (strcmp(mode, 'returnObs'))
        obsOut = obsIn;
        obsOut.stayTime(ind) = obsIn.stayTime(ind) + 1;    
        obsOut.dispRandn(ind) = randn(); % overwritten value not recorded!
    else
        fprintf('unexpected behavior in addToObsPosTime()');
    end    
    obsOut.rawPosXY = [obsIn.rawPosXY; newObs.posXY];
    obsOut.rawGridInd = [obsIn.rawGridInd; newObs.gridInd];
    obsOut.time = obsIn.time + 1;
end

function actualPix = getActualPix(posX, posY, actualImage)
    expConst = defineExperimentConstants();
    [posRow, posCol] = transformPosCm2Ind(posX, posY);
    linearInd = sub2ind([expConst.pixelsPerDim, expConst.pixelsPerDim], posRow, posCol);
    actualPix = actualImage(linearInd);
end

function [posRow, posCol] = transformPosCm2Ind(posX, posY)
    expConst = defineExperimentConstants(); % don't really need this    
    pixelsPerDim = expConst.pixelsPerDim;
    imageLengthCm = expConst.imageLengthCm;
    posCol = round( 1+(pixelsPerDim-1)*(posX/imageLengthCm+0.5) );
    posRow = round( 1-(pixelsPerDim-1)*(posY/imageLengthCm-0.5) );
end
% group end

function gaussMix = computeExpectedPixDistr(obs, param, probeGridXY)
    dispParam = param.disp;
    decisParam = param.decis;
    displayedPixVariance = calVarianceOnPoints(obs, dispParam);
    expectedProbesVariance = calVarianceOnPoints(obs, dispParam, probeGridXY);
    [logTypeLikelihoods, gaussMix] = calGPPredictions(obs, displayedPixVariance, probeGridXY, expectedProbesVariance);
    gaussMix.logWeights = calLogTypePosteriors(logTypeLikelihoods, decisParam);
end

% calVarianceOnPoints group start
function pointsVariance = calVarianceOnPoints(obs, dispParam, probeGrid)
% TODO: include addedImageVariance = 0.01^2 somewhere?
% maybe: modulateTemporalPrecision()
% maybe: modulateSpatialPrecision()
% maybe: combineSpatiotemporalPrecision()
% maybe: transformPrecision2Variance()
    if (nargin == 2)
        mode = 'obs';
    elseif (nargin == 3)
        mode = 'probe';
    end
    if (strcmp(mode, 'obs'))
        displayedPixVariance = calVarianceOnDisplayedPix(obs, dispParam);
        pointsVariance = displayedPixVariance;
    elseif (strcmp(mode, 'probe'))
        expectedProbesVariance = calVarianceOnProbes(obs, probeGrid, dispParam);
        pointsVariance = expectedProbesVariance;
    end
end

function displayedPixVariance = calVarianceOnDisplayedPix(obs, dispParam)
    displayedPixVariance = 1./obs.stayTime.*dispParam.dispNoiseSD^2;
    % fprintf('min displayed pixel variance = %f\n' ,min(displayedPixVariance));
    % fprintf('max displayed pixel variance = %f\n' ,max(displayedPixVariance));
end

function expectedProbesVariance = calVarianceOnProbes(obs, probeGrid, dispParam)
% maybe: natural perception noise with smooth dip variance across probes
    probeNum = size(probeGrid, 1);
    expectedStayTime = ones(probeNum, 1);
    expectedStayTime(obs.gridInd) = expectedStayTime(obs.gridInd) + obs.stayTime;
    expectedProbesVariance = 1./expectedStayTime.*dispParam.dispNoiseSD^2;
    % fprintf('min expected probe variance = %f\n' ,min(expectedProbesVariance));
    % fprintf('max expected probe variance = %f\n' ,max(expectedProbesVariance));
end
% group end

% calGPPredictions group start
function varargout = calGPPredictions(obs, displayedPixVariance, probeGridXY, expectedProbesVariance) % messy
%  works in cm
    if (nargin == 2)
        mode = 'likelihood';        
    elseif (nargin == 4)
        mode = 'distribution';
    end
    imageTypeDef = defineImageType();
    if (strcmp(mode, 'likelihood'))
        logTypeLikelihoods = zeros(1, imageTypeDef.typeNum);        
    elseif (strcmp(mode,'distribution'))
       probeNum = size(probeGridXY, 1);
       gaussMix.means = zeros(probeNum, imageTypeDef.typeNum);
       gaussMix.variances = zeros(probeNum, imageTypeDef.typeNum);
    end
    obsPosXY = obs.posXY;
    actualPix = obs.actualPix;
    dispRandn = obs.dispRandn;
    obsNum = obs.num;
    displayedPix = calDisplayedPix(actualPix, dispRandn, displayedPixVariance);
    for imageTypeInd = 1:imageTypeDef.typeNum
        imageTypeHyp = imageTypeDef.hyp(imageTypeInd).def;
        actualPixCov = feval(imageTypeDef.covfunc{:}, imageTypeHyp, obsPosXY);
        L = chol(actualPixCov + diag(displayedPixVariance));  % named by C. Rausmussen
        alpha = solve_chol(L, displayedPix);  % named by C. Rausmussen  
        logTypeLikelihoods(imageTypeInd) = -((displayedPix)'*alpha/2 + sum(log(diag(L))) + obsNum*log(2*pi)/2);
        if (strcmp(mode, 'distribution'))
           btweenPixAndProbesCov = feval(imageTypeDef.covfunc{:}, imageTypeHyp, obsPosXY, probeGridXY);  %obsNum-by-probeNum
           gaussMix.means(:,imageTypeInd) = btweenPixAndProbesCov' * alpha;
           V = L'\btweenPixAndProbesCov;  % named by C. Rausmussen
           gaussMix.variances(:,imageTypeInd) = max(0, ones(probeNum,1) - sum(V.*V,1)' + expectedProbesVariance );
        end
    end
    varargout{1} = logTypeLikelihoods;        
    if (strcmp(mode, 'distribution'))
        varargout{2} = gaussMix;
    end    
end

function imageTypeDef = defineImageType()
    expConst = defineExperimentConstants();
    overallVariance = 1;
    xScalePA = expConst.imageLengthCm/20;
    yScalePA = expConst.imageLengthCm/20;
    xScaleSH = expConst.imageLengthCm/6;
    yScaleSH = expConst.imageLengthCm/30;
    xScaleSV = expConst.imageLengthCm/30;
    yScaleSV = expConst.imageLengthCm/6;
    imageTypeDef.hyp(1).def = [log(xScalePA), log(yScalePA), log(overallVariance)];
    imageTypeDef.hyp(2).def = [log(xScaleSH), log(yScaleSH), log(overallVariance)];
    imageTypeDef.hyp(3).def = [log(xScaleSV), log(yScaleSV), log(overallVariance)];
    imageTypeDef.covfunc = {@covSEard};
    imageTypeDef.typeNum = 3;
end

function displayedPix = calDisplayedPix(actualPix, dispRandn, displayedPixVariance)
    displayedPix = actualPix + dispRandn.*sqrt(displayedPixVariance);
end
% group end

function logTypePosteriors = calLogTypePosteriors(logTypeLikelihoods, decisParam)
    priorTypeVector = decisParam.priorTypeVector;
    maxLogLikelihood = max(logTypeLikelihoods);
    logLikelihoodDiffs = logTypeLikelihoods - maxLogLikelihood;
    logNormalization = log( sum( priorTypeVector.*exp(logLikelihoodDiffs) ) );    
    logTypePosteriors = log(priorTypeVector) + logLikelihoodDiffs - logNormalization;
end

function [probChooseCategP, logPosteriorCategP] = calChoiceProb(logTypePosteriors, decisParam)
    softmaxRate = decisParam.softmaxRate;
    lapseRate = decisParam.lapseRate;
    logPosteriorCategP = logTypePosteriors(1);
    logPosteriorCategS = log( exp(logTypePosteriors(2)) + exp(logTypePosteriors(3)) );
    logPosteriorRatio = logPosteriorCategP - logPosteriorCategS;  % for choosing Type P
    softenedLikelihood = 1 / ( 1+exp(-softmaxRate*logPosteriorRatio) );
    probChooseCategP = (1-lapseRate)*softenedLikelihood + 0.5*lapseRate;    
end

% findPosWithMaxScore group start
function [newObs, scoreMapLinearized] = findPosWithMaxScore(gaussMix, probeGridXY, obs, param)
% ToDo variations: score can also be entropy
% ToDo variations: score can be calculated by brute force
% ToDo variations: can find pos with softmax 
    scoreMapLinearized = calScoreMap(gaussMix);
    moveCostMap = calMoveCostMap(probeGridXY, obs, param);
    costMapLinearized = scoreMapLinearized + moveCostMap; % bad naming convension cost vs score
    selectionProb = transformHardMaxScore2Prob(costMapLinearized);
    selectionInd = getIndByProbMatching(selectionProb);
    newObs.posXY = probeGridXY(selectionInd, :);
    newObs.gridInd = selectionInd;
end

function scoreMapLinearized = calScoreMap(gaussMix)
    boundeLogWeights = lowerBoundLogWeights(gaussMix.logWeights, -30);
    means = reshape( gaussMix.means, numel(gaussMix.means), 1);
    variances = reshape( gaussMix.variances, numel(gaussMix.variances), 1);
    scoreMapLinearized = a_GetBALDGM_mex(boundeLogWeights, means, variances);    
    scoreMapLinearized = max(0, scoreMapLinearized);
end

function moveCostMap = calMoveCostMap(probeGridXY, obs, param)
    obsNum = obs.num;
    lastX = obs.posXY(obsNum, 1);
    lastY = obs.posXY(obsNum, 2);
    ind = probeGridXY(:,1) == lastX & probeGridXY(:,2) == lastY;
    moveCostMap = ind*param.move.cost;
end

function boundedLogWeights = lowerBoundLogWeights(logWeights, lowerBound)
    if( min(logWeights) < lowerBound)
        boundedLogWeights = max(logWeights, lowerBound);
        logNormalization = log( sum(exp(boundedLogWeights)) );
        boundedLogWeights = boundedLogWeights - logNormalization;
    else
        boundedLogWeights = logWeights;
    end
end

function prob = transformHardMaxScore2Prob(scoreMapLinearized)
    maxValue = max(scoreMapLinearized);
    prob = zeros(numel(scoreMapLinearized),1);
    prob( scoreMapLinearized == maxValue ) = 1;
    prob = prob/sum(prob);
end

function index = getIndByProbMatching(prob)  
    cumProb = cumsum(prob);
    index = find(cumProb > rand(), 1, 'first');    
end
% group end

function figureHandle = setupSimuPlot()
    scrsz = get(0,'ScreenSize');
    figureHandle = figure('Name', 'Simu Plot Window', 'position',[scrsz(3)/10 scrsz(4)/4 scrsz(3)/1.2 scrsz(4)/2],'menubar','none');
end

% plotSimuStep group start
function plotSimuStep(obsPosXY, actualImage, scoreMapLinearized)        
    obsPosX = obsPosXY(:,1);
    obsPosY = obsPosXY(:,2);
    rowNum = sqrt( numel(scoreMapLinearized) );
    expConst = defineExperimentConstants();
    colormap(gray);

    linSpaceCm = linspace(-expConst.imageLengthCm/2, expConst.imageLengthCm/2, rowNum) ;
    subplot(1,2,1);
    imagesc(linSpaceCm, linSpaceCm, reshape(scoreMapLinearized, rowNum, rowNum));
    addPlotObsPos(obsPosX, obsPosY);
    axis square;
       
    imageLinSpaceCm = linspace( -expConst.imageLengthCm/2, expConst.imageLengthCm/2, size(actualImage,1) ) ;
    subplot(1,2,2);
    imagesc(imageLinSpaceCm, imageLinSpaceCm, actualImage);
    addPlotObsPos(obsPosX, obsPosY);
    axis square;
    
    drawnow
    pause
end

function addPlotObsPos(obsPosX, obsPosY)
    hold on;         
    plot(obsPosX, obsPosY, 'r.');
    plot(obsPosX(end), obsPosY(end), 'ro','MarkerFaceColor','g', 'MarkerSize',5);
    hold off;
end
% group end

% plotPosIndVSTime group start
function plotPosIndVsTime(recordObs, colorSpec)
    trialNum = numel(recordObs);
    for trialInd = 1:trialNum
        trialObs = recordObs(trialInd);
        [~, uniqueInd, ~] = tallyTimeWRTUnique(trialObs, 'rawPosXY');
        time = 1:trialObs.time;
        plotTimeDuration(uniqueInd, time, trialInd/trialNum/2, colorSpec);
    end
end

function [uniqueField, uniqueInd, timeDuration] = tallyTimeWRTUnique(trialObs, fieldname)
    field = trialObs.(fieldname);
    [uniqueField, ~, uniqueInd] = unique(field, 'rows', 'stable');
    timeDuration = zeros(numel(uniqueField), 1);
    for index = 1:numel(uniqueInd)
        timeDuration(uniqueInd(index)) = timeDuration(uniqueInd(index)) + 1;
    end
end

function plotTimeDuration(aField, time, fieldShift, colorSpec)
    entryNum = numel(aField);
    for ind = 1:entryNum
        if (ind < entryNum && aField(ind) == aField(ind+1))
            plot (aField(ind:ind+1) + fieldShift, time(ind:ind+1), '-', 'LineWidth', 1, 'Color', colorSpec);
        else
            plot (aField(ind) + fieldShift, time(ind), '.', 'MarkerSize', 1, 'Color', colorSpec);
        end
        hold on;
    end
end
% group end

% plotDuration group start
function plotDuration(simuObs, lineSpec, timeProbe)
    if (nargin == 2)
        mode = 'fixationIndex';
    elseif ( nargin == 3)
        mode = 'timeProbe';
    end
    trialNum = numel(simuObs);
    if (strcmp(mode, 'fixationIndex'))
        obsTime = simuObs(1).time; % assumes shape of simu(1).obsNum is constant
        fixationDuration = nan(trialNum, obsTime);
    elseif (strcmp(mode, 'timeProbe'))
        colNum = numel(timeProbe);
        fixationDurationProbe = nan(trialNum, colNum);
    end
    for trialInd = 1:trialNum
        trialObs = simuObs(trialInd);
        fixationInd = calFixationInd(trialObs);
        [trialFixationDuration, fixationNum] = calFixationDuration(fixationInd);
        if (strcmp(mode, 'fixationIndex')) 
            fixationDuration(trialInd, 1:fixationNum) = trialFixationDuration;
        elseif (strcmp(mode, 'timeProbe'))
            trialFixDurationProbe = calFixationDurationProbe(fixationInd, trialFixationDuration, timeProbe);
            fixationDurationProbe(trialInd, :) = trialFixDurationProbe;
        end
    end
    if (strcmp(mode, 'fixationIndex'))
        plotStatMat(fixationDuration, lineSpec);
    elseif (strcmp(mode, 'timeProbe'))
        plotStatMat(fixationDurationProbe, lineSpec);
    end
end

function fixationInd = calFixationInd(trialObs)
    rawGridInd = trialObs.rawGridInd;
    obsTime = trialObs.time;
    fixationInd = ones(obsTime, 1);
    for timeInd = 1:obsTime-1
       if ( rawGridInd(timeInd) == rawGridInd(timeInd+1))
           fixationInd(timeInd+1) = fixationInd(timeInd);
       else
           fixationInd(timeInd+1) = fixationInd(timeInd) + 1;
       end
    end
%     fixationInd
%     pause
end

function [fixationDuration, fixationNum] = calFixationDuration(fixationInd)
    fixationNum = max(fixationInd);
    fixationDuration = nan(fixationNum, 1);
    for ind = 1:fixationNum
        fixationDuration(ind) = sum(fixationInd == ind);
    end
%     fixationDuration
%     pause
end

function fixDurationProbe = calFixationDurationProbe(fixInd, fixDuration, timeProbe)
    timeProbeInd = floor(timeProbe);
    fixProbInd = fixInd(timeProbeInd);
    fixDurationProbe = fixDuration(fixProbInd);
%     fixDurationProbe
%     pause
end
% group end

% plotTimeVsUtilField group start
function plotTimeVsUtilField(utilMeas, field, logIt, lineSpec)
    [trialNum, colNum] = calField2MatSize(utilMeas, field);
    simuField = nan(trialNum, colNum);
    for trialInd = 1:trialNum
        simuField(trialInd, :) = utilMeas(trialInd).(field);
    end
    if (logIt == 1)
        plotStatMat(log(simuField), lineSpec);
    else
        plotStatMat(simuField, lineSpec);
    end
end

function [rowNum, colNum] = calField2MatSize(utilMeas, field)
    rowNum = numel(utilMeas);
    colNum = numel(utilMeas(1).(field)); % assumes shape of utilMeas(n).field is constant
end
% group end

% plotStatMat group start
function plotStatMat(mat, lineSpec)
% ToAdd: 25-50-75% tile fill plot
    %plotPercentiles(mat); % problematic
    plotAvgSem(mat, lineSpec);
end

function plotPercentiles(mat)
    [rowN, colN] = size(mat);
    line75Per = zeros(1, colN);
    line50Per = zeros(1, colN);
    line25Per = zeros(1, colN);
    for ind = 1:colN
        col = mat(:,ind);
        col = sort(col);       
        line75Per(ind) = col(round(rowN*0.75));
        line50Per(ind) = col(round(rowN*0.50));
        line25Per(ind) = col(round(rowN*0.25));
    end
    fill([1:1:colN, colN:-1:1], [line75Per, fliplr(line50Per)], 'r', 'FaceAlpha', 0.3);
    hold on;
    fill([1:1:colN, colN:-1:1], [line50Per, fliplr(line25Per)], 'b', 'FaceAlpha', 0.3);
    hold on;
end

function plotAvgSem(mat, lineSpec)
    colN = size(mat, 2);
    avgVec = nanmean(mat, 1);
    varVec = nanvar(mat, 1);
    numVec = sum(1-isnan(mat), 1);
    semVec = sqrt(varVec./numVec);
    errorbar(1:colN, avgVec, semVec, lineSpec);
    hold on;
end
% group end
