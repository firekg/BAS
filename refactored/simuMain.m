%%
% refactoring conventions
% function name stricture: verb-adj-noun-prep
% compute is a collection of cals

function [simuObsMat, utilMeasMat, varyParam] = simuMain(imageData)    
    param = initParam();
    
%     simulateBAS(imageData, param);

    varyParam.dispNoiseSD = [0.5, 1, 1.2, 1.5, 2];
    varyParam.moveCost = [0, exp(-7)];
    simuObsMat = struct([]);
    utilMeasMat = struct([]);
    noiseLevelNum = numel(varyParam.dispNoiseSD);
    costLevelNum = numel(varyParam.moveCost);
    for par1Ind = 1:noiseLevelNum
        for par2Ind = 1:costLevelNum
            param.disp.dispNoiseSD = varyParam.dispNoiseSD(par1Ind);
            param.move.cost = varyParam.moveCost(par2Ind);
            param.simu.trialNum = 100;
            param.simu.maxObsTime = 600;
            param.simu.maxObsNum = 100;
            param.simu.maxCertainty = 0.97;
            [simuObs, utilMeas] = simulateBAS(imageData, param);
            simuObsMat(par1Ind, par2Ind).simuObs = simuObs;
            utilMeasMat(par1Ind, par2Ind).utilMeas = utilMeas;
        end
    end
end

function param = initParam()
% can include expConst and imageTypeDef
    param.disp.dispNoiseSD = 1;
    param.percept.obsNoiseSD = 0.1;
    param.percept.forgettingRate = nan;
    param.percept.accuityFallOffRate =  nan;
    param.percept.accuityRadius = nan;
    param.decis.softmaxRate = 0;
    param.decis.lapseRate = 0;
    param.decis.priorTypeVector = [0.5, 0.25, 0.25];
    param.move.cost = 0;
    param.simu.trialNum = 10;
    param.simu.maxObsTime = 30;
    param.simu.maxObsNum = 10;
    param.simu.maxCertainty = 0.97;
end

function [simuObs, utilMeasure] = simulateBAS(imageData, param)
% toDo: I think my update functions are inefficient!!
% toDo: test results against old version
% toFix: dispNoiseSD = 100, non-smooth BAS scores
% toCheck: tiny dispNoiseSD --> no returns
    if (nargout==0)
        mode = 'watch';
    elseif (nargout ==2)
        mode = 'record';
    end
    probeGridXY = prepProbeGridXY();
    trialNum = param.simu.trialNum;
    maxObsTime = param.simu.maxObsTime;
    if (strcmp(mode, 'watch'))
        figureHandle = setupSimuPlot();
    end
    for trialInd = 1:trialNum
        [newObs, actualImage] = initSimuTrial(imageData, probeGridXY, trialInd);
        obs = initObs();
        trialUtilMeas = initTrialUtilMeas(maxObsTime);
        obsInd = 0;
        while ( meetGoCriteria(obs, obsInd, trialUtilMeas, param) ) 
            obsInd = obsInd + 1;
            obs = updateObs(obs, newObs, actualImage);
            gaussMix = computeExpectedPixDistr(obs, param, probeGridXY);
            [newObs, scoreMapFlat] = findPosWithMaxScore(gaussMix, probeGridXY, obs, param);
            if (strcmp(mode, 'watch'))
                plotSimuStep([obs.posXY; newObs.posXY], actualImage, scoreMapFlat);
                if strcmp(get(figureHandle,'currentcharacter'), 'q')
                    break
                end
            elseif (strcmp(mode, 'record'))
                trialUtilMeas = updateTrialUtilMeas(trialUtilMeas, obsInd, scoreMapFlat, obs, newObs, gaussMix.logWeights, param);
            end
        end
        if (strcmp(mode, 'watch'))
            if strcmp(get(figureHandle,'currentcharacter'), 'q')
                break
            end
        elseif (strcmp(mode, 'record'))            
            simuObs(trialInd) = obs;
            utilMeasure(trialInd).maxBas = trialUtilMeas.maxBas;
            utilMeasure(trialInd).gapBas = trialUtilMeas.gapBas;
            utilMeasure(trialInd).certainty = trialUtilMeas.certainty;
        end
        fprintf('Cost %f; SD level %f; trial %d done.\n', param.move.cost, param.disp.dispNoiseSD, trialInd);
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

function trialUtilMeas = initTrialUtilMeas(maxObsTime)
    trialUtilMeas.maxBas = nan(maxObsTime, 1);
    trialUtilMeas.gapBas = nan(maxObsTime, 1);
    trialUtilMeas.certainty = nan(maxObsTime, 1);
end

function met = meetGoCriteria(obs, obsInd, trialUtilMeas, param)
    criterium1 = (obs.time <= param.simu.maxObsTime);
    criterium2 = (obs.num <= param.simu.maxObsNum);
    if (obsInd == 0)
        certainty = 0.5;
    else
        certainty = trialUtilMeas.certainty(obsInd);
    end
    criterium3 = (certainty <= param.simu.maxCertainty);
    met = (criterium1 && criterium2 && criterium3);
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
function [newObs, scoreMapFlat] = findPosWithMaxScore(gaussMix, probeGridXY, obs, param)
% ToDo variations: score can also be entropy
% ToDo variations: score can be calculated by brute force
% ToDo variations: can find pos with softmax 
    scoreMapFlat = calScoreMap(gaussMix);
    moveCostMap = calMoveCostMap(probeGridXY, obs, param);
    costMapLin = scoreMapFlat + moveCostMap; % bad naming convension cost vs score
    selectionProb = transformHardMaxScore2Prob(costMapLin);
    selectionInd = getIndByProbMatching(selectionProb);
    newObs.posXY = probeGridXY(selectionInd, :);
    newObs.gridInd = selectionInd;
end

function scoreMapFlat = calScoreMap(gaussMix)
    boundeLogWeights = lowerBoundLogWeights(gaussMix.logWeights, -30);
    means = reshape( gaussMix.means, numel(gaussMix.means), 1);
    variances = reshape( gaussMix.variances, numel(gaussMix.variances), 1);
    scoreMapFlat = a_GetBALDGM_mex(boundeLogWeights, means, variances);    
    scoreMapFlat = max(0, scoreMapFlat);
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

function prob = transformHardMaxScore2Prob(scoreMapFlat)
    maxValue = max(scoreMapFlat);
    prob = zeros(numel(scoreMapFlat),1);
    prob(scoreMapFlat == maxValue) = 1;
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
function plotSimuStep(obsPosXY, actualImage, scoreMapFlat)        
    obsPosX = obsPosXY(:,1);
    obsPosY = obsPosXY(:,2);
    rowNum = sqrt( numel(scoreMapFlat) );
    expConst = defineExperimentConstants();
    colormap(gray);

    linSpaceCm = linspace(-expConst.imageLengthCm/2, expConst.imageLengthCm/2, rowNum) ;
    subplot(1,2,1);
    imagesc(linSpaceCm, linSpaceCm, reshape(scoreMapFlat, rowNum, rowNum));
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

function trialUtilMeasOut = updateTrialUtilMeas(trialUtilMeasIn, obsInd, scoreMapFlat, obs, newObs, logWeights, param) % too many arguments!
    trialUtilMeasOut = trialUtilMeasIn;
    trialUtilMeasOut.maxBas(obsInd) = max(scoreMapFlat); % vs scoreMapFlat(newObs.gridInd) 
    trialUtilMeasOut.gapBas(obsInd) = max(scoreMapFlat) - scoreMapFlat(obs.rawGridInd(obsInd));
    [~, logPosteriorCategP] = calChoiceProb(logWeights, param.decis);
    trialUtilMeasOut.certainty(obsInd) = max( 1-exp(logPosteriorCategP), exp(logPosteriorCategP) );
end


