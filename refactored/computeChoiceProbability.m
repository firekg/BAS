%%
function probChooseCategP = computeChoiceProbability(observations, param)
    perceptParam = param.percept;
    decisionParam = param.decision;
    perceivedPixelsVariance = calVarianceOnPoints(observations, perceptParam);
    logTypeLikelihoods = calGPPredictions(observations, perceivedPixelsVariance);    
    logTypePosteriors = calLogTypePosteriors(logTypeLikelihoods, decisionParam);
    [probChooseCategP, ~] = calChoiceProbability(logTypePosteriors, decisionParam);
end