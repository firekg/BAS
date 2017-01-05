function Di = removeTrials(D,rmInd)
% yes = 1: get correct trials
% yes = 0: get incorrect trials

ind = 1:D.Trials;
ind(rmInd) = [];

Di.Trials=length(ind);
Di.RevealPosX=D.RevealPosX(ind,:);
Di.RevealPosY=D.RevealPosY(ind,:);
Di.RevealTime=D.RevealTime(ind,:);
Di.StateAnswerTime=D.StateAnswerTime(ind,:);
Di.StateSearchTime=D.StateSearchTime(ind,:);
Di.AnswerChoice=D.AnswerChoice(ind,:);
Di.AnswerReal=D.AnswerReal(ind,:);
Di.ImageID=D.ImageID(ind,:);
Di.MaxRevealingTrial=D.MaxRevealingTrial(ind,:);

[nrow, ~] = size(Di.AnswerChoice);
fprintf('trial %d;\n', nrow);
end