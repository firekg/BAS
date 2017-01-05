function [Di, rmInd]  = correctTrials(D,yes)
% yes = 1: get correct trials
% yes = 0: get incorrect trials

ind = 1:D.Trials;
if (yes==1)
    rmInd = D.AnswerChoice~=D.AnswerReal;
    ind(rmInd) = [];
elseif (yes==0)
    rmInd = D.AnswerChoice==D.AnswerReal;
    ind(rmInd) = [];
end

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
c = sum(Di.AnswerChoice == Di.AnswerReal);
fprintf('trial %d; correct: %d\n', nrow, c);
end