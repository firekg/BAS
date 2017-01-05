function Di=a_Dtrials(D,tr)
% tr = indices of trials

Di.Trials=length(tr);
Di.RevealPosX=D.RevealPosX(tr,:);
Di.RevealPosY=D.RevealPosY(tr,:);
Di.RevealTime=D.RevealTime(tr,:);
Di.StateAnswerTime=D.StateAnswerTime(tr,:);
Di.StateSearchTime=D.StateSearchTime(tr,:);
Di.AnswerChoice=D.AnswerChoice(tr,:);
Di.AnswerReal=D.AnswerReal(tr,:);
Di.ImageID=D.ImageID(tr,:);
Di.MaxRevealingTrial=D.MaxRevealingTrial(tr,:);