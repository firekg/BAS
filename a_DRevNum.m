function Di = a_DRevNum(D,revNum)
% type: 20=patch; 6=stripy horizontal; 30=stripy vertical

ind=1:D.Trials;
ind(D.MaxRevealingTrial~=revNum)=[];

Di.Trials=length(ind);
Di.RevealPosX=D.RevealPosX(ind,:);
Di.RevealPosY=D.RevealPosY(ind,:);
Di.RevealZ=D.RevealZ(ind,:);
Di.RevealTime=D.RevealTime(ind,:);
Di.StateAnswerTime=D.StateAnswerTime(ind,:);
Di.StateSearchTime=D.StateSearchTime(ind,:);
Di.AnswerChoice=D.AnswerChoice(ind,:);
Di.AnswerReal=D.AnswerReal(ind,:);
Di.ImageID=D.ImageID(ind,:);
Di.MaxRevealingTrial=D.MaxRevealingTrial(ind,:);

end