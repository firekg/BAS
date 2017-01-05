function Ds=a_Get_AnswerReal(D,num)
% type: 0=random; 1=patch; 2=stripy horizontal; 3=stripy vertical

ind=1:D.Trials;  ind(D.AnswerReal~=num)=[];
Ds.Trials=length(ind);
Ds.RevealPosX=D.RevealPosX(ind,:);
Ds.RevealPosY=D.RevealPosY(ind,:);
Ds.RevealTime=D.RevealTime(ind,:);
Ds.StateAnswerTime=D.StateAnswerTime(ind,:);
Ds.StateSearchTime=D.StateSearchTime(ind,:);
Ds.AnswerChoice=D.AnswerChoice(ind,:);
Ds.AnswerReal=D.AnswerReal(ind,:);
Ds.ImageID=D.ImageID(ind,:);
Ds.MaxRevealingTrial=D.MaxRevealingTrial(ind,:);
Ds.RevealType=D.RevealType(ind,:);

