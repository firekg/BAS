function PLt=a_Get_PL_of_RevealType(D,type)
% type: 0=random; 1=patch; 2=stripy horizontal; 3=stripy vertical

ind=1:D.Trials;  ind(D.RevealType~=type)=[];
PLt.Trials=length(ind);
PLt.RevealPosX=D.RevealPosX(ind,:);
PLt.RevealPosY=D.RevealPosY(ind,:);
PLt.RevealTime=D.RevealTime(ind,:);
PLt.StateAnswerTime=D.StateAnswerTime(ind,:);
PLt.StateSearchTime=D.StateSearchTime(ind,:);
PLt.AnswerChoice=D.AnswerChoice(ind,:);
PLt.AnswerReal=D.AnswerReal(ind,:);
PLt.ImageID=D.ImageID(ind,:);
PLt.MaxRevealingTrial=D.MaxRevealingTrial(ind,:);
PLt.RevealType=D.RevealType(ind,:);

