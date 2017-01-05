function Dm = a_DImgVSRev(D, imgXID, revID)
% ImgXID: 20->patch; 6->stripy horizontal; 30->stripy vertical
% RevID:  1 ->patch; 2->stripy horizontal; 3 ->stripy vertical

ind=1:D.Trials;
%imgIDRemove = D.ImageID(:,2)~=imgXID;
%revIDRemove = D.RevealType~=revID;
ind(D.ImageID(:,2)~=imgXID | D.RevealType~=revID)=[];
Dm.Trials=length(ind);
Dm.RevealPosX=D.RevealPosX(ind,:);
Dm.RevealPosY=D.RevealPosY(ind,:);
Dm.RevealTime=D.RevealTime(ind,:);
Dm.StateAnswerTime=D.StateAnswerTime(ind,:);
Dm.StateSearchTime=D.StateSearchTime(ind,:);
Dm.AnswerChoice=D.AnswerChoice(ind,:);
Dm.AnswerReal=D.AnswerReal(ind,:);
Dm.ImageID=D.ImageID(ind,:);
Dm.MaxRevealingTrial=D.MaxRevealingTrial(ind,:);
Dm.RevealType=D.RevealType(ind,:);

end