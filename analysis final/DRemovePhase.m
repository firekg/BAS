function Dout = DRemovePhase(Din)

Dout = Din;

ntr = Din.Trials;
for  itr = 1:ntr
    nr = Din.MaxRevealingTrial(itr);
    x = Din.RevealPosX(itr,1:nr);
    y = Din.RevealPosY(itr,1:nr);
    xmu = mean(x);
    ymu = mean(y);
    xrp = x-xmu;
    yrp = y-ymu;
    %xrp(xrp<-10 | xrp>10) = 0;
    %yrp(yrp<-10 | yrp>10) = 0;    
    Dout.RevealPosX(itr,1:nr) = xrp;
    Dout.RevealPosY(itr,1:nr) = yrp;    
end

fprintf('max x = %f\n', max(Din.RevealPosX(:)));
fprintf('min x = %f\n', min(Din.RevealPosX(:)));
fprintf('max y = %f\n', max(Din.RevealPosY(:)));
fprintf('min y = %f\n', min(Din.RevealPosY(:)));

fprintf('max x = %f\n', max(Dout.RevealPosX(:)));
fprintf('min x = %f\n', min(Dout.RevealPosX(:)));
fprintf('max y = %f\n', max(Dout.RevealPosY(:)));
fprintf('min y = %f\n', min(Dout.RevealPosY(:)));

end