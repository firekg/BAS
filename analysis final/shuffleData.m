function Dout = shuffleData(Din)

Dout = Din;
%shuffle revealPos X, Y, and Z together
for trial = 1:Din.Trials
    mrev = Din.MaxRevealingTrial(trial);
    perm = randperm(mrev);
    x = Din.RevealPosX(trial,1:mrev);
    y = Din.RevealPosY(trial,1:mrev);
    z = Din.RevealZ(trial,1:mrev);
    Dout.RevealPosX(trial,1:mrev) = x(perm);
    Dout.RevealPosY(trial,1:mrev) = y(perm);
    Dout.RevealZ(trial,1:mrev) = z(perm);
end

end