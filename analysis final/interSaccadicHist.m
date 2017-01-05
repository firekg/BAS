function [dtHist, x] = interSaccadicHist()

load('SCdata.mat','SCAL');
load('EPdata.mat','EPAL');
load('BZdata.mat','BZAL');
revTime = [SCAL.RevealTime; EPAL.RevealTime; BZAL.RevealTime];
revTime(revTime==0) = nan;
time = revTime;
dt = diff(time, 1, 2);
x = 0.001:0.05:1.501; %following land99, 50ms bin, 0-3 s
dtHist = hist(dt(:), x);
%bar(x, dtHist);

end

