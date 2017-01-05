%% prepare data
load('BZdata.mat','BZAL');
D = BZAL;
%D = SIM_idealBAS_PBwPix(1,1);

x = D.RevealPosX;
y = D.RevealPosY;
nT = D.Trials;
nrVec = D.MaxRevealingTrial;
dis = zeros(sum(nrVec-1),1);
ang = zeros(sum(nrVec-1),1);
rec = zeros(sum(nrVec-1),4);

% get saccade distance (dis) and angle (ang)
count = 1;
for trial = 1:nT
    for rev = 1:nrVec(trial)-1
        dx = x(trial, rev+1) - x(trial, rev);
        dy = y(trial, rev+1) - y(trial, rev);
        dis(count) = sqrt(dx^2+dy^2)*27.8/20; %cm to degree
        ang(count) = atan2(dy,dx)/pi*180;
        rec(count,1) = trial;
        rec(count,2) = rev+1;
        rec(count,3) = x(trial, rev+1);
        rec(count,4) = y(trial, rev+1);
        count = count + 1;
    end
end

%% where in "time" and space are the small saccades?
% small columns: trial; rev; rev x; rev y; dis
small_ind = (dis<0.59);
small = rec(small_ind,:);
small(:,5) = dis(small_ind);

% Nothing special jumps out in terms of time and space.
% How can this happen?
% 

%% May not be a problem now. I was using the wrong cond (cond=2)