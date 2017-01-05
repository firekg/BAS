%% load data
[D,~,~] = DATAFILE_Read('data\ALDCSY20150901a.DAT');

%% plot a trial
trial = 1;

nrev = D.MaxRevealingTrial(trial);
revX = D.RevealPosX(trial,1:nrev);
revY = D.RevealPosY(trial,1:nrev);

% eye movement in search phase
Neye = floor(D.StateSearchTime(trial));
eyeX = D.FrameData.EyeTrackerEyeXY(trial,1:Neye,1);
eyeY = D.FrameData.EyeTrackerEyeXY(trial,1:Neye,2);

% eye movement in answer phase
Nans = Neye + 1 + floor(D.StateAnswerTime(trial));
ansX = D.FrameData.EyeTrackerEyeXY(trial, Neye+1:Nans, 1);
ansY = D.FrameData.EyeTrackerEyeXY(trial, Neye+1:Nans, 2);

clf;
hold on;
plot(ansX, ansY, '.', 'color' , 0.5*[1,1,1]);
plot(eyeX, eyeY, 'b.');
plot(eyeX(1), eyeY(1), 'g.');
plot(revX, revY, 'ro');

% border
plot([-10, 10], [10, 10], 'k');
plot([-10, 10], [-10, -10], 'k');
plot([10, 10], [-10, 10], 'k');
plot([-10, -10], [-10, 10], 'k');

axis([-15, 15, -15, 15]);
axis square;
