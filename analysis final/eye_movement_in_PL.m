%%
[BZPL1,~,~] = DATAFILE_Read('data\BZ20130910PLM1.DAT'); %BZPL1.FrameData=[];
[BZPL2,~,~] = DATAFILE_Read('data\BZ20130910PLM2.DAT'); %BZPL2.FrameData=[];
[BZPL3,~,~] = DATAFILE_Read('data\BZ20130910PLM3.DAT'); %BZPL3.FrameData=[];
[BZPL4,~,~] = DATAFILE_Read('data\BZ20130910PLM4.DAT'); %BZPL4.FrameData=[];
[BZPL5,~,~] = DATAFILE_Read('data\BZ20130910PLM5.DAT'); %BZPL5.FrameData=[];
[BZPL6,~,~] = DATAFILE_Read('data\BZ20130910PLM6.DAT'); %BZPL6.FrameData=[];
[BZPL7,~,~] = DATAFILE_Read('data\BZ20130910PLM7.DAT'); %BZPL7.FrameData=[];
[BZPL8,~,~] = DATAFILE_Read('data\BZ20130910PLM8.DAT'); %BZPL8.FrameData=[];
%%
 
Data = BZPL2;
for Tr = 1:100;

TarX = Data.RevealPosX;
TarY = Data.RevealPosY;
MaxRev = Data.MaxRevealingTrial;

Eye = Data.FrameData.EyeTrackerEyeXY;
%Eye(abs(Eye)>10) = nan;
EyeX = Eye(:,:,1);
EyeY = Eye(:,:,2);

clf
plot(EyeX(Tr, :), EyeY(Tr, :), '.-'); hold on;
plot(EyeX(Tr, 1:1000), EyeY(Tr, 1:1000), '.y');
plot(TarX(Tr, 1:2), TarY(Tr, 1:2), 'ro', 'MarkerFaceColor', 'r');
plot(TarX(Tr, 3:MaxRev(Tr)), TarY(Tr, 3:MaxRev(Tr)), 'go');
axis([-10 10 -10 10])

pause
end

% Don't see any correlation between revealing locations and eye
% movements....