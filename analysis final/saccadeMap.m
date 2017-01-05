function Dout = saccadeMap(Din)
%% for one heuristic

% load('SCdata.mat', 'SCAL');
% D = SCAL;

x = D.RevealPosX;
y = D.RevealPosY;
x(x==0) = nan;
y(y==0) = nan;
dx = diff(x,1,2);
dy = diff(y,1,2);

% plot(dx(:),dy(:),'.');
% plot(dx(:),dy(:),'.');

% pa = SCAL.ImageID(:,2)==20;
% sh = SCAL.ImageID(:,2)==6;
% sv = SCAL.ImageID(:,2)==30;

% plot(dx(pa,:),dy(pa,:),'.');
% plot(dx(sh,:),dy(sh,:),'.');
% plot(dx(sv,:),dy(sv,:),'.');

Dout = Din;
Dout.SaccadeMapX = dx;
Dout.SaccadeMapY = dy;

end
