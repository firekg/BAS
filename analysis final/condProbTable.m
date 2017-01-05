function T = condProbTable(D)
pa = D.ImageID(:,2)==20;
sh = D.ImageID(:,2)==6;
sv = D.ImageID(:,2)==30;

winThetaPa = D.winner(pa,:);
winThetaPa = winThetaPa(:);
nTotThetaPa = sum(~isnan(winThetaPa));
pPaPa = sum(winThetaPa==1)/nTotThetaPa;
pShPa = sum(winThetaPa==2)/nTotThetaPa;
pSvPa = sum(winThetaPa==3)/nTotThetaPa;

winThetaSh = D.winner(sh,:);
winThetaSh = winThetaSh(:);
nTotThetaSh = sum(~isnan(winThetaSh));
pPaSh = sum(winThetaSh==1)/nTotThetaSh;
pShSh = sum(winThetaSh==2)/nTotThetaSh;
pSvSh = sum(winThetaSh==3)/nTotThetaSh;

winThetaSv = D.winner(sv,:);
winThetaSv = winThetaSv(:);
nTotThetaSv = sum(~isnan(winThetaSv));
pPaSv = sum(winThetaSv==1)/nTotThetaSv;
pShSv = sum(winThetaSv==2)/nTotThetaSv;
pSvSv = sum(winThetaSv==3)/nTotThetaSv;

T = [ [pPaPa, pShPa, pSvPa];...
      [pPaSh, pShSh, pSvSh];...
      [pPaSv, pShSv, pSvSv]];
% the frist is the inferred type, second is the true type
  
end
