function a_ErrorShade(x,upp,low,colorv,a)

x=reshape(x,numel(x),1);  upp=reshape(upp,numel(upp),1);  low=reshape(low,numel(low),1);
xbox = [x; flipdim(x,1)];
ybox = [upp; flipdim(low,1)];
xbox(isnan(xbox))=[]; ybox(isnan(xbox))=[];
xbox(isnan(ybox))=[]; ybox(isnan(ybox))=[];
fill(xbox, ybox, colorv, 'FaceAlpha', a, 'LineStyle','None');

