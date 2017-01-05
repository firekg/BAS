function [ZI, ZI_sub]=b_GenGPimage(xsf,ysf)
%% generate GP images
% type = image type, for naming
% xsf = x scale factor in npix/xsf
% ysf = y scale factor in npix/ysf
% istart = image number start, for naming
% iend = image number end, for naming

npix = 770;
ngppix = npix/5;

% Initialize GP image meshgrid
px = linspace(1, npix, ngppix);
py = linspace(1, npix, ngppix);
[x1image, x2image] = meshgrid(px,py);
ximage = [reshape(x1image, ngppix^2,1), reshape(x2image, ngppix^2,1)];

px2 = linspace(1, npix, ngppix/2);
py2 = linspace(1, npix, ngppix/2);

% Initialize GP hyperparameter
covfunc = {@covSEard};
likfunc = {@likGauss};
hyp.cov = [log(npix/xsf), log(npix/ysf), log(1)];  % determines what name
hyp.lik = log(1e-2);

% Initialize Chol decomposition
K = feval(covfunc{:}, hyp.cov, ximage) + 1e-8*eye(ngppix^2);  % added 1e-6 noise to allow chol

% Initialize interpolation meshgrid
pxI = linspace(1, npix, npix);
pyI = linspace(1, npix, npix);
[XI, YI] = meshgrid(pxI,pyI);

I=zeros(npix,npix,3);
% for i=istart:iend
    yimage = chol(K)'*randn(ngppix^2,1);
    yimage = reshape(yimage,ngppix,ngppix);
    yimage_sub = imresize(yimage,1/2);
    ZI = interp2(px, py, yimage, XI, YI, 'spline');
    ZI_sub = interp2(px2, py2, yimage_sub, XI, YI, 'spline');
    temp = exp(hyp.lik)*randn(npix,npix);
    %ZI = ZI + temp;
    %ZI_sub = ZI_sub + temp;
    
%     ZI = ZI + exp(hyp.lik)*randn(npix,npix);  % add hyp.lik noise after spline
%     while(max(max(abs(ZI)))>4)  % to make sure we get constrained samples
%         yimage = chol(K)'*randn(ngppix^2,1);
%         yimage = reshape(yimage,ngppix,ngppix);
%         ZI = interp2(px,py,yimage,XI,YI, 'spline');
%         %ZI = ZI + exp(hyp.lik)*randn(npix,npix);  % add hyp.lik noise after spline
%     end
 
%     % save the actual values
%     dlmwrite(strcat('im/Type',num2str(type),'X',num2str(xsf),...
%         'Y',num2str(ysf),'Num', num2str(i), 'raw.txt'), ZI, ' ');
%     
%     %red/blue
%     IR=mat2gray(ZI,[-4,4]);  % scaled for image display
%     IG=IR*0;
%     IB=1-IR;
%     I(:,:,1)=IR;  I(:,:,2)=IG;  I(:,:,3)=IB;
%     imwrite(I, strcat('im/Type',num2str(type),'X',num2str(xsf),...
%         'Y',num2str(ysf),'Num', num2str(i), '.bmp'), 'bmp');
    %     i
% end
vmax = max([max(abs(ZI(:))),max(abs(ZI_sub(:)))]);
subplot(2,2,1); image(64*(ZI/vmax/2+0.5)); axis square;
subplot(2,2,2); image(64*(ZI_sub/vmax/2+0.5)); axis square;
subplot(2,2,3); image(64*((ZI-ZI_sub)/vmax/2+0.5)); axis square;
%figure; hist(ZI(:)-ZI_sub(:), 1e5);

end