% easy
% pa=[20,20]; sh=[6,30]; sv=[30,6];

% medium easy
% pa=[10,10]; sh=[6,15]; sv=[15,6];

% medium
pa=[7,7]; sh=[4,10]; sv=[4,10];

% sort of hard
% pa=[5,5]; sh=[3,7]; sv=[7,3];

% all small -- might not be bad
% pa=[20,20]; sh=[15,30]; sv=[15,30];

p=0.9;  v=4;  b=-log(1/p-1)/v;
sigmoid=@(x) 0.5*(1+exp(-b*x)).^(-1);

npix = 770;
ngppix = npix/10;

MASKOn=0;
appersize=30;  %in pixels
mu = [npix npix]./2;
Sigma = [appersize^2, 0; 0, appersize^2];
x1 = 1:npix; x2 = 1:npix;
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma); 
F = reshape(F,length(x2),length(x1));
MASK = mat2gray(repmat(F,[1,1,3]));

hf=figure(1);
set(gcf,'Position',[50,50,900,900]);
%set(gcf,'units','centimeters','PaperSize',[20 20]);
set(gcf,'Toolbar','None');
set(gcf,'MenuBar','None');

ni=10;
for i=1:ni
    clf;
    
    if(rand()>0.5)
        xsf=pa(1);  ysf=pa(2);  t=1;
    else
        if(rand()>0.5)
            xsf=sh(1);  ysf=sh(2);  t=2;
        else
            xsf=sv(1);  ysf=sv(2);  t=2;
        end
    end
    
    disp(sprintf('Trial=%d. Type=%d;',i, t));
        
    % Initialize GP image meshgrid
    px = linspace(1, npix, ngppix);
    py = linspace(1, npix, ngppix);
    [x1image, x2image] = meshgrid(px,py);
    ximage = [reshape(x1image,ngppix^2,1), reshape(x2image,ngppix^2,1)];
    
    % Initialize GP hyperparameter
    covfunc = {@covSEard};
    likfunc = {@likGauss};
    hyp.cov = [log(npix/xsf), log(npix/ysf), log(1)];  % determines what name
    hyp.lik = log(1e-2);
    
    % Initialize Chol decomposition
    K = feval(covfunc{:}, hyp.cov, ximage) + 1e-12*eye(ngppix^2);  % added 1e-6 noise to allow chol
        
    yimage = chol(K)'*randn(ngppix^2,1);
    yimage = reshape(yimage,ngppix,ngppix);

% % Initialize interpolation meshgrid
%     pxI = linspace(1, npix, npix);
%     pyI = linspace(1, npix, npix);
%     [XI, YI] = meshgrid(pxI,pyI);
% 
%     ZI = interp2(px,py,yimage,XI,YI, 'spline');
%     ZI = ZI + exp(hyp.lik)*randn(npix,npix);  % add hyp.lik noise after spline
%     while(max(max(abs(ZI)))>4)  % to make sure we get constrained samples
%         yimage = chol(K)'*randn(ngppix^2,1);
%         yimage = reshape(yimage,ngppix,ngppix);
%         ZI = interp2(px,py,yimage,XI,YI, 'spline');
%         ZI = ZI + exp(hyp.lik)*randn(npix,npix);  % add hyp.lik noise after spline
%     end
%     
%     %red/blue
%     IR=mat2gray(ZI,[-4,4]);  % scaled for image display
%     IG=IR*0;
%     IB=1-IR;
%     I=zeros(npix,npix,3);
%     I(:,:,1)=IR;  I(:,:,2)=IG;  I(:,:,3)=IB;    
%     image(I);
    
    figure(1);
    set(gca,'position',[0 0 1 1],'units','normalized');

    %image(BG);  %noisy background
    if(MASKOn==1)
        image(MASK);
    end
    
    width=sigmoid(yimage);
    height=sigmoid(-yimage);
    for y=1:rows(width)
        for x=1:cols(width)
            h=rectangle('Position',10*[x,y,width(y,x),height(y,x)],'Curvature',[1 1]);
            set(h,'FaceColor','k');
        end
    end
    myaa(8,'standard');

    axis off;
        
%     filename=['img' num2str(i) '.pdf'];  export_fig(filename);  % to pdf
%     print('-dpng', sprintf('%s.png',filename), '-r1200');
%     print (gcf, '-dbmp', ['img' num2str(i) '.png']); 
    
%     % allows pressing q to exit the loop
%     if strcmp(get(hf,'currentcharacter'),'q')
%         break
%     end
%     drawnow

    pause
    close(gcf);

end

%%
% 
% width=relsize*sigmoid(yimage);
% height=relsize*sigmoid(-yimage);
% 
% figure(1);
% clf;
% hold on;
% for i=1:rows(width)
%     for j=1:rows(width)
%         h=rectangle('Position',[i,j,width(i,j),height(i,j)],'Curvature',[0 0]);
%         set(h,'FaceColor','k');
%     end
% end
% hold off;
% 
