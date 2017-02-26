function RevPosM = Get_RevPosM_cm(nRev,RevID)

Type=RevID(1);
X=RevID(2);
Y=RevID(3);
Num=RevID(4);

y=imread(strcat('for exp/Type', num2str(Type),'X',num2str(X),'Y',num2str(Y),'Num', num2str(Num), '.bmp'));
y=cast(y,'double');  yimage=(y(:,:,1))/255*10-5;  %pixel value max=255; mode of extreme value~5;

BALDxD = zeros(nRev,2);

[npix,ImageSize,~,w,pin]=Set_para_in();
BALDxD(1,:) = [randn(1),randn(1)]*0.2*npix/ImageSize + npix/2 + 1;  % eye movement noise~0.3cm; scaled from cm to pixel

count=1;
while ( count<nRev )
    
    count=count+1;
    
    % BALD searcher
    [~, ~, H_BALD, z] = Get_BALDmap_window(BALDxD(1:count-1,:), yimage, pin,w);
    [~, ind] = max(H_BALD);
    BALDxnew = z(ind,:)+[randn(1),randn(1)]*0.2*npix/ImageSize;
    BALDxnew = min(npix-w+1,BALDxnew);  BALDxnew = max(1+w-1,BALDxnew);  
    BALDxD(count,:) = BALDxnew;

end

% This is the size of RevPosM: zeros(1,nRev*2);
RevPosM=reshape(BALDxD',1,nRev*2);  %[x1,y1,x2,y2,...,xn,yn]
RevPosM=(RevPosM - 1 - npix/2)*ImageSize/npix;