function [varargout] = a_Get_RevealMap(D,nvspec,mode)
% mode=1: takes D; outputs revMap
% mode=2: takes D and nvspec; outputs nv-specific revMap 
% mode=3: takes D and nvspec; outputs nv-specific x-y covariance
% mode=4: takes D and nvspec; outputs binned revMap 

npix=770; ImageSize=20;

if (mode<=2)
    % creating gaussian filter
    var=50^2; %var=(npix/ImageSize/2)^2;
    n=sqrt(var)*5;
    [x1,x2]=meshgrid(1:n,1:n);
    z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % in pixel coordinate
    GaussFilter = zeros(n*n,1);
    for i=1:n*n
        GaussFilter(i) = mvnpdf(z(i,:),[1,1]*n/2,[1,0;0,1]*var);
    end
    GaussFilter = reshape(GaussFilter,n,n);
    
    revMap=zeros(npix,npix);
    count = 0;
    for trial=1:D.Trials    
        nv = D.MaxRevealingTrial(trial);
        xDx = 1+(npix-1)*(D.RevealPosX(trial,1:nv)/ImageSize+0.5);  % cm to pixel: from inverting last two lines of a_GET_BALD
        xDy = 1-(npix-1)*(D.RevealPosY(trial,1:nv)/ImageSize-0.5);
        xD = round([xDy',xDx']); % x is col 2, y is col 1; plotting convention, opposite of GP computation convention
        if (mode==1)
            for i=1:nv
                revMap(xD(i,1),xD(i,2)) = revMap(xD(i,1),xD(i,2)) + 1;
            end
        end
        if (mode==2)
            if(nv>=nvspec)
                i=nvspec;
                revMap(xD(i,1),xD(i,2)) = revMap(xD(i,1),xD(i,2)) + 1;
            end
        end
        count = count+nv;
    end

    revMap = conv2(revMap,GaussFilter,'same');
    revMap = revMap/sum(sum(revMap));
    varargout={revMap};
end

if (mode==3)
    xD=[];
    for trial=1:D.Trials    
        nv = D.MaxRevealingTrial(trial);
        if(nv>=nvspec)
            i=nvspec;
            xDx = 1+(npix-1)*(D.RevealPosX(trial,i)/ImageSize+0.5);  % cm to pixel
            xDy = 1-(npix-1)*(D.RevealPosY(trial,i)/ImageSize-0.5);
            xDa = [xDx,xDy];  % x in col 1; y in col 2; GP computation convention
            xD = [xD; xDa ]; 
        end
    end
    covM=cov(xD);
    mV=mean(xD);
    varargout={mV,covM};
end


