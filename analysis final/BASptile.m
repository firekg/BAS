function B = BASptile()
    load('lpx.mat');
    load('expDataWPix.mat');
    
    %pv= [  lsigo,    phit,   lphis, phis0,        lbx,        lbd, lapse, prpa, window, prxpa, prCommon, shift];
    pv =[log(1.0),   nan,     nan,   nan,       log(0),        nan,     0,  1/2,      1,   1/3,      1/3,    12];
    %pw=[lbm,  km,  nm, bias];
    pw=[   0, nan, nan,    0]; % none of the last 3 are used here
    
    % make probe points z
    npix=770;
    %ImageSize=20;
    n=110;
    xx = linspace(1,npix,n);
    yy = linspace(1,npix,n);
    [x1,x2]=meshgrid(xx,yy);
    z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
    %dxx = xx(2)-xx(1);
        
    for i = 1:3 %3
        if(i==1)
            AL = SCALwPix;
            pv(1) = log(0.5);
            pv(6) = 0.34901995;
            pv(7) = 0.04471412;
            pv(12) = 16;
        elseif(i==2)  
            AL = EPALwPix;  
            pv(1) = log(0.5);
            pv(6) = 0.637376715;
            pv(7) = 0.120818875;
            pv(12) = 17;
        elseif(i==3)  
            AL = BZALwPix;  
            pv(1) = log(0.3);
            pv(6) = 0.374908933;
            pv(7) = 0.101003052;
            pv(12) = 15;
        end
        fprintf('Done loading parti %d.\n', i);
        for task = 1
            Gmap = zeros(AL.Trials,25,3); % zero makes lpx of no effect in a_Model 
            B(i,1,task) = a_Get_BALDscoreProp(AL,[],pv,pw,randn(1000,25), Gmap, z, 2.5); %BAS
            B(i,2,task) = a_Get_BALDscoreProp(AL,[],pv,pw,randn(1000,25), Gmap, z, 2.6); %MaxEnt
            fprintf('Done parti %d task %d.\n', i, task);
        end
    end

end

%% quick analysis
% 2015-09-11
load('BASptile-2015-09-11.mat');
ncomb = numel(Bcomb);
clf;
% xx = 0:1:102;

cond = 1;
%BAS
ptile = [];
for icomb = 1:ncomb
    for parti = 1:3
        temp = Bcomb(icomb).B(parti,cond).PTILExs;
        ptile = [ptile; temp];
    end
end
ptile = ptile(:);
ptile(isnan(ptile))=[];
subplot(1,2,1);
hold on;
% nh = hist(ptile, xx);
% plot(xx, nh/numel(ptile), 'k');
% [pdf,xi] = ksdensity(ptile, xx, 'width', 10, 'support', [0,100.1]);
[pdf,xi] = ksdensity(ptile, 'support', [0,101]);
% bar(xi, pdf, 'EdgeColor', 'none');
xbox = [0; xi'; 100];
ybox = [0; pdf'; 0];
fill(xbox, ybox, [0,0,1], 'FaceAlpha', 1, 'LineStyle','None');
xlabel('Score percentile');
ylabel('Probability density');
axis square
axis([0, 100, 0, 0.151]);

%MaxEnt
cond = 2;
ptile = [];
for icomb = 1:ncomb
    for parti = 1:3
        temp = Bcomb(icomb).B(parti,cond).PTILExs;
        ptile = [ptile; temp];
    end
end
ptile = ptile(:);
ptile(isnan(ptile))=[];
subplot(1,2,2);
hold on;
% nh = hist(ptile, xx);
% plot(xx, nh/numel(ptile), 'k');
% [pdf,xi] = ksdensity(ptile, xx, 'width', 10, 'support', [0,100.1]);
[pdf,xi] = ksdensity(ptile, 'support', [0,101]);
% bar(xi, pdf, 'EdgeColor', 'none');
xbox = [0; xi'; 100];
ybox = [0; pdf'; 0];
fill(xbox, ybox, [0,0,1], 'FaceAlpha', 1, 'LineStyle','None');
xlabel('Score percentile');
ylabel('Probability density');
axis square
axis([0, 100, 0, 0.151]);
 

set(gcf,'Color',[1,1,1]);


