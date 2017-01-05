function GaussFilter=a_GaussFilter(mapsize)
% mapsize: number of rows of map to be filtered
npix=770;  bin=npix/mapsize;
var=(20/bin)^2; %variance is fixed at 20 for npix=770
% var=(40/bin)^2; % bigger variance - 2015-10-17
n=ceil(sqrt(var)*6);
[x1,x2]=meshgrid(1:n,1:n);
z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % in pixel coordinate
GaussFilter = zeros(n*n,1);
for i=1:n*n
    GaussFilter(i) = mvnpdf(z(i,:),[1,1]*n/2,[1,0;0,1]*var);
end
GaussFilter = reshape(GaussFilter,n,n);

end