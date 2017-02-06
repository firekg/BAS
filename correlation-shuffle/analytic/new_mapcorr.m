function rho=mapcorr(ixMu,ixNu,dist,maskprod)


[~,M]=size(ixMu);
% [R,M]=size(ixMu);

ixMu=ixMu(:);
ixNu=ixNu(:);

MuMu=dist(ixMu,ixMu);
NuNu=dist(ixNu,ixNu);
MuNu=dist(ixMu,ixNu);


[VMu,VNu]=deal(zeros(M,1));
for i=1:M
    VMu(i)=maskprod{i,i}*MuMu(:);
    VNu(i)=maskprod{i,i}*NuNu(:);
end

C=zeros(M,M);
for i=1:M
    for j=1:M        
        C(i,j)=maskprod{i,j}*MuNu(:);
    end
end

rho=C./sqrt(VMu*VNu');

