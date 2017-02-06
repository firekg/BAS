function rho=mapcorr(ixMu,ixNu,dist,maskprod)


[~,M]=size(ixMu);
% [R,M]=size(ixMu);

ixMu=ixMu(:);
ixNu=ixNu(:);

MuMu=dist(ixMu,ixMu);
NuNu=dist(ixNu,ixNu);
MuNu=dist(ixMu,ixNu);

% vv=repmat(-1/(M*R),[M*R M]);
% for i=1:M
%     vv((i-1)*R+1:i*R,i)=((M-1)/(M*R));
% end

[VMu,VNu]=deal(zeros(M,1));
for i=1:M
    vvv=maskprod(:,:,i,i);
%     vvv=vv(:,i)*vv(:,i)';
    VMu(i)=sum(sum(vvv.*MuMu));
    VNu(i)=sum(sum(vvv.*NuNu));
end

C=zeros(M,M);
for i=1:M
    for j=1:M        
        C(i,j)=sum(sum(maskprod(:,:,i,j).*MuNu));
%         C(i,j)=sum(sum((vv(:,i)*vv(:,j)').*MuNu));
    end
end

rho=C./sqrt(VMu*VNu');
