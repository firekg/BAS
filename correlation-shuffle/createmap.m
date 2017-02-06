function M=createmap(x,y,gridx,gridy)

M=zeros(length(gridx),length(gridy));

gridx=[-inf gridx(1:end-1)+diff(gridx)/2 inf];
gridy=[-inf gridy(1:end-1)+diff(gridy)/2 inf];

[~,~,ybin]=histcounts(y,gridy);

for i=1:length(gridy)-1
   M(:,i)=histcounts(x(ybin==i),gridx); 
end

M=M/length(x);
