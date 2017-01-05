function img=revMapColor(M)

rf=M;
rf(M<=0) = (1+M(M<=0)).^2.5;  %goes to 0
rf(M>0) = 1; %stays at 1

gf=M;
gf(M<=0) = (1 + M(M<=0)).^2.0; %goes to 0
gf(M>0) = (1 - M(M>0)).^2.0;  %goes to 0

bf=M;
bf(M<=0) = 1; %stays at 1
bf(M>0) = (1-M(M>0)).^2.5; %goes to 0

img(:,:,1)=rf;
img(:,:,2)=gf;
img(:,:,3)=bf;

end