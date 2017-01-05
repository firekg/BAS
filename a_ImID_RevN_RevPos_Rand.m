function [Mi,Mn,Mp]=a_ImID_RevN_RevPos_Rand(nt,istart)
% nt = number of trials
% istart = image number start
% Mi = matrix image ID
% Mn = matrix revaling number
% Mp = matrix revaling position

ShowOrder=zeros(nt,4);  Mi=zeros(nt,4);
for i=1:nt/2    ShowOrder(i,:)=[1,20,20,i+istart];         end
for i=1:nt/4    ShowOrder(nt/2+i,:)=[2,6,30,i+istart];     end
for i=1:nt/4    ShowOrder(nt*3/4+i,:)=[2,30,6,i+istart];   end

% get revealing number & position
ni=5; nd=5; nf=25;  % number of revealings 3:3:15
nrev = ni:nd:nf;
nrevv=repmat(nrev,1,nt/length(nrev));
Mn=zeros(nt,1);

% mix
rr=randperm(nt);
for i=1:nt
    Mi(i,:)=ShowOrder(rr(i),:);
    Mn(i)=nrevv(rr(i));
end

% Mp=rand(nt,nf*2)*20-10;  % uniform between (-10,10)

Mp=randn(1,nt*nf*2)*(20/3);  % s.d. = image size/3
while (sum(abs(Mp)>10)>0)
    ind = abs(Mp)>10;
    ng = sum(ind);
    Mp(ind)=randn(1,ng);    
end
Mp=reshape(Mp,nt,nf*2);

end