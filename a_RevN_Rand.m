function [Mn]=a_RevN_Rand(nt,ni,nd,nf)

nrev = ni:nd:nf;
nrevv=repmat(nrev,1,nt/length(nrev));
rr=randperm(nt);
Mn=zeros(nt,1);
for i=1:nt
    Mn(i)=nrevv(rr(i));
end

end