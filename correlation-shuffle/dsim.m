function L=dsim(nx,ny,alpha)

% nbarx=sum(nx);
% nbary=sum(ny);

D=@(nz,beta) gammaln(sum(nz)+1)+gammaln(sum(zeros(size(nz))+beta))-gammaln(sum(nz+beta))+ ...
    + sum(gammaln(nz+beta)-gammaln(beta)-gammaln(nz+1));

L=D(nx+ny,alpha+nx+ny)-D(nx,alpha+nx)-D(ny,alpha+ny);




