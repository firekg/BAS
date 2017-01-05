function isoperiR=a_isoperiR(D)
%%
ntrial = D.Trials;
%ntrial=800;
nrmax  = 25;
isoperiR  = NaN(ntrial,nrmax);

for trial=1:ntrial;        
    nr=D.MaxRevealingTrial(trial);
    %nr=25;
    dis=zeros(nr);
    x=[];  y=[];
    x=D.RevealPosX(trial,1:nr);
    y=D.RevealPosY(trial,1:nr);

  % get distance matrix (dis) and angle matrix (ang)
    for i=1:nr
    for j=1:nr
        dx=x(i)-x(j);  dy=y(i)-y(j);
        dis(i,j)= sqrt(dx^2+dy^2);
    end
    end

    for i=3:nr
        isoperiR(trial,i)  = a_get_isoperiR(x(1:i),y(1:i),dis(1:i,1:i));
    end
        
end

end

%% isoperimetric ratio of convex hull normalized by circle
function isoperiR = a_get_isoperiR(x,y,dis)

    [hull, hullarea] = convhull(x,y);
    peri=0;
    for k=1:length(hull)-1
        peri=peri+dis(hull(k),hull(k+1));
    end
    isoperiR = hullarea/peri^2*4*pi;

end