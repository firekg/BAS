function DIM=a_DIM(D)

for trial=1:D.Trials
    
    Type = D.ImageID(trial,1);
    X    = D.ImageID(trial,2);    
    Y    = D.ImageID(trial,3);
    Num  = D.ImageID(trial,4);
    DIM(trial) = load(strcat('im/Type',num2str(Type),'X',num2str(X),'Y',num2str(Y),'Num', num2str(Num),'raw.mat'));
    
end

end

