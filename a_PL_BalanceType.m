function Dt=a_PL_BalanceType(D)
% input: PLR or PLB or PLaB

fields = fieldnames(D);
for i = 1:numel(fields)
    Dt.(fields{i})=[];
end    
Dt.Trials=0;

for nr=5:5:25

    D1 = a_Get_RevealNumber(D,nr);
    Dt1 = a_Get_AnswerReal(D1,1);
    Dt2 = a_Get_AnswerReal(D1,2);
    if(Dt1.Trials > Dt2.Trials)       
        dt = Dt1.Trials - Dt2.Trials; 
        for i = 2:numel(fields)  %start from 2 to skip Dt.Trials
            temp = Dt1.(fields{i}); 
            temp(1:2:dt*2-1,:)=[];
            Dt1.(fields{i})=temp;
        end
        Dt1.Trials = Dt2.Trials;
    elseif(Dt2.Trials > Dt1.Trials)
        dt = Dt2.Trials - Dt1.Trials; 
        for i = 2:numel(fields)   %start from 2 to skip Dt.Trials
            temp = Dt2.(fields{i}); 
            temp(1:2:dt*2-1,:)=[];
            Dt2.(fields{i})=temp;
         end
        Dt2.Trials = Dt1.Trials;
    end
    Dtc = a_Combine_Data(Dt1,Dt2);
    Dt = a_Combine_Data(Dt,Dtc);

end
    
end
