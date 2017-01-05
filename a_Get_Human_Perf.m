function [perf, sd, time, x]=a_Get_Human_Perf(D)

RevN=D.MaxRevealingTrial;
P=(D.AnswerChoice==D.AnswerReal);
t=D.StateAnswerTime/1000;
uni_revn = unique(RevN);
rn =length(uni_revn);

perf=zeros(rn,1);
sd=zeros(rn,1);
time=zeros(rn,1);

for i=1:rn
   correct =sum(P(RevN==uni_revn(i)));
   total = length(P(RevN==uni_revn(i)));
   perf(i) = correct/total;
   sd(i)=sqrt(perf(i)*(1-perf(i))/total);
   time(i) = mean(t(RevN==uni_revn(i)));
end
x=uni_revn;