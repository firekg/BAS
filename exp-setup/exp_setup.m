%% setup: training phase
[Mo,~,~]=Rand_ImOrder_RevN_RevPos();
dlmwrite('TrImageLoadOrderT6.txt', Mo, ' ');

%% setup: active learning phase
% get image order
[Mo,Mn,~]=Rand_ImOrder_RevN_RevPos();
dlmwrite('ALImageLoadOrderT6.txt', Mo, ' ');
dlmwrite('ALRevealingNumberT6.txt', Mn, ' ');

%% setup: passive learning phase: random, BALD, and antiBALD revealing

% random revealing
[ImOrder_r,RevN_r,RevPos_r]=Rand_ImOrder_RevN_RevPos();
RevTyp_r=zeros(200,1); %0

% BALD patchy
[ImOrder_p, RevN_p, RevID_p]=Gen_Rev_Sequence(1,20,20);
nt=size(ImOrder_p,1);  RevPos_p=zeros(nt,30);
for i=1:nt
    RevPos_p(i,1:RevN_p(i)*2) = Get_RevPosM_cm(RevN_p(i),RevID_p(i,:));
end
RevTyp_p=ones(200,1); %1

% BALD stripy horizontal
[ImOrder_sh, RevN_sh, RevID_sh]=Gen_Rev_Sequence(2,6,30);
nt=size(ImOrder_sh,1);  RevPos_sh=zeros(nt,30);
for i=1:nt
    RevPos_sh(i,1:RevN_sh(i)*2) = Get_RevPosM_cm(RevN_sh(i),RevID_sh(i,:));
end
RevTyp_sh=ones(200,1)*2; %2

% BALD stripy vertical
[ImOrder_sv, RevN_sv, RevID_sv]=Gen_Rev_Sequence(2,30,6);
nt=size(ImOrder_sv,1);  RevPos_sv=zeros(nt,30);
for i=1:nt
    RevPos_sv(i,1:RevN_sv(i)*2) = Get_RevPosM_cm(RevN_sv(i),RevID_sv(i,:));
end
RevTyp_sv=ones(200,1)*3; %3

% combine them
nt=800;  rr=randperm(nt);
ImOrder = [ImOrder_r; ImOrder_p; ImOrder_sh; ImOrder_sv];
RevN    = [   RevN_r;    RevN_p;    RevN_sh;    RevN_sv];
RevPos  = [ RevPos_r;  RevPos_p;  RevPos_sh;  RevPos_sv];
RevTyp  = [ RevTyp_r;  RevTyp_p;  RevTyp_sh;  RevTyp_sv];
PL=Gen_Fake_PL_DATA(ImOrder,RevN,RevPos,RevTyp);

% mix them up
IO=zeros(nt,4);  RN=zeros(nt,1);  RP=zeros(nt,30);  RT=zeros(nt,1);
for i=1:nt
    IO(i,:) = ImOrder(rr(i),:);
    RN(i,:) = RevN(rr(i),:);
    RP(i,:) = RevPos(rr(i),:);
    RT(i,:) = RevTyp(rr(i),:);
end

% save them as text
dlmwrite('PLImageLoadOrderT2a.txt', IO(1:400,:), ' ');
dlmwrite('PLRevealingNumberT2a.txt', RN(1:400,:), ' ');
dlmwrite('PLRevealingPositionT2a.txt', RP(1:400,:), ' ');
dlmwrite('PLRevealingTypeT2a.txt', RT(1:400,:), ' ');

dlmwrite('PLImageLoadOrderT2b.txt', IO(401:800,:), ' ');
dlmwrite('PLRevealingNumberT2b.txt', RN(401:800,:), ' ');
dlmwrite('PLRevealingPositionT2b.txt', RP(401:800,:), ' ');
dlmwrite('PLRevealingTypeT2b.txt', RT(401:800,:), ' ');

%% test: does antiBALD revealing make sense?
[ImOrder_aB, RevN_aB, RevID_aB, RevTyp_aB]=Gen_TESTantiBALD_Sequence();
nt=size(ImOrder_aB,1);  RevPos_aB=zeros(nt,30);
for i=1:nt
    RevPos_aB(i,1:RevN_aB(i)*2) = Get_RevPosM_cm(RevN_aB(i),RevID_aB(i,:));
end
PL=Gen_Fake_PL_DATA(ImOrder_aB, RevN_aB, RevPos_aB, RevTyp_aB);

%% test: do the data make sense?
PLR=Get_PL_of_RevealType(PL,0);
[PLB,PLaB] = Get_PL_of_Bald_antiBald(PL);
AL=PLR;
PERF=Get_Performances(AL,PLR,PLB,PLaB);
%% test: plot performance
plot_performances(PERF)


