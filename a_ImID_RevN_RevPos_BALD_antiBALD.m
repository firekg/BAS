%%
% 800 trials in total: 200 random; 200 BALD; 400 anti-BLAD
%  
% 200 random    --> 100pa + 50sh + 50sv 
%               --> Image ID: 301:400pa, 301:350sh, 301:350sv
%   
% 200 BALD      --> 100pa + 50sh + 50sv
%               --> Image ID: 1:100pa, 1:50sh, 1:50sv 
%                  (Image ID = Reveal ID)
% 
% 400 anti-BALD --> Reveal ID / Image ID:
%             (100) 101:200pa / 51:100sh,  51:100sv
%             (150) 101:250sh / 101:200pa, 351:400sv 
%             (150) 101:250sv / 201:300pa, 351:400sh 

% the BALD & antiBALD is equivalent to:
% 
% 200pa --> Reveal ID / Image ID
%           1:200pa   / 1:100pa, 51:100sh, 51:100sv
% 
% 200sh --> Reveal ID / Image ID
%           1:50sh, 101:250sh / 1:50sh, 101:200pa, 351:400sv
% 
% 200sv --> Reveal ID / Image ID
%           1:50sv, 101:250sv / 1:50sv, 201:300pa, 351:400sh

% 200 random
[ImID_r,RevN_r,RevPos_r]=a_ImID_RevN_RevPos_Rand(200,300);
RevTyp_r=zeros(200,1); %0

load('im/BALDrev_w3_noise05_n110_nr25_full.mat')

%  200 pa
ImID_pa = [ ones(100,1), ones(100,1)*20, ones(100,1)*20, (1:100)';...
            ones(50,1)*2, ones(50,1)*6, ones(50,1)*30, (51:100)';...
            ones(50,1)*2, ones(50,1)*30, ones(50,1)*6, (51:100)'];
RevN_pa = repmat((5:5:25)',40,1);
RevID_pa = [ ones(200,1), ones(200,1)*20, ones(200,1)*20, (1:200)'];
RevPos_pa = zeros(200,50);
RevPos_pa(:,1:2:49) = BALDrevTpa.W_BALDx(:,1:200)';
RevPos_pa(:,2:2:50) = BALDrevTpa.W_BALDy(:,1:200)';
RevTyp_pa=ones(200,1); %1

%  200 sh
ImID_sh = [ ones(50,1)*2, ones(50,1)*6, ones(50,1)*30, (1:50)';...
            ones(100,1), ones(100,1)*20, ones(100,1)*20, (101:200)';...
            ones(50,1)*2, ones(50,1)*30, ones(50,1)*6, (351:400)'];
RevN_sh = repmat((5:5:25)',40,1);
RevID_sh = [ ones(50,1)*2, ones(50,1)*6, ones(50,1)*30, (1:50)';...
             ones(150,1)*2, ones(150,1)*6, ones(150,1)*30, (101:250)'];
RevPos_sh = zeros(200,50);
RevPos_sh(1:50,1:2:49) = BALDrevTsh.W_BALDx(:,1:50)';
RevPos_sh(1:50,2:2:50) = BALDrevTsh.W_BALDy(:,1:50)';
RevPos_sh(51:200,1:2:49) = BALDrevTsh.W_BALDx(:,101:250)';
RevPos_sh(51:200,2:2:50) = BALDrevTsh.W_BALDy(:,101:250)';
RevTyp_sh=ones(200,1)*2; %2

%  200 sv
ImID_sv = [ ones(50,1)*2, ones(50,1)*30, ones(50,1)*6, (1:50)';...
            ones(100,1), ones(100,1)*20, ones(100,1)*20, (201:300)';...
            ones(50,1)*2, ones(50,1)*6, ones(50,1)*30, (351:400)'];
RevN_sv = repmat((5:5:25)',40,1);
RevID_sv = [ ones(50,1)*2, ones(50,1)*30, ones(50,1)*6, (1:50)';...
             ones(150,1)*2, ones(150,1)*30, ones(150,1)*6, (101:250)'];
RevPos_sv = zeros(200,50);
RevPos_sv(1:50,1:2:49) = BALDrevTsv.W_BALDx(:,1:50)';
RevPos_sv(1:50,2:2:50) = BALDrevTsv.W_BALDy(:,1:50)';
RevPos_sv(51:200,1:2:49) = BALDrevTsv.W_BALDx(:,101:250)';
RevPos_sv(51:200,2:2:50) = BALDrevTsv.W_BALDy(:,101:250)';
RevTyp_sv=ones(400,1)*3; %3

% combine
nt=800;  rr=randperm(nt);
ImID   = [   ImID_r;    ImID_pa;    ImID_sh;    ImID_sv];
RevN   = [   RevN_r;    RevN_pa;    RevN_sh;    RevN_sv];
RevPos = [ RevPos_r;  RevPos_pa;  RevPos_sh;  RevPos_sv];
RevTyp = [ RevTyp_r;  RevTyp_pa;  RevTyp_sh;  RevTyp_sv];

% trim
RevPos(RevPos>9.9)=9.9;
RevPos(RevPos<-9.9)=-9.9;

% mix
II=zeros(nt,4);  RN=zeros(nt,1);  RP=zeros(nt,50);  RT=zeros(nt,1);
for i=1:nt
    II(i,:) = ImID(rr(i),:);
    RN(i,:) = RevN(rr(i),:);
    RP(i,:) = RevPos(rr(i),:);
    RT(i,:) = RevTyp(rr(i),:);
end

% save as 8 medium experiments, each 100 trials
for i=1:8
    dlmwrite(strcat('cfile/PLImageLoadOrderM',num2str(i),'.txt'), II((i-1)*100+1:i*100,:), ' ');
    dlmwrite(strcat('cfile/PLRevealingNumberM',num2str(i),'.txt'), RN((i-1)*100+1:i*100,:), ' ');
    dlmwrite(strcat('cfile/PLRevealingPositionM',num2str(i),'.txt'), RP((i-1)*100+1:i*100,:), ' ');
    dlmwrite(strcat('cfile/PLRevealingTypeM',num2str(i),'.txt'), RT((i-1)*100+1:i*100,:), ' ');
end


%%
% % 400 random
% [ImID_r,RevN_r,RevPos_r]=a_ImID_RevN_RevPos_Rand(400,600);
% RevTyp_r=zeros(400,1); %0
% 
% % the BALD & antiBALD is equivalent to:
% % 
% % 400pa --> Reveal ID / Image ID
% %           1:400pa   / 1:200pa, 101:200sh, 101:200sv
% % 
% % 400sh --> Reveal ID / Image ID
% %           1:100sh, 201:500sh / 1:100sh, 201:400pa, 701:800sv
% % 
% % 400sv --> Reveal ID / Image ID
% %           1:100sv, 201:500sv / 1:100sv, 401:600pa, 701:800sh
% 
% % 1600 trials in total: 400 random; 400 BALD; 800 anti-BLAD
% %  
% % 400 random    --> 200pa + 100sh + 100sv 
% %               --> Image ID: 601:800pa, 601:700sh, 601:700sv
% %   
% % 400 BALD      --> 200pa + 100sh + 100sv
% %               --> Image ID: 1:200pa, 1:100sh, 1:100sv 
% %                  (Image ID = Reveal ID)
% % 
% % 800 anti-BALD --> Reveal ID / Image ID:
% %             (200) 201:400pa / 101:200sh, 101:200sv
% %             (300) 201:500sh / 201:400pa, 701:800sv 
% %             (300) 201:500sv / 401:600pa, 701:800sh 
% 
% load('im/BALDrev_w3_noise05_n110_nr25_full.mat')
% 
% %  400 pa
% ImID_pa = [ ones(200,1), ones(200,1)*20, ones(200,1)*20, (1:200)';...
%             ones(100,1)*2, ones(100,1)*6, ones(100,1)*30, (101:200)';...
%             ones(100,1)*2, ones(100,1)*30, ones(100,1)*6, (101:200)'];
% RevN_pa = repmat((5:5:25)',80,1);
% RevID_pa = [ ones(400,1), ones(400,1)*20, ones(400,1)*20, (1:400)'];
% RevPos_pa = zeros(400,50);
% RevPos_pa(:,1:2:49) = BALDrevTpa.W_BALDx(:,1:400)';
% RevPos_pa(:,2:2:50) = BALDrevTpa.W_BALDy(:,1:400)';
% RevTyp_pa=ones(400,1); %1
% 
% %  400 sh
% ImID_sh = [ ones(100,1)*2, ones(100,1)*6, ones(100,1)*30, (1:100)';...
%             ones(200,1), ones(200,1)*20, ones(200,1)*20, (201:400)';...
%             ones(100,1)*2, ones(100,1)*30, ones(100,1)*6, (701:800)'];
% RevN_sh = repmat((5:5:25)',80,1);
% RevID_sh = [ ones(100,1)*2, ones(100,1)*6, ones(100,1)*30, (1:100)';...
%              ones(300,1)*2, ones(300,1)*6, ones(300,1)*30, (201:500)'];
% RevPos_sh = zeros(400,50);
% RevPos_sh(1:100,1:2:49) = BALDrevTsh.W_BALDx(:,1:100)';
% RevPos_sh(1:100,2:2:50) = BALDrevTsh.W_BALDy(:,1:100)';
% RevPos_sh(101:400,1:2:49) = BALDrevTsh.W_BALDx(:,201:500)';
% RevPos_sh(101:400,2:2:50) = BALDrevTsh.W_BALDy(:,201:500)';
% RevTyp_sh=ones(400,1)*2; %2
% 
% %  400 sv
% ImID_sv = [ ones(100,1)*2, ones(100,1)*30, ones(100,1)*6, (1:100)';...
%             ones(200,1), ones(200,1)*20, ones(200,1)*20, (401:600)';...
%             ones(100,1)*2, ones(100,1)*6, ones(100,1)*30, (701:800)'];
% RevN_sv = repmat((5:5:25)',80,1);
% RevID_sv = [ ones(100,1)*2, ones(100,1)*30, ones(100,1)*6, (1:100)';...
%              ones(300,1)*2, ones(300,1)*30, ones(300,1)*6, (201:500)'];
% RevPos_sv = zeros(400,50);
% RevPos_sv(1:100,1:2:49) = BALDrevTsv.W_BALDx(:,1:100)';
% RevPos_sv(1:100,2:2:50) = BALDrevTsv.W_BALDy(:,1:100)';
% RevPos_sv(101:400,1:2:49) = BALDrevTsv.W_BALDx(:,201:500)';
% RevPos_sv(101:400,2:2:50) = BALDrevTsv.W_BALDy(:,201:500)';
% RevTyp_sv=ones(400,1)*3; %3
% 
% % combine
% nt=1600;  rr=randperm(nt);
% ImID   = [   ImID_r;    ImID_pa;    ImID_sh;    ImID_sv];
% RevN   = [   RevN_r;    RevN_pa;    RevN_sh;    RevN_sv];
% RevPos = [ RevPos_r;  RevPos_pa;  RevPos_sh;  RevPos_sv];
% RevTyp = [ RevTyp_r;  RevTyp_pa;  RevTyp_sh;  RevTyp_sv];
% 
% % trim
% RevPos(RevPos>9.9)=9.9;
% RevPos(RevPos<-9.9)=-9.9;
% 
% % mix
% II=zeros(nt,4);  RN=zeros(nt,1);  RP=zeros(nt,50);  RT=zeros(nt,1);
% for i=1:nt
%     II(i,:) = ImID(rr(i),:);
%     RN(i,:) = RevN(rr(i),:);
%     RP(i,:) = RevPos(rr(i),:);
%     RT(i,:) = RevTyp(rr(i),:);
% end
% 
% % save as 16 medium experiments, each 100 trials
% for i=1:16
%     dlmwrite(strcat('cfile/PLImageLoadOrderM',num2str(i),'.txt'), II((i-1)*100+1:i*100,:), ' ');
%     dlmwrite(strcat('cfile/PLRevealingNumberM',num2str(i),'.txt'), RN((i-1)*100+1:i*100,:), ' ');
%     dlmwrite(strcat('cfile/PLRevealingPositionM',num2str(i),'.txt'), RP((i-1)*100+1:i*100,:), ' ');
%     dlmwrite(strcat('cfile/PLRevealingTypeM',num2str(i),'.txt'), RT((i-1)*100+1:i*100,:), ' ');
% end


