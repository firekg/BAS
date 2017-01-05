%% Image & revealing prep:

% Generates GP image, stored in folder 'Exp_a\im\'
a_GenGPimage(1,20,20,1,800);
a_GenGPimage(2,6,30,1,800);
a_GenGPimage(2,30,6,1,800);

%
% Generate BALD revealings from images.  To be used in BALD & antiBALD
BALDrevTpa=a_Get_BALD(1,20,20,1,800); %BALD revealing Type patchy
BALDrevTsh=a_Get_BALD(2,6,30,1,800);  %BALD revealing Type stripy horizontal
BALDrevTsv=a_Get_BALD(2,30,6,1,800);  %BALD revealing Type stripy vertical
% Save in 'im/BALDrev_w7_noise05_n70.mat'
% Save in 'im/BALDrev_w3_noise05_n110.mat'
% Save in 'im/BALDrev_w1_noise01_n155.mat'
% Save in 'im/BALDrev_w1_noise01_n155_nr25.mat'
% Save in 'im/BALDrev_w3_noise05_n110_nr25.mat'
% Save in 'im/BALDrev_w7_noise05_n70_nr25.mat'
% Save in 'im/BALDrev_w3_noise05_n110_nr25_full.mat'
save('im/BALDrev_w3_noise05_n110_nr25_full.mat')

%% AL short (40 trials) x 40
% Produce 40 short experiments, each 40 trials, always 25 revealings.
% No image repeated.
% Use up all 800 pa, 400 sh, 400 sv
ntrial=40;
for ex=1:40
    % training
    [Mi,~,~]=a_ImID_RevN_RevPos_Rand(ntrial,(ex-1)*ntrial/2);
    dlmwrite(strcat('cfile/TrImageLoadOrderS',num2str(ex),'.txt'), Mi, ' ');
    % active learning
    [Mi,~,~]=a_ImID_RevN_RevPos_Rand(ntrial,(ex-1)*ntrial/2);
    Mn=25*ones(ntrial,1);
    dlmwrite(strcat('cfile/ALImageLoadOrderS',num2str(ex),'.txt'), Mi, ' ');
    dlmwrite(strcat('cfile/ALRevealingNumberS',num2str(ex),'.txt'), Mn, ' ');
end

%% AL medium (100 trials) x 16
% Produce 16 medium experiment, each 100 trials, variable revealings.
% Uses 1:800 pa
ntrial=100;
for ex=1:16
    % training
    [Mi,~,~]=a_ImID_RevN_RevPos_Rand(ntrial,(ex-1)*ntrial/2);
    dlmwrite(strcat('cfile/TrImageLoadOrderM',num2str(ex),'.txt'), Mi, ' ');
    % active learning
    [Mi,~,~]=a_ImID_RevN_RevPos_Rand(ntrial,(ex-1)*ntrial/2);
    [Mn]=a_RevN_Rand(ntrial,5,5,25);
    dlmwrite(strcat('cfile/ALImageLoadOrderM',num2str(ex),'.txt'), Mi, ' ');
    dlmwrite(strcat('cfile/ALRevealingNumberM',num2str(ex),'.txt'), Mn, ' ');
end

% All outputs stored in folder 'Exp_a\cfile\'

%% PL medium (100 trials) x 16
% Produce 16 medium experiments, each 100 trials, variable revealings.
% Mix of BALD, antiBALD, and random.
a_ImID_RevN_RevPos_BALD_antiBALD;

% All outputs stored in folder 'Exp_a\cfile\'


