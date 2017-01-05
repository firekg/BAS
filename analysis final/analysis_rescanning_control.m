% analysis rescanning control

%% load no rescanning control 2015-10-14
[AJAL1,~,~] = DATAFILE_Read('data\AJ20151016ALRCM5.DAT');
[AJAL2,~,~] = DATAFILE_Read('data\AJ20151016ALRCM6.DAT');
[AJAL3,~,~] = DATAFILE_Read('data\AJ20151016ALRCM7.DAT');
[AJAL4,~,~] = DATAFILE_Read('data\AJ20151017ALRCM8.DAT');
[AJAL5,~,~] = DATAFILE_Read('data\AJ20151017ALRCM9.DAT');
[AJAL6,~,~] = DATAFILE_Read('data\AJ20151017ALRCM10.DAT');
AJAL1.FrameData=[];
AJAL2.FrameData=[];
AJAL3.FrameData=[];
AJAL4.FrameData=[];
AJAL5.FrameData=[];
AJAL6.FrameData=[];
AJAL = a_Combine_Data(AJAL1,AJAL2,AJAL3,AJAL4,AJAL5,AJAL6);
% saved as AJdata.mat

%%
[ATAL1,~,~] = DATAFILE_Read('data\AT20151017ALRCM1.DAT');
[ATAL2,~,~] = DATAFILE_Read('data\AT20151017ALRCM2.DAT');
[ATAL3,~,~] = DATAFILE_Read('data\AT20151017ALRCM3.DAT');
[ATAL4,~,~] = DATAFILE_Read('data\AT20151019ALRCM4.DAT');
[ATAL5,~,~] = DATAFILE_Read('data\AT20151019ALRCM5.DAT');
[ATAL6,~,~] = DATAFILE_Read('data\AT20151019ALRCM6.DAT');
ATAL1.FrameData=[];
ATAL2.FrameData=[];
ATAL3.FrameData=[];
ATAL4.FrameData=[];
ATAL5.FrameData=[];
ATAL6.FrameData=[];
ATAL = a_Combine_Data(ATAL1,ATAL2,ATAL3,ATAL4,ATAL5,ATAL6);

%%
[KHAL1,~,~] = DATAFILE_Read('data\KH20151025ALRCM1.DAT');
[KHAL2,~,~] = DATAFILE_Read('data\KH20151025ALRCM2.DAT');
[KHAL3,~,~] = DATAFILE_Read('data\KH20151026ALRCM3.DAT');
[KHAL4,~,~] = DATAFILE_Read('data\KH20151026ALRCM4.DAT');
[KHAL5,~,~] = DATAFILE_Read('data\KH20151026ALRCM5.DAT');
[KHAL6,~,~] = DATAFILE_Read('data\KH20151026ALRCM6.DAT');
KHAL1.FrameData=[];
KHAL2.FrameData=[];
KHAL3.FrameData=[];
KHAL4.FrameData=[];
KHAL5.FrameData=[];
KHAL6.FrameData=[];
KHAL = a_Combine_Data(KHAL1,KHAL2,KHAL3,KHAL4,KHAL5,KHAL6);

%% perf

PERF = zeros(5,4,5,2);
mPERF = zeros(25,4,10,2);
INFO = zeros(25,4,10,2);

nparti = 4;
% for parti = 1:nparti
%     if parti == 1
%         load('AJdata.mat', 'AJAL');
%         D = AJAL;
%         %sym = 'b-o';
%     elseif parti == 2
%         load('ATdata.mat', 'ATAL');
%         D = ATAL;
%         %sym = 'b-s';
%     elseif parti == 3
%         load('KHdata.mat', 'KHAL');
%         D = KHAL;
%         %sym = 'b-^';
%     elseif parti == 4        
%         D = a_Combine_Data(AJAL, ATAL, KHAL);
%     end
%     [PERF(:,parti,1,1), PERF(:,parti,1,2), ~, ~] = a_Get_Human_Perf(D);
%     %correct = sum(D.AnswerChoice == D.AnswerReal)/D.Trials;
%     %fprintf('parti %d, fraction correct = %f\n', parti, correct);
%     %hold on;
%     %[PerfAL, sdAL, ~, x] = a_Get_Human_Perf(D);
%     %errorbar(x,PerfAL,sdAL,sym);
% end

for parti = 1:nparti
    if parti == 1
        load('SCdata.mat', 'SCAL');
        D = SCAL;
        %sym = 'b-o';
    elseif parti == 2
        load('EPdata.mat', 'EPAL');
        D = EPAL;
        %sym = 'b-s';
    elseif parti == 3
        load('BZdata.mat', 'BZAL');
        D = BZAL;
        %sym = 'b-^';
    elseif parti == 4        
        D = a_Combine_Data(SCAL, EPAL, BZAL);
    end
    [PERF(:,parti,1,1), PERF(:,parti,1,2), ~, ~] = a_Get_Human_Perf(D);
    %correct = sum(D.AnswerChoice == D.AnswerReal)/D.Trials;
    %fprintf('parti %d, fraction correct = %f\n', parti, correct);
    %hold on;
    %[PerfAL, sdAL, ~, x] = a_Get_Human_Perf(D);
    %errorbar(x,PerfAL,sdAL,sym);
end

% plot([0,30],[0.5,0.5],':','Color',[0.5,0.5,0.5]);
% plot([0,30],[1,1],':','Color',[0.5,0.5,0.5]);
% 
% ylabel('Fraction correct)');
% xlabel('Revealing number');
% axis([0 27.5 0.31 1.05]);
% legend('P1','P2','P3','P4','P5','P6');
% 
% set(gcf,'Color',[1,1,1]);
% % save figure
% % export_fig('Perf.pdf')

%% density maps
Drp = DRemovePhase(D);
[REV, RrevMap, RdrevMap] = revmap_sim(D);
plot_rev_maps(RrevMap, RdrevMap);

%% correlation
% clear;
% load('ATdata.mat');
% D = ATAL;

[REV, ~, ~] = revmap_sim(Drp);
subj = 1;
analysis_hinton;

CORR = zeros(25,3,2,2,4);
parti = 1;
cond = 1;
muw = mean(within,2);
sdw = std(within,[],2);
CORR(:,1,1,cond,parti) = muw;
CORR(:,2,1,cond,parti) = muw+sdw;
CORR(:,3,1,cond,parti) = muw-sdw;
muw = mean(between,2);
sdw = std(between,[],2);
CORR(:,1,2,cond,parti) = muw;
CORR(:,2,2,cond,parti) = muw+sdw;
CORR(:,3,2,cond,parti) = muw-sdw;

CORR25 = zeros(3,2,2,4);
plot_corr(CORR, CORR25);

%% correlation 25
% clear;
% load('AJdata.mat');
% D = AJAL;

[REV, ~, ~] = revmap_sim(D);
subj = 1;
analysis_hinton_rev25;

CORR25 = zeros(3,2,2,4);
parti=1;
cond=1;
within = sort(within);
n = numel(within);
CORR25(1,1,cond,parti) = mean(within);
CORR25(2,1,cond,parti) = within(round(n*0.95));
CORR25(3,1,cond,parti) = within(round(n*0.05));
between = sort(between);
n = numel(between);
CORR25(1,2,cond,parti) = mean(between);
CORR25(2,2,cond,parti) = between(round(n*0.95));
CORR25(3,2,cond,parti) = between(round(n*0.05));
disp(sprintf('Parti%d-cond%d same mean: %f.', parti, cond, mean(within)));
disp(sprintf('Parti%d-cond%d diff mean: %f.', parti, cond, mean(between)));
test = within<0;
temp = sum(test);
disp(sprintf('Parti%d-cond%d same: %f.', parti, cond, temp/1000));
test = between>0;
temp = sum(test);
disp(sprintf('Parti%d-cond%d; diff: %f.', parti, cond, temp/1000));
test = within>between;
temp = sum(test);
disp(sprintf('Parti%d-cond%d same>diff: %f.', parti, cond, temp/1000));


