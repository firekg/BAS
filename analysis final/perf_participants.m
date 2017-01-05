function [PERF, Pmix, Bias, BiasMix, mPERF, mPmix, lpxPERF, lpxPmix, INFO, INFOC] = perf_participants(PROP)
%% can be used with input PROP from prop_sim_parfor
%  or can be used without input to load PROP from prop_sim
% CONTAINS mPerf_comb_prop

% PERF has 4 layers: revealing number x5(1:5:25); parti x4;
% condiitons x5 (AL,PLR,PLB,PLaB,PLBaBrn); moments x2 (mean,SD)
PERF = zeros(5,4,5,2);
% Pmix has 4 layers: revealing number x5(1:5:25); parti x4; 
% condiitons x9 (PA, SV, SH matches and mismatches); moments x2 (mean,SD)
Pmix = zeros(5,4,9,2);
% Bias has 4 layers: rev number; parti x4; condition x5; p-or-s x2
Bias = zeros(5,4,5,2);
% BiasMix has 4 layers: rev number; parti x4; condition x9; p-or-s x2
BiasMix = zeros(5,4,9,2);

% load exp data
load('SCdata.mat', 'SCAL', 'SCPLRbt', 'SCPLBbt', 'SCPLaBbt');
load('EPdata.mat', 'EPAL', 'EPPLRbt', 'EPPLBbt', 'EPPLaBbt');
load('BZdata.mat', 'BZAL', 'BZPLRbt', 'BZPLBbt', 'BZPLaBbt'); 

% loading expPROP
nprop = 5; %20;
nsimu = 10; %10;
% mPERFC and INFOC 5 layers: rev num; parti; 
% condition (AL,PLR, PLB, PLaB, nBAS, pBAS, heu1, heu2, heu3, AL-noPriorBias); moments; simu num
mPERFC = zeros(25,4,10,2, nprop*nsimu);
INFOC = zeros(25,4,10,2, nprop*nsimu);
% mPmixC 5 layers: rev number; parti; condiitons; moments; simu num 
mPmixC = zeros(25,4,9,2, nprop*nsimu);

% load lpx
load('lpx.mat');
lpxPERF  = zeros(25,4,5,2); % condition 1 should be skipped
lpxPmix = zeros(25,4,9,2);

aBasId = 4; % anti-BAS ID
nBasId = 8; % noisy BAS ID
%pBasId = 9; % prob-matching BAS ID

for parti=1:4
    % for PERF
    if(parti==1)
        AL=SCAL;  PLR=SCPLRbt;  PLB=SCPLBbt;  PLaB=SCPLaBbt;
    elseif(parti==2)
        AL=EPAL;  PLR=EPPLRbt;  PLB=EPPLBbt;  PLaB=EPPLaBbt;
    elseif(parti==3)
        AL=BZAL;  PLR=BZPLRbt;  PLB=BZPLBbt;  PLaB=BZPLaBbt;
    elseif(parti==4)  
        AL=a_Combine_Data(SCAL,EPAL,BZAL);
        PLR=a_Combine_Data(SCPLRbt,EPPLRbt,BZPLRbt);
        PLB=a_Combine_Data(SCPLBbt,EPPLBbt,BZPLBbt);
        PLaB=a_Combine_Data(SCPLaBbt,EPPLaBbt,BZPLaBbt);
    end    
    [PERF(:,parti,1,1), PERF(:,parti,1,2), ~, ~] = a_Get_Human_Perf(AL);
    [PERF(:,parti,2,1), PERF(:,parti,2,2), ~, ~] = a_Get_Human_Perf(PLR);
    [PERF(:,parti,3,1), PERF(:,parti,3,2), ~, ~] = a_Get_Human_Perf(PLB);
    [PERF(:,parti,4,1), PERF(:,parti,4,2), ~, ~] = a_Get_Human_Perf(PLaB);
    [Bias(:,parti,1,1), Bias(:,parti,1,2)] = getHumanChoiceBias(AL);
    [Bias(:,parti,2,1), Bias(:,parti,2,2)] = getHumanChoiceBias(PLR);
    [Bias(:,parti,3,1), Bias(:,parti,3,2)] = getHumanChoiceBias(PLB);
    [Bias(:,parti,4,1), Bias(:,parti,4,2)] = getHumanChoiceBias(PLaB);
    
    % for lpx
    if(parti<4)
        lpxPLR = lpx(parti,1);
        lpxPLB = lpx(parti,2);
        lpxPLaB = lpx(parti,3);
    elseif(parti==4)  
        lpxPLR = a_Combine_Data(lpx(1,1), lpx(2,1), lpx(3,1));
        lpxPLB = a_Combine_Data(lpx(1,2), lpx(2,2), lpx(3,2));
        lpxPLaB = a_Combine_Data(lpx(1,3), lpx(2,3), lpx(3,3));
    end    
    [lpxPERF(:,parti,2,1), lpxPERF(:,parti,2,2)] = modelPerfProp(lpxPLR.mla, PLR.AnswerReal, 'hard');
    [lpxPERF(:,parti,3,1), lpxPERF(:,parti,3,2)] = modelPerfProp(lpxPLB.mla, PLB.AnswerReal, 'hard');
    [lpxPERF(:,parti,4,1), lpxPERF(:,parti,4,2)] = modelPerfProp(lpxPLaB.mla, PLaB.AnswerReal, 'hard');
    
    % for Pmix
    PLBaB = a_Combine_Data(PLB, PLaB);
    for isimu = 1:nsimu;
        % load PROP
        %load(['expPROP_data/expPROP',num2str(isimu),'.mat']);
        %load(['expPROP_data/expPROPLowNoise',num2str(isimu),'.mat']);
        %load(['expPROP_data/expPROPLowNoiseNewPrior',num2str(isimu),'.mat']);
        %load(['expPROP_data/expPROPRevNoise',num2str(isimu),'.mat']);
        %load(['expPROP_data/expPROPRevNoiseNonRandObs.mat']); % only one
        %load(['expPROP_data/expPROPBiasPA0d2',num2str(isimu),'.mat']);
        %load(['expPROP_data/expPROPRefitNoise',num2str(isimu),'.mat']);
        %PROP = addFieldTrials(PROP);
        % used with input PROP generated with parfor
        prop = PROP(isimu).prop;
        prop = addFieldTrials(prop);
        for iprop = 1:nprop;
            isp =(isimu-1)*nprop+iprop;
            if (parti < 4)
                %mAL = PROP(parti,1,iprop);  
                %mPLR = PROP(parti,2,iprop);  
                %mPLB = PROP(parti,3,iprop);  
                %mPLaB = PROP(parti,13,iprop);
                % used with input PROP generated with parfor
                mAL = prop(parti,1,iprop);
                mPLR = prop(parti,2,iprop);
                mPLB = prop(parti,3,iprop);
                mPLaB = prop(parti,aBasId,iprop);
%                 mBAS = prop(parti,nBasId,iprop);
%                 %pBAS = prop(parti,pBasId,iprop);
%                 mHeu1 = prop(parti,5,iprop);
%                 mHeu2 = prop(parti,6,iprop);
%                 mHeu3 = prop(parti,7,iprop);
            else
                %mAL = a_Combine_Data(PROP(1,1,iprop), PROP(2,1,iprop), PROP(3,1,iprop));
                %mPLR = a_Combine_Data(PROP(1,2,iprop), PROP(2,2,iprop), PROP(3,2,iprop));
                %mPLB  =a_Combine_Data(PROP(1,3,iprop), PROP(2,3,iprop), PROP(3,3,iprop));
                %mPLaB = a_Combine_Data(PROP(1,13,iprop), PROP(2,13,iprop), PROP(3,13,iprop));
                % used with input PROP generated with parfor
                mAL = a_Combine_Data(prop(1,1,iprop), prop(2,1,iprop), prop(3,1,iprop));
                mPLR = a_Combine_Data(prop(1,2,iprop), prop(2,2,iprop), prop(3,2,iprop));
                mPLB  =a_Combine_Data(prop(1,3,iprop), prop(2,3,iprop), prop(3,3,iprop));
                mPLaB = a_Combine_Data(prop(1,aBasId,iprop), prop(2,aBasId,iprop), prop(3,aBasId,iprop));
%                 mBAS = a_Combine_Data(prop(1,nBasId,iprop), prop(2,nBasId,iprop), prop(3,nBasId,iprop));
%                 %pBAS = a_Combine_Data(prop(1,pBasId,iprop), prop(2,pBasId,iprop), prop(3,pBasId,iprop));
%                 mHeu1 = a_Combine_Data(prop(1,5,iprop), prop(2,5,iprop), prop(3,5,iprop));
%                 mHeu2 = a_Combine_Data(prop(1,6,iprop), prop(2,6,iprop), prop(3,6,iprop));
%                 mHeu3 = a_Combine_Data(prop(1,7,iprop), prop(2,7,iprop), prop(3,7,iprop));
            end
            
%             [mPERFC(:,parti,1,1,isp), mPERFC(:,parti,1,2,isp)] = modelPerfProp(mAL.softmla, AL.AnswerReal, 'soft');
%             [mPERFC(:,parti,2,1,isp), mPERFC(:,parti,2,2,isp)] = modelPerfProp(mPLR.softmla, PLR.AnswerReal, 'soft');
%             [mPERFC(:,parti,3,1,isp), mPERFC(:,parti,3,2,isp)] = modelPerfProp(mPLB.softmla, PLB.AnswerReal, 'soft');
%             [mPERFC(:,parti,4,1,isp), mPERFC(:,parti,4,2,isp)] = modelPerfProp(mPLaB.softmla, PLaB.AnswerReal, 'soft');
%             [mPERFC(:,parti,5,1,isp), mPERFC(:,parti,5,2,isp)] = modelPerfProp(mBAS.softmla, AL.AnswerReal, 'soft');
%             %[mPERFC(:,parti,6,1,isp), mPERFC(:,parti,6,2,isp)] = modelPerfProp(pBAS.softmla, AL.AnswerReal, 'soft');
%             [mPERFC(:,parti,7,1,isp), mPERFC(:,parti,7,2,isp)] = modelPerfProp(mHeu1.softmla, AL.AnswerReal, 'soft');
%             [mPERFC(:,parti,8,1,isp), mPERFC(:,parti,8,2,isp)] = modelPerfProp(mHeu2.softmla, AL.AnswerReal, 'soft');
%             [mPERFC(:,parti,9,1,isp), mPERFC(:,parti,9,2,isp)] = modelPerfProp(mHeu3.softmla, AL.AnswerReal, 'soft');
            
            [INFOC(:,parti,1,1,isp), INFOC(:,parti,1,2,isp)] = calInfo(mAL.mla);
            [INFOC(:,parti,2,1,isp), INFOC(:,parti,2,2,isp)] = calInfo(mPLR.mla);
            [INFOC(:,parti,3,1,isp), INFOC(:,parti,3,2,isp)] = calInfo(mPLB.mla);
            [INFOC(:,parti,4,1,isp), INFOC(:,parti,4,2,isp)] = calInfo(mPLaB.mla);
%             [INFOC(:,parti,5,1,isp), INFOC(:,parti,5,2,isp)] = calInfo(mBAS.mla);
%             %[INFOC(:,parti,6,1,isp), INFOC(:,parti,6,2,isp)] = calInfo(pBAS.mla);
%             [INFOC(:,parti,7,1,isp), INFOC(:,parti,7,2,isp)] = calInfo(mHeu1.mla);
%             [INFOC(:,parti,8,1,isp), INFOC(:,parti,8,2,isp)] = calInfo(mHeu2.mla);
%             [INFOC(:,parti,9,1,isp), INFOC(:,parti,9,2,isp)] = calInfo(mHeu3.mla);
        
%             % for Pmix
%             mPLBaB = a_Combine_Data(mPLB, mPLaB);
%             lpxPLBaB = a_Combine_Data(lpxPLB, lpxPLaB);
%             count = 0;
%             for imgXID = [20, 6, 30]
%                 for revID = 1:3
%                     [Dm, ind] = DImgVsRev(PLBaB, imgXID, revID);
%                     count = count + 1;
%                     if (isp==1) % these only need to be computed once
%                         [Pmix(:,parti,count,1), Pmix(:,parti,count,2), ~, ~] = a_Get_Human_Perf(Dm);
%                         [BiasMix(:,parti,count,1), BiasMix(:,parti,count,2)] = getHumanChoiceBias(Dm);
%                         [lpxPmix(:,parti,count,1), lpxPmix(:,parti,count,2)] = ...
%                             modelPerfProp(lpxPLBaB.mla(ind,:), Dm.AnswerReal, 'hard');
%                     end
%                     [mPmixC(:,parti,count,1,isp), mPmixC(:,parti,count,2,isp)] = ...
%                         modelPerfProp(mPLBaB.softmla(ind,:), Dm.AnswerReal, 'soft');
%                 end
%             end
%             
%             % for BAS & true anti-BAS renormalized perf
%             if (isp==1)
%                 [PERF(:,parti,5,1), PERF(:,parti,5,2)] = getBaBRenor(Pmix(:,parti,:,:));
%                 [lpxPERF(:,parti,5,1), lpxPERF(:,parti,5,2)] = getBaBRenor(lpxPmix(:,parti,:,:));
%             end
%             [mPERFC(:,parti,5,1,isp), mPERFC(:,parti,5,2,isp)] = getBaBRenor(mPmixC(:,parti,:,:,isp));

        end %iprop
    end %isimu        
    fprintf('Done parti %d.\n', parti);
end %iparti

mPERF = mean(mPERFC, 5);
mPmix = mean(mPmixC, 5);
INFO = mean(INFOC, 5);
%save('PERF.mat','PERF');

end

%%
function [Dm, ind] = DImgVsRev(D, imgXID, revID)
% ImgXID: 20->patch; 6->stripy horizontal; 30->stripy vertical
% RevID:  1 ->patch; 2->stripy horizontal; 3 ->stripy vertical

ind=1:D.Trials;
ind(D.ImageID(:,2)~=imgXID | D.RevealType~=revID)=[];
Dm.Trials=length(ind);
Dm.RevealPosX=D.RevealPosX(ind,:);
Dm.RevealPosY=D.RevealPosY(ind,:);
Dm.RevealTime=D.RevealTime(ind,:);
Dm.StateAnswerTime=D.StateAnswerTime(ind,:);
Dm.StateSearchTime=D.StateSearchTime(ind,:);
Dm.AnswerChoice=D.AnswerChoice(ind,:);
Dm.AnswerReal=D.AnswerReal(ind,:);
Dm.ImageID=D.ImageID(ind,:);
Dm.MaxRevealingTrial=D.MaxRevealingTrial(ind,:);
Dm.RevealType=D.RevealType(ind,:);

end
%% copied from mPerf_comb_prop--bad style!
function [mu, std] = modelPerfProp(mla, ansReal, mode)

    if ( strcmp(mode, 'hard') )
        correct = calCorret(mla, ansReal);
    elseif ( strcmp(mode, 'soft') )
        correct = calProb(mla, ansReal);
    end
    mu = nanmean(correct);  
    std = sqrt(nanvar(correct)./sum(~isnan(correct)));
    
end
%% copied from mPerf_comb_prop--bad style!
function correct = calCorret(m, answ)
    nrev = 25;
    [nrow, ncol] = size(m);
    choice = nan(nrow, ncol);
    
    thres = 1e-3;
    choice(m > 0.5+thres) = 1;
    choice(m < 0.5-thres) = 2;    
    ind = (m > 0.5-thres) & (m < 0.5+thres);
    n = sum(sum(ind==1));
    choice(ind) = round(rand(n,1)) + 1;
    
    ansM = repmat(answ, 1, nrev);
    correct = (choice == ansM) + 0; %convert to double
    correct(isnan(m)) = nan;
end
%% also for modelPerfProp
function correct = calProb(p, answ)
    nrev = 25;
    [nrow, ncol] = size(p);
    choice = nan(nrow, ncol);
    randAns = rand(nrow, ncol);
    
    choice(randAns < p) = 1;
    choice(randAns > p) = 2;    
    
    ansM = repmat(answ, 1, nrev);
    correct = (choice == ansM) + 0; %convert to double
    correct(isnan(p)) = nan;
end

%%
function Dout = addFieldTrials(Din)
    Dout = Din;
    [d1, d2, d3] = size(Din);
    for i1 = 1:d1
        for i2 = 1:d2
            for i3 = 1:d3
                [ntrial, ~] = size(Din(i1,i2,i3).mla);
                Dout(i1,i2,i3).Trials = ntrial;
            end
        end
    end
end
%% copy from a_Get_Human_Performance
function [p,s] = getHumanChoiceBias(D)

    RevN = D.MaxRevealingTrial;
    uni_revn = unique(RevN);
    rn =length(uni_revn);
    
    p = zeros(rn,1);
    s = zeros(rn,1);

    pc = D.AnswerChoice==1;
    sc = D.AnswerChoice==2;

    for i=1:rn
       pTemp = sum(pc(RevN==uni_revn(i)));
       sTemp = sum(sc(RevN==uni_revn(i)));
       tot = pTemp + sTemp;
       p(i) = pTemp/tot;
       s(i) = sTemp/tot;
    end
    
end

%% for BAS + true anti-BAS renormalized
function [mu, sd] = getBaBRenor(PmixSub)
    si = size(PmixSub);
    xn = si(1);
    % BAS involves 1,5,9 in Pmix count
    % true anti-BAS involves 2,3,4,7 in Pmix count
    muM = zeros(xn,7); %mean perf
    sdM = zeros(xn,7); %sd perf
    count = 0;
    for ind = [1,2,3,4,5,7,9]
        count = count+1;
        muM(:,count) = PmixSub(:,1,ind,1);
        sdM(:,count) = PmixSub(:,1,ind,2);
    end
    wM = 1/8*ones(xn,7); %renormalzation weights
    wM(:,1) = 1/4;
    mu = sum(wM.*muM,2);
    sd = sqrt(sum((wM.*sdM).^2,2));
end
%%
function [infoMu, infoSD] = calInfo(mla)
    m = mla;  
    M= -log(0.5) +m.*log(m) +(1-m).*log(1-m);  
    M= -M/log(0.5);
    infoMu = nanmean(M);  
    infoSD = sqrt(nanvar(M)./sum(~isnan(M))); 
end
%% To get the std on the slope, first get the 95%-int by:
% %[b,bint] = regress(PERF(:,4,1,1),[(5:5:25)',ones(5,1)])
% % then to get std, do
% % (bint(1,2) - b(1))/1.96
% % for AL, I got slope = 0.0113, std = 0.0026
% % for BAS, I got slope = 0.0081, std = 0.0030
% %  so t = (0.0113 - 0.0081)/sqrt(0.0026^2+0.0030^2) = 0.8061

%% diagnostic plots
% P = PERF;
% Px = Pmix;
% x = 5:5:25;
% clf;
% % plot PERF
% for parti = 1:4
%     subplot(4,4, 4*(parti-1)+1);
%     hold on;
%     errorbar(x,P(:,parti,2,1),P(:,parti,2,2),'-bs','MarkerFaceColor',[0,0,1],'MarkerSize',5);
%     errorbar(x,P(:,parti,1,1),P(:,parti,1,2),'-ro','MarkerFaceColor',[1,0,0],'MarkerSize',5);
%     %errorbar(x,P(:,parti,3,1),P(:,parti,3,2),'-kd','MarkerFaceColor',[0,0,0],'MarkerSize',5);
%     %errorbar(x,P(:,parti,4,1),P(:,parti,4,2),'-d','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5],'MarkerSize',5);    
%     errorbar(x,P(:,parti,5,1),P(:,parti,5,2),'-kd','MarkerFaceColor',[0,0,0],'MarkerSize',5);
%     axis([0 27 0 1]);
%     title(['P', num2str(parti)]);
%     ylabel('P(correct)');
%     xlabel('Revealing number');
% end
% 
% % plot PERF_ImgVsRev
% for parti = 1:4
%     for img = 1:3
%         for rev = 1:3
%             subplot(4,4, 4*(parti-1)+1+img);
%             hold on;
%             cond = rev + 3*(img-1);
%             % marker for rev Type
%             if (rev==1) sym = '-o'; end
%             if (rev==2) sym = '->'; end
%             if (rev==3) sym = '-^'; end
%             if (img == rev)
%                 cl = [0,0,0];
%             else
%                 cl = 0.5*[1,1,1];
%             end
%             errorbar(x,Px(:,parti,cond,1),Px(:,parti,cond,2), sym, 'color', cl);   
%             axis([0 27 0 1])
%             ylabel('P(correct)');
%             xlabel('Revealing number');
%         end
%     end
% end
% subplot(4,4,2);  title('Img = PA');
% subplot(4,4,3);  title('Img = SH');
% subplot(4,4,4);  title('Img = SV');
% 
% set(gcf,'Color',[1,1,1]);
% 
% %% diagnostic plot for choice biases
% B = Bias;
% Bx = BiasMix;
% x = 5:5:25;
% clf;
% % plot PERF
% for parti = 1:4
%     subplot(4,4, 4*(parti-1)+1);
%     hold on;
%     plot(x,B(:,parti,2,1),'-bs','MarkerFaceColor',[0,0,1],'MarkerSize',5);
%     plot(x,B(:,parti,1,1),'-ro','MarkerFaceColor',[1,0,0],'MarkerSize',5);
%     plot(x,B(:,parti,3,1),'-kd','MarkerFaceColor',[0,0,0],'MarkerSize',5);
%     plot(x,B(:,parti,4,1),'-d','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5],'MarkerSize',5);    
% %     plot(x,B(:,parti,2,2),'--bs','MarkerFaceColor',[0,0,1],'MarkerSize',5);
% %     plot(x,B(:,parti,1,2),'--ro','MarkerFaceColor',[1,0,0],'MarkerSize',5);
% %     plot(x,B(:,parti,3,2),'--kd','MarkerFaceColor',[0,0,0],'MarkerSize',5);
% %     plot(x,B(:,parti,4,2),'--d','Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5],'MarkerSize',5);    
%     axis([0 27 0 1]);
%     title(['P', num2str(parti)]);
%     ylabel('P(choice = PA)');
%     xlabel('Revealing number');
% end
% 
% % plot PERF_ImgVsRev
% for parti = 1:4
%     for img = 1:3
%         for rev = 1:3
%             subplot(4,4, 4*(parti-1)+1+img);
%             hold on;
%             cond = rev + 3*(img-1);
%             % marker for rev Type
%             if (rev==1) sym = 'o'; end
%             if (rev==2) sym = '>'; end
%             if (rev==3) sym = '^'; end
%             if (img == rev)
%                 cl = [0,0,0];
%             else
%                 cl = 0.5*[1,1,1];
%             end
%             plot(x,Bx(:,parti,cond,1), [sym, '-'], 'color', cl);   
%             %plot(x,Bx(:,parti,cond,2), [sym, '--'], 'color', cl);   
%             axis([0 27 0 1])
%             ylabel('P(choice = PA)');
%             xlabel('Revealing number');
%         end
%     end
% end
% subplot(4,4,2);  title('Img = PA');
% subplot(4,4,3);  title('Img = SH');
% subplot(4,4,4);  title('Img = SV');
% 
% set(gcf,'Color',[1,1,1]);
% 
% %% diagnostic plot for model performance
% P = mPERF;
% Px = mPmix;
% % P = lpxPERF;
% % Px = lpxPmix;
% x = 1:25;
% clf;
% % plot PERF
% for parti = 1:4
%     subplot(4,4, 4*(parti-1)+1);
%     hold on;
%     errorbar(x,P(:,parti,2,1),P(:,parti,2,2),'-b','MarkerFaceColor',[0,0,1]);
%     errorbar(x,P(:,parti,1,1),P(:,parti,1,2),'-r','MarkerFaceColor',[1,0,0]);
%     errorbar(x,P(:,parti,3,1),P(:,parti,3,2),'-k','MarkerFaceColor',[0,0,0],'LineWidth',2);
%     errorbar(x,P(:,parti,4,1),P(:,parti,4,2),'-','Color',[0.5,0.5,0.5],'LineWidth',2);    
%     %errorbar(x,P(:,parti,5,1),P(:,parti,5,2),'-k','MarkerFaceColor',[0,0,0],'LineWidth',2);
%     axis([0 27 0 1]);
%     title(['P', num2str(parti)]);
%     ylabel('P(correct)');
%     xlabel('Revealing number');
% end
% 
% % plot PERF_ImgVsRev
% for parti = 1:4
%     for img = 1:3
%         for rev = 1:3
%             subplot(4,4, 4*(parti-1)+1+img);
%             hold on;
%             cond = rev + 3*(img-1);
%             % marker for rev Type
%             if (rev==1) sym = '-o'; end
%             if (rev==2) sym = '->'; end
%             if (rev==3) sym = '-^'; end
%             if (img == rev)
%                 cl = [0,0,0];
%             else
%                 cl = 0.5*[1,1,1];
%             end
%             errorbar(x,Px(:,parti,cond,1),Px(:,parti,cond,2), sym, 'color', cl);   
%             axis([0 27 0 1])
%             ylabel('P(correct)');
%             xlabel('Revealing number');
%         end
%     end
% end
% subplot(4,4,2);  title('Img = PA');
% subplot(4,4,3);  title('Img = SH');
% subplot(4,4,4);  title('Img = SV');
% 
% set(gcf,'Color',[1,1,1]);
% 
% %% plot only average perf both human and model -- this is very much like plot_mperf
% 
% P = PERF;
% x = 5:5:25;
% 
% mP = mPERF;
% mx = 1:25;
% 
% clf;
% % plot PERF
% for parti = 1:4
%     subplot(1, 4, parti);
%     hold on;
%     errorbar(x,P(:,parti,2,1),P(:,parti,2,2),'bs','MarkerFaceColor',[0,0,1]);
%     errorbar(x,P(:,parti,1,1),P(:,parti,1,2),'ro','MarkerFaceColor',[1,0,0]);
%     errorbar(x,P(:,parti,3,1),P(:,parti,3,2),'kd','MarkerFaceColor',[0,0,0]);
%     errorbar(x,P(:,parti,4,1),P(:,parti,4,2),'^','Color',0.5*[1,1,1],'MarkerFaceColor',0.5*[1,1,1]);    
%     plot(mx,mP(:,parti,2,1),'b-');
%     plot(mx,mP(:,parti,1,1),'r-');
%     plot(mx,mP(:,parti,3,1),'k-');
%     plot(mx,mP(:,parti,4,1),'-','Color',0.5*[1,1,1]);    
%     axis([0 27 0.4 1]);
%     title(['P', num2str(parti)]);
%     ylabel('P(correct)');
%     xlabel('Revealing number');
%     axis square
% end
% 
% set(gcf,'Color',[1,1,1]);

