function [REV, RrevMap, RdrevMap] = revmap_sim(SIM)
%% make a matrix of structure to store all data
% % REV has 3 layers: 
% % patterns x4 (PA,SH,SV,ALL); 
% % participants x4 (SC,EP,BZ,AVG);
% % conditions x4 (AL, BAS, mnBAS or pnBAS)
% ncond = 2;

% if 0
%     nSIM = 10;
%     %load(['SIM_data/SIM1.mat']);
%     %load(['newNoisePROP_data/newNoisePROP1.mat']);
%     SIMC = SIM(:,:,1);
%     for iSIM = 2:nSIM
%         %load(['SIM_data/SIM',num2str(iSIM),'.mat']);
%         %load(['newNoisePROP_data/newNoisePROP',num2str(iSIM),'.mat']);
%         SIMtemp = SIM(:,:,iSIM);
%         for irow = 1:3
%             for jcol = 2
%                 SIMC(irow,jcol) = a_Combine_Data(SIMC(irow,jcol), SIMtemp(irow,jcol));
%             end
%         end
%     end
%     SIM = SIMC;
% end

% mnBAS = 2; %2 is mBAS, 3 is nBAS
% % pnBAS = 3; % motor noise + prob matching
% FlagCorrect = 1;

% load('SCdata.mat', 'SCAL', 'SCPLBbt');  
% REV(1,1,1)= a_DimageID(SCAL,20);
% REV(2,1,1)= a_DimageID(SCAL,6);
% REV(3,1,1)= a_DimageID(SCAL,30);  
% REV(1,1,2)=a_DimageID(SIM(1,mnBAS),20);   
% REV(2,1,2)=a_DimageID(SIM(1,mnBAS),6);   
% REV(3,1,2)=a_DimageID(SIM(1,mnBAS),30);  

% % <2015-09-10 correct-incorrect
% [REV(1,1,1), rmInd1] = correctTrials( a_DimageID(SCAL,20), FlagCorrect);
% [REV(2,1,1), rmInd2]= correctTrials( a_DimageID(SCAL,6), FlagCorrect);
% [REV(3,1,1), rmInd3]= correctTrials( a_DimageID(SCAL,30), FlagCorrect);
% REV(1,1,2) = removeTrials( a_DimageID(SIM(1,mnBAS),20), rmInd1);   
% REV(2,1,2) = removeTrials( a_DimageID(SIM(1,mnBAS),6),  rmInd2);   
% REV(3,1,2) = removeTrials( a_DimageID(SIM(1,mnBAS),30), rmInd3);   
% %>

% load('EPdata.mat', 'EPAL', 'EPPLBbt');
% REV(1,2,1)= a_DimageID(EPAL,20);
% REV(2,2,1)= a_DimageID(EPAL,6);
% REV(3,2,1)= a_DimageID(EPAL,30);  
% REV(1,2,2)=a_DimageID(SIM(2,mnBAS),20);   
% REV(2,2,2)=a_DimageID(SIM(2,mnBAS),6);   
% REV(3,2,2)=a_DimageID(SIM(2,mnBAS),30);  

% % <2015-09-10 correct-incorrect
% [REV(1,2,1), rmInd1] = correctTrials( a_DimageID(EPAL,20), FlagCorrect);
% [REV(2,2,1), rmInd2] = correctTrials( a_DimageID(EPAL,6), FlagCorrect);
% [REV(3,2,1), rmInd3] = correctTrials( a_DimageID(EPAL,30), FlagCorrect);  
% REV(1,2,2) = removeTrials( a_DimageID(SIM(2,mnBAS),20), rmInd1);   
% REV(2,2,2) = removeTrials( a_DimageID(SIM(2,mnBAS),6),  rmInd2);   
% REV(3,2,2) = removeTrials( a_DimageID(SIM(2,mnBAS),30), rmInd3);  
% %>

% load('BZdata.mat', 'BZAL', 'BZPLBbt');
% REV(1,3,1)= a_DimageID(BZAL,20);
% REV(2,3,1)= a_DimageID(BZAL,6);
% REV(3,3,1)= a_DimageID(BZAL,30);  
% REV(1,3,2)=a_DimageID(SIM(3,mnBAS),20);   
% REV(2,3,2)=a_DimageID(SIM(3,mnBAS),6);   
% REV(3,3,2)=a_DimageID(SIM(3,mnBAS),30);  

% % <2015-09-10 correct-incorrect
% [REV(1,3,1), rmInd1] = correctTrials( a_DimageID(BZAL,20), FlagCorrect);
% [REV(2,3,1), rmInd2] = correctTrials( a_DimageID(BZAL,6), FlagCorrect);
% [REV(3,3,1), rmInd3] = correctTrials( a_DimageID(BZAL,30), FlagCorrect);  
% REV(1,3,2) = removeTrials( a_DimageID(SIM(3,mnBAS),20), rmInd1);   
% REV(2,3,2) = removeTrials( a_DimageID(SIM(3,mnBAS),6),  rmInd2);   
% REV(3,3,2) = removeTrials( a_DimageID(SIM(3,mnBAS),30), rmInd3);  
% %>

%% combining pattern as the 4th pattern
% for parti=1:3
% for k=1:ncond
%     zz1=REV(1,parti,k);    zz2=REV(2,parti,k);    zz3=REV(3,parti,k);
%     REV(4,parti,k)=a_Combine_Data(zz1,zz2,zz3);    
% end
% end
% % combining participants, the 4th participant
% for i=1:4
% for k=1:ncond
%     REV(i,4,k)=a_Combine_Data(REV(i,1,k),REV(i,2,k),REV(i,3,k));    
% end
% end
% %save('REVf.mat','REV');
% fprintf('REV assignment complete.\n');
% % RrevMap = [];
% % RdrevMap = [];


%% 2015-10-26 rewritten
% REV has 3 layers: 
% patterns x4 (PA,SH,SV,ALL); 
% participants x5 (SC,EP,BZ,AJ,AT);
% conditions x2 (AL,BAS)

% mnBAS = 2;
% nparti = 4;
% ncond = 2;
% FlagCorrect = 0;
% 
% % for main paper
% for cond = 1:ncond
%     for parti = 1:nparti-1
%         if parti ==1
%             load('SCdata.mat', 'SCAL');
%             D = SCAL;
%         elseif parti == 2
%             load('EPdata.mat', 'EPAL');
%             D = EPAL;
%         elseif parti == 3
%             load('BZdata.mat', 'BZAL');
%             D = BZAL;
%         end
%         [D, rmInd] = correctTrials(D, FlagCorrect); %for correct/incorrect
%         if cond==2 %overwrite
%             D = SIM(parti,mnBAS);
%             %D = removeTrials(SIM(parti,mnBAS), rmInd); %for correct/incorrect
%         end
%         Drp = DRemovePhase(D);
%         REV(1,parti,cond)= a_DimageID(Drp,20);
%         REV(2,parti,cond)= a_DimageID(Drp,6);
%         REV(3,parti,cond)= a_DimageID(Drp,30);
%     end
% end

%% 2015-10-28 for rescanning control
nparti = 4;
ncond = 2;

% for rescanning control
for cond = 1:ncond
    for parti = 1:nparti-1
        if cond == 1
            if parti == 1
                load('AJdata.mat', 'AJAL');
                D = AJAL;
            elseif parti == 2
                load('ATdata.mat', 'ATAL');
                D = ATAL;
            elseif parti == 3
                load('KHdata.mat', 'KHAL');
                D = KHAL;
            end
        elseif cond==2
            if parti == 1
                load('SCdata.mat', 'SCAL');
                D = SCAL;
            elseif parti == 2
                load('EPdata.mat', 'EPAL');
                D = EPAL;
            elseif parti == 3
                load('BZdata.mat', 'BZAL');
                D = BZAL;
            end
        end
        Drp = DRemovePhase(D);
        REV(1,parti,cond)= a_DimageID(Drp,20);
        REV(2,parti,cond)= a_DimageID(Drp,6);
        REV(3,parti,cond)= a_DimageID(Drp,30);
    end
end


%% combining pattern as the 4th pattern
for parti = 1:nparti-1
    for cond = 1:ncond
        zz1=REV(1,parti,cond);
        zz2=REV(2,parti,cond);
        zz3=REV(3,parti,cond);
        REV(4,parti,cond) = a_Combine_Data(zz1,zz2,zz3);
    end
end

% combining participants, the 4th participant
for i=1:4
for k=1:ncond
    REV(i,nparti,k)=a_Combine_Data(REV(i,1,k),REV(i,2,k),REV(i,3,k));    
end
end

fprintf('REV assignment complete.\n');

%% ###################################################
% caluclate revealing densities (difference and mean)
% define gaussian filter
npix=770; ImageSize=20; bin=1;
GaussFilter=a_GaussFilter(npix);

RrevMap = zeros(770,770,4,nparti,ncond,25);
RdrevMap = zeros(770,770,3,nparti,ncond,25);
for rev = 1:25
    revMap=zeros(770,770,4,nparti,ncond);
    for k=1:ncond
    for parti=1:nparti
    for i=1:4
        zz=REV(i,parti,k);
        %revMap(:,:,i,parti,k) = a_Get_RevealMap_mod2(zz,1:rev,1:zz.Trials,GaussFilter,770,2); %for normal map
        %revMap(:,:,i,parti,k) = a_Get_RevealMap_mod2(zz,1:25,1:zz.Trials,GaussFilter,770,2); %for non-center-removed map
        %revMap(:,:,i,parti,k) = a_Get_RevealMap_mod2(zz,1:25,1:zz.Trials,GaussFilter,770,3); %for center-removed map
        %revMap(:,:,i,j,k) = a_Get_RevealMap_mod2(zz,1:rev,1:zz.Trials,GaussFilter,770,3); %for center-removed order map
        revMap(:,:,i,parti,k) = a_Get_RevealMap_mod2(zz,rev,1:zz.Trials,GaussFilter,770,2); %for fixation-specific maps
        %revMap(:,:,i,j,k) = a_Get_RevealMap_mod2(zz,rev,1:zz.Trials,GaussFilter,770,3); %for center-removed fixation-specific maps
    %     imagesc(revMap(:,:,i,j,k));
    %     pause
    end
    end
    end

    %drevMap has 5 layers: 770; 770; pattern; parti; cond
    % caluclate difference densities
    drevMap=zeros(770,770,3,nparti,ncond);
    for k=1:ncond
    for parti=1:nparti
    for i=1:3
        drevMap(:,:,i,parti,k) = revMap(:,:,i,parti,k)-revMap(:,:,4,parti,k);
        %imagesc(drevMap(:,:,i,j,k));
        %pause
    end
    end
    end

    % scale maps to [-1,1]
    drevV=abs(drevMap(:));
    fprintf('max diff map = %f.\n', max(drevV));
    drevMap=drevMap/max(drevV);

    revV=abs(revMap(:));
    fprintf('max mean map = %f.\n', max(revV));
    revMap=revMap/max(revV);
    
    RrevMap(:,:,:,:,:,rev) = revMap;
    RdrevMap(:,:,:,:,:,rev) = drevMap;
    fprintf('revMap %d/25.\n', rev);
end
    

end