function [REV, RrevMap, RdrevMap] = revmap_sim_heuristics(SIM)
%% make a matrix of structure to store all data
% REV has 3 layers: 
% patterns x4 (PA,SH,SV,ALL); 
% participants x4 (SC,EP,BZ,AVG);
% conditions x2 (AL, heuristics)

if (nargin==0)
    nSIM = 10;
    load('heuristics_data/gungho/heuristicSIMnPROP1.mat', 'SIM');
    SIMC = SIM;
    for iSIM = 2:nSIM
        load(['heuristics_data/gungho/heuristicSIMnPROP',num2str(iSIM),'.mat'], 'SIM');
        for irow = 1:3
            for jcol = 2
                SIMC(irow,jcol) = a_Combine_Data(SIMC(irow,jcol), SIM(irow,jcol));
            end
        end
    end
    SIM = SIMC;
end

heuInd = 9; %2 is mBAS, 4 is hueristics, 9?

load('SCdata.mat', 'SCAL');  
REV(1,1,1)=a_DimageID(SCAL,20);
REV(2,1,1)=a_DimageID(SCAL,6);
REV(3,1,1)=a_DimageID(SCAL,30);  
REV(1,1,2)=a_DimageID(SIM(1,heuInd),20);   
REV(2,1,2)=a_DimageID(SIM(1,heuInd),6);   
REV(3,1,2)=a_DimageID(SIM(1,heuInd),30);  

load('EPdata.mat', 'EPAL');
REV(1,2,1)=a_DimageID(EPAL,20);
REV(2,2,1)=a_DimageID(EPAL,6);
REV(3,2,1)=a_DimageID(EPAL,30);  
REV(1,2,2)=a_DimageID(SIM(2,heuInd),20);   
REV(2,2,2)=a_DimageID(SIM(2,heuInd),6);   
REV(3,2,2)=a_DimageID(SIM(2,heuInd),30);  
load('BZdata.mat', 'BZAL');
REV(1,3,1)=a_DimageID(BZAL,20);
REV(2,3,1)=a_DimageID(BZAL,6);
REV(3,3,1)=a_DimageID(BZAL,30);  
REV(1,3,2)=a_DimageID(SIM(3,heuInd),20);   
REV(2,3,2)=a_DimageID(SIM(3,heuInd),6);   
REV(3,3,2)=a_DimageID(SIM(3,heuInd),30);  

% combining pattern as the 4th pattern
for iparti=1:3
for icond=1:2
    zz1=REV(1,iparti,icond);    zz2=REV(2,iparti,icond);    zz3=REV(3,iparti,icond);
    REV(4,iparti,icond)=a_Combine_Data(zz1,zz2,zz3);    
end
end
% combining participants, the 4th participant
for ipattern=1:4
for icond=1:2
    REV(ipattern,4,icond)=a_Combine_Data(REV(ipattern,1,icond),REV(ipattern,2,icond),REV(ipattern,3,icond));    
end
end
%save('REVf.mat','REV');
fprintf('REV assignment complete.\n');

% ###################################################
% caluclate revealing densities (difference and mean)
% define gaussian filter
npix=770; ImageSize=20; bin=1;
GaussFilter=a_GaussFilter(npix);

RrevMap = zeros(770,770,4,4,2,25);
RdrevMap = zeros(770,770,3,4,2,25);
for rev = 1:25
    revMap=zeros(770,770,4,4,2);
    for icond=1:2
    for iparti=1:4
    for ipattern=1:4
        zz=REV(ipattern,iparti,icond);
        %revMap(:,:,i,j,k) = a_Get_RevealMap_mod2(zz,1:rev,1:zz.Trials,GaussFilter,770,2); %for normal map
        revMap(:,:,ipattern,iparti,icond) = a_Get_RevealMap_mod2(zz,1:25,1:zz.Trials,GaussFilter,770,3); %for center removed map
    %     imagesc(revMap(:,:,i,j,k));
    %     pause
    end
    end
    end

    %drevMap has 5 layers: 770; 770; pattern; parti; cond
    % caluclate difference densities
    drevMap=zeros(770,770,3,4,2);
    for icond=1:2
    for iparti=1:4
    for ipattern=1:3
        drevMap(:,:,ipattern,iparti,icond) = revMap(:,:,ipattern,iparti,icond)-revMap(:,:,4,iparti,icond);
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