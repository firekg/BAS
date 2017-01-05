function [PROP] = prop_sim(SIM)
%% PROP has 2 layers: participants x3 (SC,EP,BZ); conditions x5 (ALprop, PLRprop, BASprop, mBASprop, nBASprop)
    load('lpx.mat');
    load('expDataWPix.mat');
    
    %pv= [  lsigo,    phit,   lphis, phis0,        lbx,        lbd, lapse, prpa, window, prxpa, prCommon, shift];
    pv =[log(1.0),   nan,     nan,   nan,       log(0),        nan,     0,  1/2,      1,   1/3,      1/3,    12];
    %pw=[lbm,  km,  nm, bias];
    pw=[   0, nan, nan,    0]; % none of the last 3 are used here
    
    % make probe points z
    npix=770;
    %ImageSize=20;
    n=110;
    xx = linspace(1,npix,n);
    yy = linspace(1,npix,n);
    [x1,x2]=meshgrid(xx,yy);
    z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
    %dxx = xx(2)-xx(1);
        
    for i = 1:3 %3
        if(i==1)
            AL = SCALwPix;
            PLR = SCPLRwPix;  
            PLB = SCPLBwPix;  
            PLaB = SCPLaBwPix;
            pv(1) = log(0.5);
            pv(6) = 0.34901995;
            pv(7) = 0.04471412;
            pv(12) = 16;
        elseif(i==2)  
            AL = EPALwPix;  
            PLR = EPPLRwPix;  
            PLB = EPPLBwPix;  
            PLaB = EPPLaBwPix;  
            pv(1) = log(0.5);
            pv(6) = 0.637376715;
            pv(7) = 0.120818875;
            pv(12) = 17;
        elseif(i==3)  
            AL = BZALwPix;  
            PLR = BZPLRwPix;  
            PLB = BZPLBwPix;  
            PLaB = BZPLaBwPix;  
            pv(1) = log(0.3);
            pv(6) = 0.374908933;
            pv(7) = 0.101003052;
            pv(12) = 15;
        end
        fprintf('Done loading parti %d.\n', i);
        for task = 1:5 %20
            Gmap = zeros(AL.Trials,25,3); % zero makes lpx of no effect in a_Model 

%             PROP(i,1,task) = a_Get_BALDscoreProp(AL,[],pv,pw,randn(1000,25),Gmap,z,2);  %AL prop            
%             %Gmap = lpx(i,1).lik; not used anyway now
%             PROP(i,2,task) = a_Get_BALDscoreProp(PLR,[],pv,pw,randn(1000,25),Gmap,z,2); %PLR prop            
%             %Gmap = lpx(i,2).lik;
%             PROP(i,3,task) = a_Get_BALDscoreProp(PLB,[],pv,pw,randn(1000,25),Gmap,z,2); %PLB prop
%             %Gmap = lpx(i,3).lik;
%             PROP(i,4,task) = a_Get_BALDscoreProp(PLaB,[],pv,pw,randn(1000,25),Gmap,z,2);  %PLaB prop
%             PROP(i,5,task) = a_Get_BALDscoreProp(SIM(i,1),[],pv,pw,randn(1000,25),Gmap,z,2);  %posterior-indep & order-dep
%             PROP(i,6,task) = a_Get_BALDscoreProp(SIM(i,4),[],pv,pw,randn(1000,25),Gmap,z,2);  %posterior-dep & order-indep
%             PROP(i,7,task) = a_Get_BALDscoreProp(SIM(i,5),[],pv,pw,randn(1000,25),Gmap,z,2);  %posterior-dep & order-dep
%             PROP(i,8,task) = a_Get_BALDscoreProp(SIM(i,2),[],pv,pw,randn(1000,25),Gmap,z,2);  %nBAS prop            
%             %PROP(i,9,task) = a_Get_BALDscoreProp(SIM(i,3),[],pv,pw,randn(1000,25),Gmap,z,2);  %pBAS prop
            
            % for checking stuff
%             PROP(i,1,task) = a_Get_BALDscoreProp(PLB,[],pv,pw,randn(1000,25),Gmap,z,2); %PLB prop
%             PROP(i,2,task) = a_Get_BALDscoreProp(SIM(i,1),[],pv,pw,randn(1000,25),Gmap,z,2);  %ideal BAS with prior bias
%             PROP(i,3,task) = a_Get_BALDscoreProp(SIM(i,2),[],pv,pw,randn(1000,25),Gmap,z,2);  %ideal BAS without prior bias

            % 2015-09-22 shuffled AL and PLB
            AL25 = a_DRevNum(AL, 25);
            PROP(i,1,task) = a_Get_BALDscoreProp(AL25,[],pv,pw,randn(1000,25),Gmap,z,2);  %AL prop
            AL25s = shuffleData(AL25);
            PROP(i,2,task) = a_Get_BALDscoreProp(AL25s,[],pv,pw,randn(1000,25),Gmap,z,2);  %AL-shuffled prop
            PLB25 = a_DRevNum(PLB, 25);
            PROP(i,3,task) = a_Get_BALDscoreProp(PLB25,[],pv,pw,randn(1000,25),Gmap,z,2);  %PLB prop
            PLB25s = shuffleData(PLB25);
            PROP(i,4,task) = a_Get_BALDscoreProp(PLB25s,[],pv,pw,randn(1000,25),Gmap,z,2);  %PLB-shuffled prop

            fprintf('Done parti %d task %d.\n', i, task);
        end
    end

end

%% 2015-07-02 with or without Prior bias or on ideal BAS
%             PROP(i,1,task) = a_Get_BALDscoreProp(PLB,[],pv,pw,randn(1000,25),Gmap,z,2); %PLB prop
%             PROP(i,2,task) = a_Get_BALDscoreProp(SIM(i,1),[],pv,pw,randn(1000,25),Gmap,z,2);  %ideal BAS with prior bias
%             PROP(i,3,task) = a_Get_BALDscoreProp(SIM(i,2),[],pv,pw,randn(1000,25),Gmap,z,2);  %ideal BAS without prior bias            


%% older version -- 2015-04-22
% function [PROP] = prop_sim(SIM)
% %% PROP has 2 layers: participants x3 (SC,EP,BZ); conditions x5 (ALprop, PLRprop, BASprop, mBASprop, nBASprop)
%     load('lpx.mat');
%     load('expDataWPix.mat');
%     %pv= [ lsigo,  phit,   lphis, phis0,        lbx,        lbd, lapse, prpa, w];
%     pv =[log(1),   nan,     nan,   nan,        log(1),        nan,     0,  0.5, 1];
%     %pw=[lbm,  km,  nm, bias];
%     pw=[   0, nan, nan,    0]; % none of the last 3 are used here
%     % make probe points z
%     npix=770;  ImageSize=20;  n=110;  xx = linspace(1,npix,n);  yy = linspace(1,npix,n);
%     [x1,x2]=meshgrid(xx,yy);  z=[reshape(x1,n*n,1),reshape(x2,n*n,1)]; % z in pixel x,y coordinate
%     dxx = xx(2)-xx(1);
%         
%     for i=1:3 %3
%         if(i==1)      
%             %load('SCdata.mat', 'SCAL', 'SCPLRbt', 'SCPLBbt', 'SCPLaBbt');  
%             AL=SCAL;
%             PLR=SCPLRbt;  
%             PLB=SCPLBbt;  
%             PLaB=SCPLaBbt;  
%             DIMAL=a_DIM(SCAL);
%             DIMPLR=a_DIM(SCPLRbt);
%             DIMPLB=a_DIM(SCPLBbt);
%             DIMPLaB=a_DIM(SCPLaBbt);
%             %pv(1) = log(1.1);  % from max likelihood fit
%             %pv(1) = log(0.0687);  % from color difference exp
%             %pv(1) = log(0.17);  % used to generate revealing locations
%             pv(1) = log(1.3);  % from max likelihood fit to AL only
%         elseif(i==2)  
%             %load('EPdata.mat', 'EPAL', 'EPPLRbt', 'EPPLBbt', 'EPPLaBbt');  
%             AL=EPAL;  
%             PLR=EPPLRbt;  
%             PLB=EPPLBbt;  
%             PLaB=EPPLaBbt;  
%             DIMAL=a_DIM(EPAL);  
%             DIMPLR=a_DIM(EPPLRbt);
%             DIMPLB=a_DIM(EPPLBbt);
%             DIMPLaB=a_DIM(EPPLaBbt);
%             %pv(1) = log(1.2);  % from max likelihood fit
%             %pv(1) = log(0.0687);  % from color difference exp
%             %pv(1) = log(0.17);  % used to generate revealing locations
%             pv(1) = log(0.9);  % from max likelihood fit to AL only
%         elseif(i==3)  
%             %load('BZdata.mat', 'BZAL', 'BZPLRbt', 'BZPLBbt', 'BZPLaBbt');  
%             AL=BZAL;  
%             PLR=BZPLRbt;  
%             PLB=BZPLBbt;  
%             PLaB=BZPLaBbt;  
%             DIMAL=a_DIM(BZAL);  
%             DIMPLR=a_DIM(BZPLRbt);
%             DIMPLB=a_DIM(BZPLBbt);
%             DIMPLaB=a_DIM(BZPLaBbt);
%             %pv(1) = log(0.9); % from max likelihood fit
%             %pv(1) = log(0.0687);  % from color difference exp
%             %pv(1) = log(0.17);  % used to generate revealing locations
%             pv(1) = log(0.9);  % from max likelihood fit to AL only
%         end
%         disp(sprintf('Done loading parti %d.', i));
%         for task=1 %20
%             Gmap = zeros(AL.Trials,3); % zero makes lpx of no effect in a_Model 
%             PROP(i,1,task) = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %AL prop
%             Gmap = lpx(i,1).lik;
%             PROP(i,2,task) = a_Get_BALDscoreProp(PLR,DIMPLR,pv,pw,randn(1000,25),Gmap,z,2); %PLR prop
%             %PROP(i,XXX) = a_Get_BALDscoreProp(SIM(i,1),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %BAS simu prop
%             Gmap = lpx(i,2).lik;
%             PROP(i,3,task) = a_Get_BALDscoreProp(PLB,DIMPLB,pv,pw,randn(1000,25),Gmap,z,2); %PLB prop
%             %PROP(i,4,task) = a_Get_BALDscoreProp(SIM(i,2),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %BAS + motor noise prop
%             
%             %ffAL = feedForwardRevPos(AL);
%             %PROP(i,5,task) = a_Get_BALDscoreProp(ffAL,DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %posterior-indep & order-dep
%             %PROP(i,6,task) = a_Get_BALDscoreProp(SIM(i,4),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %posterior-dep & order-indep
%             %PROP(i,7,task) = a_Get_BALDscoreProp(SIM(i,5),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %posterior-dep & order-dep
%             %PROP(i,8,task) = a_Get_BALDscoreProp(SIM(i,6),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);
%             %PROP(i,9,task) = a_Get_BALDscoreProp(SIM(i,7),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);
%             %PROP(i,10,task) = a_Get_BALDscoreProp(SIM(i,8),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);
%             %PROP(i,11,task) = a_Get_BALDscoreProp(SIM(i,9),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);
%             
%             %PROP(i,12, task) = a_Get_BALDscoreProp(SIM(i,3),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %BAS + motor + selection noise prop
%             Gmap = lpx(i,3).lik;
%             PROP(i,13, task) = a_Get_BALDscoreProp(PLaB,DIMPLaB,pv,pw,randn(1000,25),Gmap,z,2);  %PLaB prop
%             disp(sprintf('Done parti %d task %d.', i, task));
%         end
%     end
% 
% end
%% Attempts for parfor
%         parfor task=1:20
%             A1{task} = a_Get_BALDscoreProp(AL,DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %AL prop
%             %disp(sprintf('Done AL.'));
%             A2{task} = a_Get_BALDscoreProp(PLR,DIMPLR,pv,pw,randn(1000,25),Gmap,z,2); %PLR prop
%             %disp(sprintf('Done Random.'));
%             %PROP(i,3) = a_Get_BALDscoreProp(SIM(i,1),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %BAS simu prop
%             A3{task} = a_Get_BALDscoreProp(PLB,DIMPLB,pv,pw,randn(1000,25),Gmap,z,2); %PLB prop
%             %disp(sprintf('Done BAS.'));
%             A4{task} = a_Get_BALDscoreProp(SIM(i,2),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %BAS + motor noise prop
%             %disp(sprintf('Done mBAS.'));
%             %PROP(i,5) = a_Get_BALDscoreProp(SIM(i,3),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %BAS + motor + selection noise prop
%             %PROP(i,5,task) = PROP(i,4);
%             A5{task} = a_Get_BALDscoreProp(SIM(i,4),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %BAS + motor noise + softmax log weights
%             %disp(sprintf('Done nBAS.'));
%             %disp(sprintf('parti %d PROP done.',i));
%             disp(sprintf('Done parti %d task %d.', i, task));
%         end
%         for task=1:20
%             PROP(i,1,task) = A1{task};
%             %disp(sprintf('Done AL.'));
%             PROP(i,2,task) = A2{task};
%             %disp(sprintf('Done Random.'));
%             %PROP(i,3) = a_Get_BALDscoreProp(SIM(i,1),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %BAS simu prop
%             PROP(i,3,task) = A3{task};
%             %disp(sprintf('Done BAS.'));
%             PROP(i,4,task) = A4{task};
%             %disp(sprintf('Done mBAS.'));
%             %PROP(i,5) = a_Get_BALDscoreProp(SIM(i,3),DIMAL,pv,pw,randn(1000,25),Gmap,z,2);  %BAS + motor + selection noise prop
%             %PROP(i,5,task) = PROP(i,4);
%             PROP(i,5,task) = A5{task};
%             %disp(sprintf('Done nBAS.'));
%             %disp(sprintf('parti %d PROP done.',i));            
%         end

