function xxD=test_motor_noise()
clf

xxD = zeros(100,2);
xDp = [350,300];
xDn = [350,630];
for ii=1:100
    xxD(ii,:)=MotorNoise(xDp, xDn, 770, 20);
end
plot(xxD(:,1), xxD(:,2), 'r.'); hold on;
plot(mean(xxD(:,1)), mean(xxD(:,2)), 'bo');
plot( [xDp(1), xDn(1)], [xDp(2), xDn(2)], 'k-');
axis([0 770 0 770])
axis square

% simulated noise
dis = norm(xDn-xDp,1)/770*20;
bias = norm(mean(xxD) - xDn,1)/770*20;
xxDcm = xxD/770*20;
covDcm = cov(xxDcm(:,1), xxDcm(:,2));
fprintf('Simu: dis = %f cm; bias = %f cm; x std = %f cm; y std = %f cm;\n',...
    dis, -bias, sqrt(covDcm(1,1)), sqrt(covDcm(2,2)));

% theoretical noise
sig_R = 0.13*dis + 0.29;  
sig_a = 0.011*dis + 0.27;
bias_T = -0.23*dis + 0.26;
fprintf('Theory: dis = %f cm; bias = %f cm; x std = %f cm; y std = %f cm;\n',...
    dis, bias_T, sig_a, sig_R);

end
%%
function xDnew=MotorNoise(xDp,xDn,npix,ss)

    %step1: from pixel to cm
    d=42; %distance from eye to fixation cross in cm
    xDp(1) = (xDp(1)-1)/(npix-1)*ss-ss/2;  xDp(2) = (1-xDp(2))/(npix-1)*ss+ss/2;
    xDn(1) = (xDn(1)-1)/(npix-1)*ss-ss/2;  xDn(2) = (1-xDn(2))/(npix-1)*ss+ss/2;
    %step2: define fixation vectors
    fp=[xDp(1),xDp(2),d]; %fixation previous in x,y,z, cm
    fn=[xDn(1),xDn(2),d]; %fixation now in x,y,x, cm
    %step3: find angle between fixation vectors and motor noise
    % R=acos( (fp*fn')/sqrt(fp*fp')/sqrt(fn*fn') )*180/pi;  %R in degree
    Rcm = norm(xDn-xDp);  % R in cm
    [sig_R, sig_a]=sig_endpoint(Rcm);
    %step4: find difference vector and its perpendicular vector
    dv=fn-fp;  dv(3)=[];  dv=dv/sqrt(dv*dv'); %unit difference vector
    du(1)=-dv(2);  du(2)=dv(1); %unit perpendicular vector
    %step5: compute covariance matrix in x,y
    U=[dv',du']; L=[sig_R,0;0,sig_a].^2;
    CM=U*L*U';
    %step6: sample noise from covariance
    mn = mvnrnd([0,0],CM);
    %step7: add bias
    % errTan: intercept = 0.2643, slope = -0.2284 (all in cm)
    bias = -0.23*Rcm + 0.26;
    bn = bias*dv;
    %step8: from cm to pixel
    xDnew=xDn+mn+bn;
    xDnew=max(-10,xDnew);  xDnew=min(10,xDnew); %bounds
    xDnew(1) = 1+(npix-1)*(xDnew(1)/ss+0.5);
    xDnew(2) = 1-(npix-1)*(xDnew(2)/ss-0.5);
    %mn = mn*npix/ss;
    
end
%%
function [sig_R, sig_a]=sig_endpoint(R)
mode=3;  
if (mode==1) % vanBeer07
    m_R=0.025;  a_R=0.15;   b_R=0.11; %endpoint radius parameters    
    m_a=0.015;  a_a=0.044;  b_a=0.12; %endpoint angle parameters
    sig_R = R*(m_R+a_R*exp(-b_R*R));
    sig_a = R*(m_a+a_a*exp(-b_a*R));    
elseif (mode==2)
    % Krappmann98 + vanBeers07 
    % for amp: 0.03*R + 0.23 + 1.51
    sig_R = (20/26.8)*(0.03*(R*26.8/20) + 0.23 + 1.51);  
    % for dir: 0.014*R + 0.085 + 0.626
    sig_a = (20/26.8)*(0.014*(R*26.8/20) + 0.085 + 0.626);
elseif (mode==3)
    % from own experimental measurement (both input and output is in cm)
    % sdTan: intercept = 0.3345, slope = 0.1011 (all in cm)
    % sdNor: intercept = 0.2664, slope = 0.004257 (all in cm)
    % sig_R = 0.10*R + 0.33;  
    % sig_a = 0.0043*R + 0.27;
    sig_R = 0.13*R + 0.29;  
    sig_a = 0.011*R + 0.27;
end
end
