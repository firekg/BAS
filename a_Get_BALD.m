function D=a_Get_BALD(type,xsf,ysf,istart,iend)
%%
% type=2;
% xsf=30;
% ysf=6;
% istart=10;
% iend=10;

% see a_GenGPimage for input meaning

% scrsz = get(0,'ScreenSize');
% hf=figure('Name', 'Simulation Plot Window', 'position',[scrsz(3)/10 scrsz(4)/4 scrsz(3)/1.2 scrsz(4)/2],'menubar','none');

% npix=770;  pin=[ npix/20, npix/20, npix/6, npix/30, log(0.5)];  w=7;  n=70;
npix=770;  pin=[ npix/20, npix/20, npix/6, npix/30, log(0.5)];  w=3;  n=110;  nr=25;  %px=linspace(1,npix,n);
% npix=770;  pin=[ npix/20, npix/20, npix/6, npix/30, log(0.1)];  w=1;  n=155;

D.pin=pin;
D.w=w;
D.n=n;

% record mla, mlb, BALDscore, BALDx, BALDy
d0=iend-istart+1; % d0 = iamge ID;  nr = revealing number
D.W_mla       = zeros(nr,d0);
D.W_mlb       = zeros(nr,d0);
D.W_BALDscore = zeros(nr,d0);
D.W_BALDx     = zeros(nr,d0);
D.W_BALDy     = zeros(nr,d0);

for im=istart:iend
    TT = load(strcat('im/Type',num2str(type),'X',num2str(xsf),'Y',num2str(ysf),'Num', num2str(im),'raw.mat'));
    yimage = TT.yimage;
    BALDxD = [0,0] + npix/2 + 1;  % in pixel
    % simulation loop starts
    tol=nr;
    count=1;
    while ( count<tol )    
        count=count+1;    
        % BALD caculation
        [llslpa, llsls, H_BALD, z] = a_Get_BALDmap_window(BALDxD, yimage, pin, w, 0.5, 1, 0, n, [1,npix], [1,npix]);
        [C, ind] = max(H_BALD);
        BALDxnew = z(ind,:); %z is in (x,y)
        BALDxD = [BALDxD; BALDxnew];
        %record
        D.W_mla(count-1,im)     = exp(llslpa);
        D.W_mlb(count-1,im)     = exp(llsls);
        D.W_BALDscore(count,im) = C;
        % reshape is to break rows into n subrows and put side by side
        % z: [hor row 1 of image; hor row 2 of image; ...] --> so is H_BALD
        
%         % plots
%         colormap(gray);
% 
%          subplot(1,2,1); imagesc(px, px, reshape(H_BALD,n,n)); hold on;         
%             plot( BALDxD(:,1), BALDxD(:,2), 'r.'); 
%             plot( BALDxnew(1), BALDxnew(2), 'ro','MarkerFaceColor','g', 'MarkerSize',5); hold off;
% 
%          subplot(1,2,2); imagesc(yimage); hold on;         
%             plot( BALDxD(:,1), BALDxD(:,2), 'r.'); 
%             plot( BALDxnew(1), BALDxnew(2), 'ro','MarkerFaceColor','g', 'MarkerSize',5);  hold off;
% 
%         % allows pressing q to exit the loop
%         if strcmp(get(hf,'currentcharacter'),'q')
%             break
%         end
% 
%         drawnow
%         pause
        
    end
    D.W_BALDx(:,im) = BALDxD(:,1);
    D.W_BALDy(:,im) = BALDxD(:,2);
    im
end

% from pixel to cm for display in experiment
% row[1]=10 y; column[1]= -10 x;
ss=20;  % screen size = [-10,10] cm
D.W_BALDx = (D.W_BALDx-1)/(npix-1)*ss-ss/2;
D.W_BALDy = (1-D.W_BALDy)/(npix-1)*ss+ss/2;

end
