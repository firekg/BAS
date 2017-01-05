function [opt_para_val, x0_para] = opt_lpr(lpr)

    %initialize
    seedn = 1; % number of random starting points
    %lb = [lbd, lapse];
    % no lapse
%     lb = [log(0.01),  0]; % lower bounds
%     ub = [log(10),  1e-5]; % upper bounds
    % yes lapse
    lb = [log(0.01),  0]; % lower bounds
    ub = [log(10),    1]; % upper bounds
    options = optimset('Algorithm','interior-point');
    options = optimset(options,'Display','off');
    options = optimset(options,'GradObj','on');
    options = optimset(options,'TolX',1e-6,'TolFun',1e-6);    
    
    %start optimize
    opt_para_val = zeros(seedn,3);
    x0_para = zeros(seedn,2);
    for i = 1:seedn
        x0 = rand(1,2).*(ub-lb)+lb;
        f = @(x)optfunc_lpr(x, lpr);
        [x,fval,exitflag,output] = fmincon(f,x0,[],[],[],[],lb,ub,[],options);
        opt_para_val(i,:) = [x,fval];
        x0_para(i,:) = [x0];
    end

end


