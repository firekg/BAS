function [val, dval] = optfunc_lpr(x,lpr)

lbd = x(1);  bd = exp(lbd);
lapse = x(2);

softlik = 1./(1+exp(-bd*lpr));
lapse_softlik = (1-lapse)*softlik + 0.5*lapse;

dlog = 1./sum(lapse_softlik, 4);
dlbd = (1-lapse).*softlik./(1+exp(bd*lpr)).*lpr;
dval_lbd = dlog.*sum(dlbd, 4); % gradient wrt lapse
dlapse =  0.5-softlik;
dval_lapse = dlog.*sum(dlapse, 4); % gradient wrt lapse

% first marginalize over simulation (4th layer in optimize_para)
% output1: negative sum log marginalized posterior
val = -sum(log(mean(lapse_softlik,4)));
% output2: gradient of sum of log of sum
dval(1) = -bd*sum(dval_lbd);
dval(2) = -sum(dval_lapse);

% no marginalization: when input is lpr(:):
% temp = (1-lapse)./lapse_softlik.*softlik./(1+exp(bd*lpr)); % precompute
% dval_lbd = temp.*lpr; % gradient wrt bd
% dval_lapse = (0.5-softlik)./lapse_softlik; % gradient wrt lapse

% output1: negative sum log posterior of all trials
% val = -sum(log(lapse_softlik));
% output2: gradient of sum of log
% dval(1) = -bd*sum(dval_lbd);
% dval(2) = -sum(dval_lapse);


end