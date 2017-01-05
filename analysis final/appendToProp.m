function appendToProp(path, name)
%example path: 'heuristics_data/heuristicSIMnPROP'
%example name: 'heuristicPROP'

nsimu = 10;
for isimu = 1:nsimu
    load([path, num2str(isimu), '.mat'], 'PROP');
    PROP = rmfield(PROP, 'ml');
    tempPROP = PROP;
    load(['PROP_data/PROP',num2str(isimu),'.mat'], 'PROP');
    nRepl = 3; %for heuristics prop --> 4; for newNoise --> 3 
    for ii = 1:nRepl
        tempPROP(:,ii,:) = PROP(:,ii,:);        
    end
    PROP = tempPROP;
    clear tempPROP;
    eval(sprintf('save appendedPROP%i', isimu));
end

end