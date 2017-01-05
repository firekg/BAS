function logNorMap = norRevMap(RrevMap)

    npatt = 3;
    nrev = 25;
    cond = 2; %BAS
    parti = 4;
    norMap = zeros(770,770,npatt,nrev);

    for rev = 1:nrev
        for patt = 1:npatt
            temp = RrevMap(:,:,patt,parti,cond,rev);
            norFac = sum(temp(:));
            norMap(:,:,patt,rev) = temp/norFac;
        end
    end
    
    thresh = -1e2;
    logNorMap = log(norMap);
    logNorMap(logNorMap < thresh) = thresh;
    
end
