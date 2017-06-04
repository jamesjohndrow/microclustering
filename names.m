clear; 
resample = false;

datf = readtable('Data/names/first.csv');
fndist = datf.freq_f;
fd = datf.first;

datl = readtable('Data/names/last.csv');
lndist = datl.freq_l;
ld = datl.last;

if resample
    ls = randsample(length(lndist),sum(lndist),true,lndist./sum(lndist));
    fs = randsample(length(fndist),sum(lndist),true,fndist./sum(fndist));
    fl = [fs ls];
    clear ls fs;

    [~,~,ic] = unique(fl,'rows');
    sC = accumarray(ic,1);

    save('nmsamp.mat','sC');
else
   load('nmsamp.mat'); 
end
ms = nmmoments(sC);

nratio = sum(sC).^2./(sum(sC.^2));
