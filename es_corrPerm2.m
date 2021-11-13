function [pval, crit, crit_max, r_obs, r_null] = es_corrPerm(dataX,dataY,nperms,tail,thresh)

defaultStream=RandStream.getGlobalStream;

% nperms = 1000;
% tail = 1;
% thresh = .05;
% dataX = rand(20,100);
% dataY = rand(20,100);

nobs = size(dataX,1);
nvariables = size(dataX,2);

perms=repmat(1:nobs,nperms,1);
perms=perms(reshape(randperm(nperms*nobs),nperms,nobs));

r_obs(1,:) = dataX(:,:)'*dataY;

for p=1:nperms
    r_null(p,:) = dataX(perms(p,:),:)'*dataY;
end

if tail == 1
    maskpos = r_null>repmat(r_obs,nperms,1);
    pval = nansum(maskpos,1)/nperms;
    crit = prctile(r_null,100-round(thresh*100));
    crit_max = max(crit);
elseif tail == -1
    maskneg = r_null<repmat(r_obs,nperms,1);
    pval = nansum(maskneg,1)/nperms;
    crit = prctile(r_null,round(thresh*100));
    crit_max = min(crit);
end

fprintf('\nPermutation test complete!\n');     

