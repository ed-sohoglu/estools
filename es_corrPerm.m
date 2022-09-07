function [pval, r_obs, r_null, crit] = es_corrPerm(dataX,dataY,nperms,tail,thresh)

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

for v=1:nvariables
    
    if v == 1
        fprintf('\n%d Variables...\n',nvariables);     
        fprintf('\nRunning permutation test for variable %d...\n',v);
    else
        if ~rem(v-100,1000); fprintf('\nRunning permutation test for variable %d...\n',v); end
    end
    
    r_obs(1,v) = fisherTransform(corr(dataX(:,v),dataY(:,v),'Type','Spearman','Rows','Pairwise'));
    
    for p=1:nperms
                
        r_null(p,v) = fisherTransform(corr(dataX(perms(p,:),v),dataY(:,v),'Type','Spearman','Rows','Pairwise'));
        
    end
    
    if tail==1
        pval(v) = length(find(r_null(:,v)>r_obs(1,v)))/nperms;
        crit(v) = prctile(r_null(:,v),100-round(thresh*100));
    elseif tail==-1
        pval(v) = length(find(r_null(:,v)<r_obs(1,v)))/nperms;
        crit(v) = prctile(r_null(:,v),round(thresh*100));
    end
        
    
end

fprintf('\nPermutation test complete!\n');     

