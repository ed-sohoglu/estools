function output = es_avgFactors(input,dim,factors)
% Average over levels of certain factors in a factorial design
% input = data(subjects X conditions), dim(1),dim(2) etc. = dimensions of
% factorial design, factors = indices of factors over which you want to
% average

if ~exist('factors','var')
    error('Need to specify factors argument');
end

if ~exist('dim','var')
    error('Need to specify dim argument');
end

if ndims(dim) < 5
    dim = [dim repmat(1,1,5-length(dim))];
end

count = 0;
for a=1:dim(1)
    for b=1:dim(2)
        for c=1:dim(3)
            for d=1:dim(4)
                for e=1:dim(5)
                    count = count + 1;
                    tmp(:,a,b,c,d,e) = input(:,count);    
                end
            end
        end
    end
end

for f=1:length(factors)
    tmp = mean(tmp,factors(f)+1);
end

tmp = squeeze(tmp);

dim = size(tmp);
dim = dim(2:end);

if ndims(dim) < 5
    dim = [dim repmat(1,1,5-length(dim))];
end

count = 0;
for a=1:dim(1)
    for b=1:dim(2)
        for c=1:dim(3)
            for d=1:dim(4)
                for e=1:dim(5)
                    count = count + 1;
                    output(:,count) = tmp(:,a,b,c,d,e);    
                end
            end
        end
    end
end

