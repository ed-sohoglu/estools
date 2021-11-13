function output = es_reshape(input,dim)
% Reshape vector into multidimensional matrix
% input = data (1 X n), dim(1),dim(2) etc. = dimensions of target matrix

if ~exist('dim','var')
    error('Need to specify dim argument');
end
    
if length(dim) < 5
    dim = [dim repmat(1,1,5-length(dim))];
end

count = 0;
for a=1:dim(1)
    for b=1:dim(2)
        for c=1:dim(3)
            for d=1:dim(4)
                for e=1:dim(5)
                    count = count + 1;
                    output(a,b,c,d,e) = input(count);    
                end
            end
        end
    end
end

