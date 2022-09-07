function z = es_zscoremat(data)
% Z-score over ALL dimensions of data matrix (different to matlab's zscore() which is
% peformed in vectorized fashion i.e. over specified dimensions like rows or columns

z = (data-mean(data(:)))./std(data(:));