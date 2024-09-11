% Given input matrix of arbitrary dimension, perform block averaging and
% subsampling according to input 'blocksize' array
%
% N. Laxague 2024
%
function [out_array] = block_average_and_subsample(in_array,blocksize)

array_size = size(in_array);

array_size_truncated = floor(array_size(:)'./blocksize(:)');

array_truncated = trimdata(in_array,array_size_truncated.*blocksize(:)');

reshape_sizes = reshape([blocksize(:) array_size_truncated(:)]',1,[]);

array_reshaped = reshape(array_truncated,reshape_sizes);

out_array = squeeze(mean(array_reshaped,1:2:length(array_reshaped),'omitnan'));