function out_variance = block_variance(in_array,blocksize)

[rows,cols] = size(in_array);

out_rows = floor(rows/blocksize);
out_cols = floor(cols/blocksize);

in_array = in_array(1:out_rows*blocksize,1:out_cols*blocksize);

reshaped_array = reshape(in_array,[blocksize out_rows blocksize out_cols]);

out_variance = squeeze(var(reshaped_array,[],[1 3],'omitnan'));