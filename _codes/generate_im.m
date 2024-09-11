% Takes in array of double, outputs as scaled unsigned 8-bit integer
% N. Laxague, 2019
function out_im = generate_im(in_array,cmap,low_val,high_val)

% size of array
[s1,s2] = size(in_array);

% scale colormap
cmap = cmap/max(max(cmap));
cmap_8bit = cmap*256;

% scale array
scaled_array = in_array - low_val;
scaled_array = scaled_array/(high_val-low_val)*256;
scaled_array(scaled_array<1) = 1;
scaled_array(scaled_array>256) = 256;
scaled_array = floor(scaled_array);

% compute R, G, & B channel values
scaled_vec = reshape(scaled_array,s1*s2,1);
out_R = reshape(cmap_8bit(scaled_vec,1),s1,s2);
out_G = reshape(cmap_8bit(scaled_vec,2),s1,s2);
out_B = reshape(cmap_8bit(scaled_vec,3),s1,s2);

% build RGB array
out_RGB = ones(s1,s2,3);
out_RGB(:,:,1) = out_R;
out_RGB(:,:,2) = out_G;
out_RGB(:,:,3) = out_B;
out_im = uint8(out_RGB);