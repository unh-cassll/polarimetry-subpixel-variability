
function [out_I0,out_I45,out_I90,out_I135] = compute_sparse_intensities(I0_stack,I45_stack,I90_stack,I135_stack,kernel_type)

[rows,cols,nframes] = size(I0_stack);

switch kernel_type

    % 2x2 with bilinear interpolation; Polaris Sensor Technologies
    case '2x2'

        num_blocks_rows = floor(rows / 2);
        num_blocks_cols = floor(cols / 2);
        I0_stack = I0_stack(1:2*num_blocks_rows,1:2*num_blocks_cols,:);
        I45_stack = I45_stack(1:2*num_blocks_rows,1:2*num_blocks_cols,:);
        I90_stack = I90_stack(1:2*num_blocks_rows,1:2*num_blocks_cols,:);
        I135_stack = I135_stack(1:2*num_blocks_rows,1:2*num_blocks_cols,:);

        % identify indices for each intensity within a superpixel
        inds0 = [0 0; 0 1];
        inds45 = [0 1; 0 0];
        inds90 = [1 0; 0 0];
        inds135 = [0 0; 1 0];

        if nframes > 1

            % repeat 2x2 matrices to full frame size
            inds0 = logical(repmat(inds0,[[rows cols]/2 nframes]));
            inds45 = logical(repmat(inds45,[[rows cols]/2 nframes]));
            inds90 = logical(repmat(inds90,[[rows cols]/2 nframes]));
            inds135 = logical(repmat(inds135,[[rows cols]/2 nframes]));

            % grab individual polarized light intensities from each frame
            I0_frame = reshape(I0_stack(inds0),[[rows cols]/2 nframes]);
            I45_frame = reshape(I45_stack(inds45),[[rows cols]/2 nframes]);
            I90_frame = reshape(I90_stack(inds90),[[rows cols]/2 nframes]);
            I135_frame = reshape(I135_stack(inds135),[[rows cols]/2 nframes]);

            % bilinear interpolation to full frame size
            I0_frame=imresize3(I0_frame,[rows cols nframes],'linear');
            I45_frame=imresize3(I45_frame,[rows cols nframes],'linear');
            I90_frame=imresize3(I90_frame,[rows cols nframes],'linear');
            I135_frame=imresize3(I135_frame,[rows cols nframes],'linear');

        else

            % repeat 2x2 matrices to full frame size
            inds0 = logical(repmat(inds0,[rows cols]/2));
            inds45 = logical(repmat(inds45,[rows cols]/2));
            inds90 = logical(repmat(inds90,[rows cols]/2));
            inds135 = logical(repmat(inds135,[rows cols]/2));

            % grab individual polarized light intensities from each frame
            I0_frame = reshape(I0_stack(inds0),[rows cols]/2);
            I45_frame = reshape(I45_stack(inds45),[rows cols]/2);
            I90_frame = reshape(I90_stack(inds90),[rows cols]/2);
            I135_frame = reshape(I135_stack(inds135),[rows cols]/2);

            % bilinear interpolation to full frame size
            I0_frame=imresize(I0_frame,[rows cols],'bilinear');
            I45_frame=imresize(I45_frame,[rows cols],'bilinear');
            I90_frame=imresize(I90_frame,[rows cols],'bilinear');
            I135_frame=imresize(I135_frame,[rows cols],'bilinear');

        end

        % Register - measurements to initial positions
        % shift 45 right 1
        I45_frame(:,cols,:) = []; % remove last column
        I45_frame=horzcat(I45_frame(:,1,:), I45_frame); % duplicate first col
        % shift 135 down 1
        I135_frame(rows,:,:) = []; % remove last row
        I135_frame=vertcat(I135_frame(1,:,:), I135_frame); % dupicate first row
        % shift 0 right 1 and down 1
        I0_frame(:,cols,:) = []; % remove last col
        I0_frame=horzcat(I0_frame(:,1,:),I0_frame); % duplicate first col
        I0_frame(rows,:,:) = []; % remove last row
        I0_frame=vertcat(I0_frame(1,:,:),I0_frame); % duplicate first row

        I0_frame = I0_frame(2:2:end,2:2:end,:);
        I45_frame = I45_frame(2:2:end,2:2:end,:);
        I90_frame = I90_frame(2:2:end,2:2:end,:);
        I135_frame = I135_frame(2:2:end,2:2:end,:);

        % [I0_frame] = block_average_and_subsample(I0_frame,2);
        % [I45_frame] = block_average_and_subsample(I45_frame,2);
        % [I90_frame] = block_average_and_subsample(I90_frame,2);
        % [I135_frame] = block_average_and_subsample(I135_frame,2);

        % Ratliff et al. (2009)
        % using 12-pixel kernel approximation, G. Duvarci 2024
        % ON BACK BURNER UNTIL CONVOLUTION APPROACH IS FINALIZED
    case '4x4'

        num_blocks_rows = floor(rows / 4);
        num_blocks_cols = floor(cols / 4);
        I0_stack = I0_stack(1:4*num_blocks_rows,1:4*num_blocks_cols);
        I45_stack = I45_stack(1:4*num_blocks_rows,1:4*num_blocks_cols);
        I90_stack = I90_stack(1:4*num_blocks_rows,1:4*num_blocks_cols);
        I135_stack = I135_stack(1:4*num_blocks_rows,1:4*num_blocks_cols);

        % KERNELS
        A = 0.4086;
        B = 0.2957;
        H0 = [0 0 0 0; 0 A 0 B; 0 0 0 0; 0 B 0 0]; %inds0
        H45 = [0 B 0 0; 0 0 0 0; 0 A 0 B; 0 0 0 0]; %inds45
        H90 = [0 0 B 0; 0 0 0 0; B 0 A 0; 0 0 0 0]; %inds90
        H135 = [0 0 0 0; B 0 A 0; 0 0 0 0; 0 0 B 0]; %inds135

        H0 = repmat(H0,[1 1 num_blocks_rows num_blocks_cols]);
        H45 = repmat(H45,[1 1 num_blocks_rows num_blocks_cols]);
        H90 = repmat(H90,[1 1 num_blocks_rows num_blocks_cols]);
        H135 = repmat(H135,[1 1 num_blocks_rows num_blocks_cols]);

        im0 = reshape(I0_stack,[4 num_blocks_rows 4 num_blocks_cols]);
        im0 = permute(im0,[1 3 2 4]);
        I0_frame = squeeze(sum(double(im0).*H0,[1 2]));

        im45 = reshape(I45_stack,[4 num_blocks_rows 4 num_blocks_cols]);
        im45 = permute(im45,[1 3 2 4]);
        I45_frame = squeeze(sum(double(im45).*H45,[1 2]));

        im90 = reshape(I90_stack,[4 num_blocks_rows 4 num_blocks_cols]);
        im90 = permute(im90,[1 3 2 4]);
        I90_frame = squeeze(sum(double(im90).*H90,[1 2]));

        im135 = reshape(I135_stack,[4 num_blocks_rows 4 num_blocks_cols]);
        im135 = permute(im135,[1 3 2 4]);
        I135_frame = squeeze(sum(double(im135).*H135,[1 2]));

end

out_I0 = I0_frame;
out_I45 = I45_frame;
out_I90 = I90_frame;
out_I135 = I135_frame;
