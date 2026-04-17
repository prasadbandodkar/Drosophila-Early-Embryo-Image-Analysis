function [Img, bwlabel] = segment_embryo_image(I)
% Segment the embryo from the background.
%
% Inputs:
%   I:       Image (2D or 3D, any numeric class)
% Outputs:
%   Img:     Segmented image (background zeroed, same size/class as I)
%   bwlabel: Binary mask (logical), true where embryo is detected
%
% Pipeline (mirrors CVImage._segment_embryo_image in cvimage.py):
%   1.   Max-project to 2D
%   2.   Replicate-pad and downsample to 256 px for speed
%   3.   Normalize intensity to uint8 [0,255]
%   4.   Gaussian smoothing (sigma=3.5, equivalent to 21x21 kernel)
%   4.5. Illumination correction (flat-field + CLAHE)
%   5.   Otsu thresholding
%   6.   Morphological cleaning  (open 3x3, close 21x21)
%   7.   Fill interior           (close 21x21)
%   8.   Keep largest connected component
%   9.   Crop padding, upsample mask to original resolution
%   10.  Apply mask to original image

% =====================================================================
% STEP 1: Collapse to 2D and record original size
% =====================================================================
Itmp           = max(I, [], 3);
[orig_h, orig_w] = size(Itmp);

% =====================================================================
% OPTIMIZATION: Add replicate padding then downsample to 256 px
% =====================================================================
pad       = 44;
Img_pad   = padarray(Itmp, [pad, pad], 'replicate');
[ph, pw]  = size(Img_pad);
seg_size  = 256;
scale     = seg_size / max(ph, pw);
Img_small = imresize(Img_pad, round([ph, pw] * scale));
pad_small = round(pad * scale);

% =====================================================================
% STEP 2: Normalize to uint8 [0, 255]
% =====================================================================
Itmp = im2uint8(mat2gray(Img_small));

% =====================================================================
% STEP 3: Gaussian smoothing (sigma = 3.5, matches 21x21 in Python)
% =====================================================================
Itmp = imgaussfilt(Itmp, 3.5);

% =====================================================================
% STEP 3.5: Illumination correction (flat-field + CLAHE)
% =====================================================================
Itmp = correct_illumination(Itmp);

% =====================================================================
% STEP 4: Otsu thresholding
% =====================================================================
bw = imbinarize(Itmp);

% =====================================================================
% STEP 5: Morphological cleaning — open (3x3) then close (21x21)
% =====================================================================
se_small = strel('disk', 1);   % radius 1  ≈ 3x3 ellipse
se_large = strel('disk', 10);  % radius 10 ≈ 21x21 ellipse
bw = imopen( bw, se_small);
bw = imclose(bw, se_large);

% =====================================================================
% STEP 6: Fill interior — second large closing bridges the nuclear ring
% =====================================================================
bw = imclose(bw, se_large);

% =====================================================================
% STEP 7: Keep largest connected component (8-connectivity)
% =====================================================================
CC = bwconncomp(bw, 8);
if CC.NumObjects > 1
    areas    = cellfun(@numel, CC.PixelIdxList);
    [~, idx] = max(areas);
    bw = false(size(bw));
    bw(CC.PixelIdxList{idx}) = true;
end

% =====================================================================
% STEP 8: Crop replicate-padding then upsample to original resolution
% =====================================================================
if pad_small > 0
    bw = bw(pad_small+1:end-pad_small, pad_small+1:end-pad_small);
end
bwlabel = imresize(uint8(bw), [orig_h, orig_w], 'nearest') > 0;

% =====================================================================
% STEP 9: Apply mask to original image (handles 2D and 3D)
% =====================================================================
mask = cast(bwlabel, class(I));
if ndims(I) == 3
    mask = repmat(mask, 1, 1, size(I, 3));
end
Img = I .* mask;

end


% -----------------------------------------------------------------------
function out = correct_illumination(Img)
% Flat-field correction followed by CLAHE.
% Mirrors CVImage._correct_illumination in cvimage.py.
%
% Stage 1: Large-kernel Gaussian estimates the illumination gradient.
%   Corrected = (Image - dark) / (background - dark) * mean(background)
% Stage 2: CLAHE for local contrast enhancement.

[h, w] = size(Img);
Ifl    = double(Img);

% --- kernel size: half of longest dimension, clamped to [51, 201], odd ---
ksize = floor(max(h, w) / 2);
ksize = max(51, min(ksize, 201));
if mod(ksize, 2) == 0
    ksize = ksize + 1;
end
sigma = ksize / 3.0;

% --- flat-field correction ---
background  = imgaussfilt(Ifl, sigma, 'FilterSize', ksize);
dark_frame  = prctile(Ifl(:), 1);
bg_mean     = mean(background(:));
denominator = max(background - dark_frame, 1.0);
corrected   = ((Ifl - dark_frame) ./ denominator) .* bg_mean;
corrected   = uint8(min(max(corrected, 0), 255));

% --- CLAHE (clipLimit=0.02 ≈ cv2 clipLimit=2.0) ---
out = adapthisteq(corrected, 'ClipLimit', 0.02, 'NumTiles', [4, 4]);

end
