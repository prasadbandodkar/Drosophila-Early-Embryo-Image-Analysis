function cfg = config(opts)

arguments
    opts.path_raw           (1,1) string  = '/Volumes/X1/Archive/Projects/dgrp/Data_lsm'
    opts.filetype           (1,1) string  = 'lsm'
    opts.subfolders         (1,1) logical = true
    opts.path_data          (1,1) string  = './temp'
    opts.overwrite_data     (1,1) logical = true
    opts.substring          (1,1) string  = 'BcdKrEve'


    opts.channel_names      (1,:) string  = ["Eve", "Kr", "DIC", "Dapi", "Bcd"]
    opts.channel_types      (1,:) string  = ["mRNA", "mRNA", "DIC", "Dapi", "mRNA"]
    opts.skip_channels      (1,:) logical = [false, false, true, true, true]
    opts.data_type          (1,1) string  = 'dgrp'
    opts.image_section      (1,1) string  = 'sagittal'
    
    
    opts.run_parallel       (1,1) logical = false
    opts.num_workers        (1,1) double  = 4
    opts.yesplot            (1,1) logical = true
    opts.depthInEmbryo      (1,1) double  = 30
    opts.npts               (1,1) double  = 2000
    opts.bkgrndthresh       (1,1) double  = 0.15
    opts.peakthresh         (1,1) double  = 0.1
    opts.debug              (1,1) logical = true
    opts.info               (1,1) logical = true
    opts.image_type         (1,1) string  = 'cross-section'     % 'cross-section' or 'sagittal'
end

if opts.debug
    opts.info = true;
end

if isempty(opts.path_data)
    opts.path_data = opts.path_raw;
end

cfg = opts;

end
