function cfg = stagingconfig(opts)

arguments
    opts.path_raw           (1,1) string  = '/Volumes/X2/Projects/phd/staging/data/fixed_raw/Carrell_Sophie/Confocal Images/2012-06-22 cact_and_yw/yw 1'
    opts.filetype           (1,1) string  = 'lsm'
    opts.subfolders         (1,1) logical = false
    opts.path_data          (1,1) string  = './data'
    opts.overwrite_data     (1,1) logical = true
    opts.save_analysis      (1,1) logical = true
    opts.substring          (1,1) string  = ''


    opts.channel_names      (1,:) string  = ["Dapi", "Dl", "DIC", "sna"]
    opts.channel_types      (1,:) string  = ["Dapi", "nuclear protein", "DIC", "mRNA"]
    opts.skip_channels      (1,:) logical = [true, true, true, false]
    opts.data_type          (1,1) string  = 'staging'
    opts.image_section      (1,1) string  = 'cross-section'
    
    
    opts.run_parallel       (1,1) logical = false
    opts.num_workers        (1,1) double  = 4
    opts.yesplot            (1,1) logical = true
    opts.depthInEmbryo      (1,1) double  = 40
    opts.npts               (1,1) double  = 200
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
