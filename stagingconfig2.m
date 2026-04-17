function cfg = stagingconfig2(opts)

arguments
    opts.path_raw           (1,:) string  = ["/Volumes/X2/Projects/phd/staging/data/fixed_raw/Reeves_Greg/Postdoc Images/2011-02-11_yw_sna_dl_H3.mdb", ...
                                             "/Volumes/X2/Projects/phd/staging/data/fixed_raw/Reeves_Greg/Postdoc Images/2011-02-17_yw_sna_dl_H3_more.mdb", ...
                                             "/Volumes/X2/Projects/phd/staging/data/fixed_raw/Reeves_Greg/Postdoc Images/2011-04-01_yw_sna_dl_H3_25C_more.mdb", ...
                                             "/Volumes/X2/Projects/phd/staging/data/fixed_raw/Reeves_Greg/Postdoc Images/2011-07-12_yw_sna_dl_H3_July2011.mdb"]
    opts.filetype           (1,1) string  = 'lsm'
    opts.subfolders         (1,1) logical = false
    opts.path_data          (1,1) string  = './data'
    opts.num_file_parts     (1,1) double  = 4      % number of file parts to use for filename
    % opts.subfolder          (1,1) string  = ''
    opts.overwrite_data     (1,1) logical = false
    opts.save_analysis      (1,1) logical = true
    opts.substring          (1,1) string  = ''


    opts.channel_names      (1,:) string  = ["sna", "Dl", "Dapi"]
    opts.channel_types      (1,:) string  = ["mrna", "nuclear protein", "Dapi"]
    opts.skip_channels      (1,:) logical = [false, false, true]
    opts.data_type          (1,1) string  = 'staging'
    opts.image_section      (1,1) string  = 'cross-section'
    
    
    opts.run_parallel       (1,1) logical = false
    opts.num_workers        (1,1) double  = 4

    opts.yesplot            (1,1) logical = true
    opts.depthInEmbryo      (1,1) double  = 18.36
    opts.npts               (1,1) double  = 300
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
