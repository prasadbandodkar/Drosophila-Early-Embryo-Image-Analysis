function g = gene(opts)

    arguments
        opts.name     (1,1) string = ""
        opts.channel_type (1,1) string = 'mRNA'
        opts.channel_id (1,1) int32  = 0
    end

    % Secondary Data
    opts.t   = [];
    opts.raw = [];
    opts.s   = [];

    % Tertiary Data
	switch opts.channel_type
		case 'nuclei'
			opts.channel_type = 'nuclei';
		case 'mRNA'
			opts.channel_type = 'non-nuclear protein or mRNA';
			opts.sPeaks       = [];
            opts.LeftBorder   = [];
            opts.RightBorder  = [];
            opts.Width        = [];
		case 'nuclear protein'
			opts.channel_type = 'nuclear protein';
		case 'intronic probe'
			opts.channel_type = 'intronic probe';
		case 'N/A'
			opts.channel_type = 'N/A';
	end

    % Metadata
    opts.meta = struct();

    g = opts;

end
