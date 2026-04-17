function data = analyze(data, cfg)

    arguments
        data (1,1) struct
        cfg  struct
    end

    filename = data.filename;

    % Read file using bfmatlab
    try
        r       = bfGetReader(char(filename));
        omeMeta = r.getMetadataStore();
    catch ME
        error('Error reading file: %s', ME.message);
    end

    % Extract basic metadata
    H  = r.getSizeY();
    W  = r.getSizeX();
    C  = r.getSizeC();
    Z  = r.getSizeZ();
    T  = r.getSizeT();
    scalings(1) = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER));  % um per pixel
    scalings(2) = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER));
    scalings(3) = double(omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER));
    dtype       = char(loci.formats.FormatTools.getPixelTypeString(r.getPixelType()));
    if cfg.debug
        fprintf('Image dimensions: Height: %d x Width: %d x Channels: %d x Slices: %d\n', H, W, C, Z);
        fprintf('Scalings: %g x %g x %g\n', scalings(1), scalings(2), scalings(3));
        fprintf('Image dtype: %s\n', dtype);
    end

    % Extract configuration
    depthInImage = round(cfg.depthInEmbryo / scalings(1));
    fprintf('Depth in image: %d\n', depthInImage);
    npts         = cfg.npts;

    % Default genes if data struct is unpopulated
    if isempty(fieldnames(data.genes))
        if cfg.debug
            fprintf('Data struct is unpopulated. Creating default genes.\n');
        end
        cidx_def = find(~cfg.skip_channels);
        for c = 1:numel(cidx_def)
            gname = char(cfg.channel_names(cidx_def(c)));
            data.genes.(gname) = structs.gene( ...
                name=cfg.channel_names(cidx_def(c)),         ...
                channel_type=cfg.channel_types(cidx_def(c)), ...
                channel_id=int32(cidx_def(c))                ...
            );
        end
    end

    % Active image channel indices from gene structs
    gene_names = fieldnames(data.genes);
    nC         = numel(gene_names);
    cidx       = cellfun(@(g) double(data.genes.(g).channel_id), gene_names);

    % Background estimation for all channels (used for border finding)
    bg_all = zeros(1, C);
    for c_all = 1:C
        I0              = bfGetPlane(r, r.getIndex(0, c_all-1, 0) + 1);
        [bg_all(c_all), ~] = estimate_background(I0, dtype);
    end
    bg = bg_all(cidx);   % background for requested channels only

    % Pre-allocate gene data arrays
    for c = 1:nC
        gname                  = gene_names{c};
        data.genes.(gname).t   = zeros(npts+1, Z);
        data.genes.(gname).raw = zeros(npts+1, Z);
        data.genes.(gname).s   = zeros(npts+1, 1);
    end

    % Per-z loop
    Imax      = zeros(H, W, 'like', bfGetPlane(r, 1));
    Xp        = zeros(npts+1, Z);
    Yp        = zeros(npts+1, Z);
    ArcLength = zeros(1, Z);
    for z = 1:Z

        if cfg.debug
            fprintf('Processing slice %d/%d\n', z, Z);
        end

        % Step 1: Load all channels and find borders
        Iz_all = zeros(H, W, C, 'like', bfGetPlane(r, 1));
        for c_all = 1:C
            I               = bfGetPlane(r, r.getIndex(z-1, c_all-1, 0) + 1);
            Iz_all(:,:,c_all) = imsubtract(I, bg_all(c_all));
        end
        Iz_max       = max(Iz_all, [], 3);
        Imax         = max(Imax, Iz_max);
        [xp, yp]     = find_borders(Iz_max, npts);
        Xp(:,z)      = xp';
        Yp(:,z)      = yp';
        ArcLength(z) = sum(sqrt(diff(xp).^2 + diff(yp).^2));

        % Step 2: Measure gene expression for requested channels only
        Iz                = Iz_all(:,:,cidx);
        cfg.plot_filename = sprintf('slice%d.png', z);
        cfg.plot_title    = sprintf('Domain Measure | Slice %d', z);
        [t, raw, s]       = domain_measure(Iz, xp, yp, depthInImage, npts, cfg);
        for c = 1:nC
            gname = gene_names{c};
            data.genes.(gname).t(:,z)   = t(:,c);
            data.genes.(gname).raw(:,z) = raw(:,c);
        end

    end

    % s is constant across z and genes — store once
    for c = 1:nC
        data.genes.(gene_names{c}).s = s;
    end

    % AP axis angle from max-projection
    [xp, yp] = find_borders(Imax, npts);
    phi = functions.ellipse_fit(xp, yp);
    if phi > pi/2, phi = phi - pi; end
    phi = rad2deg(phi);

    % Assemble metadata
    meta.scalings = scalings;
    meta.dtype    = dtype;
    meta.rho      = round(scalings(3) / scalings(1));
    meta.Xp       = Xp;
    meta.Yp       = Yp;
    meta.bg       = bg_all;   % background per file channel (1 x C), indexed by channel_id

    % Update data struct fields
    data.H                 = H;
    data.W                 = W;
    data.nC                = nC;
    data.Z                 = Z;
    data.T                 = T;
    data.npts              = npts;
    data.depthInEmbryo     = cfg.depthInEmbryo;
    data.depthInImage      = depthInImage;
    data.phi               = phi;
    data.ArcLength = ArcLength;
    data.metadata  = meta;

    r.close();

end
