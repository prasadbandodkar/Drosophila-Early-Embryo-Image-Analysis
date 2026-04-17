function d = data(opts)

    arguments
        opts.filename          (1,1) string = ""
        opts.H                 (1,1) double = 0    % image width  (getSizeX)
        opts.W                 (1,1) double = 0    % image height (getSizeY)
        opts.nC                (1,1) double = 0    % active channel count
        opts.Z                 (1,1) double = 0    % z slices
        opts.T                 (1,1) double = 0    % time points
        opts.npts              (1,1) double = 0    % perimeter sample points
        opts.depthInEmbryo     (1,1) double = 0
        opts.depthInImage      (1,1) double = 0
        opts.phi               (1,1) double = 0    % AP axis rotation angle (deg)
        opts.ArcLength         (1,:) double = []   % perimeter length per z
        opts.iPoles            double       = []   % pole locations per z
        opts.radiusofcurvature double       = []   % curvature per z
        opts.genes             (1,1) struct = struct()
        opts.metadata          (1,1) struct = struct()
    end

    d = opts;

end
