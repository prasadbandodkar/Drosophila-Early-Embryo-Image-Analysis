function data = fitting(data, cfg)

arguments
    data (1,1) struct
    cfg  struct
end


bkgrndthresh = cfg.bkgrndthresh;
peakthresh   = cfg.peakthresh;

% Iterate over non-skipped genes
gene_names = cellstr(cfg.channel_names(~cfg.skip_channels));
nG         = numel(gene_names);

% Compute pole indices once across all slices, and find the ventral
% midline at each z using all channel profiles simultaneously.
Xp      = data.metadata.Xp;
Yp      = data.metadata.Yp;
Z       = size(Xp, 2);
IPoles  = zeros(Z, 2);

% Gather all-channel T and S0 for find_midline
npts   = size(data.genes.(gene_names{1}).t, 1) - 1;
s      = linspace(-1, 1, npts + 1)';
T_all  = zeros(npts + 1, Z, nG);
S0_all  = cell(nG, 1);
rbr_all = zeros(nG, 1);
for g = 1:nG
    gn           = gene_names{g};
    T_all(:,:,g) = data.genes.(gn).t;
    [sPeaks_g, ~, rbr_g] = get_gene_info(gn, cfg);
    S0_all{g}    = sPeaks_g(:);   % canonical peak locations for this channel
    rbr_all(g)   = rbr_g;
end

% Compute pole indices once across all slices
% If sagittal, then use the curvature of the embryo to find the poles
% If cross-section, then use the gene expression profiles to find the midline
for z = 1:Z
    if strcmp(cfg.image_section, 'sagittal')
        [ant, post]  = find_poles(Xp(:,z), Yp(:,z));
        IPoles(z, :) = [ant.idx, post.idx];
    elseif strcmp(cfg.image_section, 'cross-section')
        % Background-subtract each channel before midline detection
        T_z = zeros(npts + 1, nG);
        for g = 1:nG
            T_z(:,g) = subtract_background(T_all(:, z, g), rbr_all(g));
        end
        [ven_z, dor_z] = find_midline(s, T_z, S0_all, nan(1, nG));
        IPoles(z, :)  = [ven_z.idx, dor_z.idx];
    end
end
data.IPoles = IPoles;


% Compute fit data for each gene
for g = 1:nG
    if cfg.info
        fprintf('Fitting gene %d/%d: %s\n', g, nG, gene_names{g});
    end

    gname = gene_names{g};
    T     = data.genes.(gname).t;   % (npts+1) x Z

    % Get prior peak positions, peak count, and background ratio for this gene
    [sPeaks, nPeaks, rbr] = get_gene_info(gname, cfg);

    % Fit the gene
    [fitdata, meta1, meta2] = fit_channel(T, rbr, gname, nPeaks, sPeaks, IPoles, cfg);

    % Store fit results back into gene struct
    data.genes.(gname).sPeaks    = fitdata.sPeaks;
    data.genes.(gname).LeftBorder = fitdata.LeftBorder;
    data.genes.(gname).RightBorder = fitdata.RightBorder;
    data.genes.(gname).Width     = fitdata.Width;

    % Per-side fit metadata
    fn = fieldnames(meta1.side1);
    for j = 1:numel(fn)
        data.genes.(gname).meta.side1.(fn{j}) = meta1.side1.(fn{j});
        data.genes.(gname).meta.side2.(fn{j}) = meta1.side2.(fn{j});
    end

    % Global fit metadata
    fn = fieldnames(meta2.side1);
    for j = 1:numel(fn)
        data.metadata.side1.(fn{j}) = meta2.side1.(fn{j});
        data.metadata.side2.(fn{j}) = meta2.side2.(fn{j});
    end
end



if cfg.debug
    r_dbg = bfGetReader(char(data.filename));
    nPts  = size(data.metadata.Xp, 1) - 1;   % number of perimeter points (excl. wrap)

    for g = 1:nG
        gname = gene_names{g};
        cid   = double(data.genes.(gname).channel_id);
        LB    = data.genes.(gname).LeftBorder;    % (Z x nPeaks)
        RB    = data.genes.(gname).RightBorder;
        SP    = data.genes.(gname).sPeaks;
        nPk   = size(LB, 2);

        for z = 1:Z
            % Load this gene's channel image for slice z
            I = bfGetPlane(r_dbg, r_dbg.getIndex(z-1, cid-1, 0) + 1);
            I = imsubtract(I, data.metadata.bg(cid));

            xp    = data.metadata.Xp(:, z);   % (nPts+1) x 1, circular (first==last)
            yp    = data.metadata.Yp(:, z);
            apole = round(IPoles(z, 1));
            ppole = round(IPoles(z, 2));

            % Reconstruct the perimeter index arrays for each side,
            % mirroring get_side's index arithmetic (operates on 1:nPts).
            if apole < ppole
                idx1 = (apole+1 : ppole)';
                idx2 = [(apole-1:-1:1), (nPts:-1:ppole)]';
            else
                idx1 = [(apole+1:nPts), (1:ppole)]';
                idx2 = (apole-1:-1:ppole)';
            end

            % Interpolate s∈[0,1] → pixel (x,y) along a side's index array
            s_to_x = @(sv, idx) interp1(linspace(0,1,numel(idx))', xp(idx), sv, 'linear', 'extrap');
            s_to_y = @(sv, idx) interp1(linspace(0,1,numel(idx))', yp(idx), sv, 'linear', 'extrap');

            fig = figure('Visible', 'off');
            imagesc(I); colormap(gray); axis image; hold on;

            % Embryo contour
            plot([xp; xp(1)], [yp; yp(1)], 'w-', 'LineWidth', 0.8);

            % Pole markers (triangle up = pole 1, triangle down = pole 2)
            plot(xp(apole), yp(apole), 'y^', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
            plot(xp(ppole), yp(ppole), 'yv', 'MarkerSize', 10, 'MarkerFaceColor', 'y');

            % Left border (green square), right border (blue square), peak (red circle)
            % drawn on both sides of the perimeter
            for pk = 1:nPk
                for sd = 1:2
                    if sd == 1
                        idx = idx1;
                    else
                        idx = idx2;
                    end
                    if isempty(idx), continue; end

                    xlb = s_to_x(LB(z,pk), idx);  ylb = s_to_y(LB(z,pk), idx);
                    xrb = s_to_x(RB(z,pk), idx);  yrb = s_to_y(RB(z,pk), idx);
                    xsp = s_to_x(SP(z,pk), idx);  ysp = s_to_y(SP(z,pk), idx);

                    plot(xlb, ylb, 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
                    plot(xrb, yrb, 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
                    plot(xsp, ysp, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'y');

                    text(xlb, ylb, sprintf('  LB%d: %.3f', pk, LB(z,pk)), ...
                        'Color', 'g', 'FontSize', 7, 'VerticalAlignment', 'middle');
                    text(xrb, yrb, sprintf('  RB%d: %.3f', pk, RB(z,pk)), ...
                        'Color', 'r', 'FontSize', 7, 'VerticalAlignment', 'middle');
                    text(xsp, ysp, sprintf('  P%d: %.3f',  pk, SP(z,pk)), ...
                        'Color', 'y', 'FontSize', 7, 'VerticalAlignment', 'middle');
                end
            end

            title(sprintf('%s | Z = %d', gname, z));
            fname = sprintf('debug_%s_z%02d.png', gname, z);
            saveas(fig, fullfile(cfg.path_data, fname));
            close(fig);
        end
    end

    r_dbg.close();
end


end


% -------------------------------------------------------------------------
function [sPeaks, nPeaks, rbr] = get_gene_info(gname, cfg)
% Load prior peak positions and background ratio for a gene from mats/.

    files  = os.get_files("./src/+geneaverages", "mat", true, gname);
    nPeaks = 0;
    rbr    = 0;
    sPeaks = [];
    for i = 1:numel(files)
        if cfg.debug
            fprintf('\tLoading gene info for %s from %s\n', gname, files(i));
        end
        d           = load(files(i));
        sPeaks(i,:) = d.s_peak;
        nPeaks      = nPeaks + 1;
        rbr         = rbr + (d.rbr - rbr) / i;   % running mean
    end
end


% -------------------------------------------------------------------------
function [fitdata, meta1, meta2] = fit_channel(T, rbr, gname, nPeaks, sPeaks, IPoles, cfg)
% Fit expression peaks for each z-slice of a gene.

nSlices = size(T, 2);

for i = 1:nSlices
    t      = subtract_background(T(:,i), rbr);
    iPoles = IPoles(i,:);

    if cfg.debug
        fprintf('\tFitting slice %d, gene %s\n', i, gname);
    end

    % Side 1: anterior → posterior or dorsal → ventral
    t1 = get_side(t, iPoles, 1);
    cfg.plot_filename = sprintf('%s_slice%d_side1_peaks.png', gname, i);
    cfg.plot_title    = sprintf('%s | Slice %d | Side 1', gname, i);
    [H, I, bcoor] = find_peaks_guided(t1, nPeaks, sPeaks, cfg);
    [data1, metadata1] = fit_peaks(t1, gname, nPeaks, H, I, bcoor, cfg);

    % Side 2: posterior → anterior or ventral → dorsal
    t2 = get_side(t, iPoles, 2);
    cfg.plot_filename = sprintf('%s_slice%d_side2_peaks.png', gname, i);
    cfg.plot_title    = sprintf('%s | Slice %d | Side 2', gname, i);
    [H, I, bcoor] = find_peaks_guided(t2, nPeaks, sPeaks, cfg);
    [data2, metadata2] = fit_peaks(t2, gname, nPeaks, H, I, bcoor, cfg);

    % Combine both sides by averaging
    fitdata.sPeaks(i,:)    = (data1.sPeaks    + data2.sPeaks)    / 2;
    fitdata.LeftBorder(i,:) = (data1.LeftBorder + data2.LeftBorder) / 2;
    fitdata.RightBorder(i,:) = (data1.RightBorder + data2.RightBorder) / 2;
    fitdata.Width(i,:)     = (data1.Width     + data2.Width)     / 2;

    % Per-side raw fit results
    fn = fieldnames(data1);
    for j = 1:numel(fn)
        meta1.side1.(fn{j}) = data1.(fn{j});
        meta1.side2.(fn{j}) = data2.(fn{j});
    end

    % Per-side fit metadata
    fn = fieldnames(metadata1);
    for j = 1:numel(fn)
        meta2.side1.(fn{j}) = metadata1.(fn{j});
        meta2.side2.(fn{j}) = metadata2.(fn{j});
    end
end

end


% -------------------------------------------------------------------------
function t_out = subtract_background(t, rbr)
% Rolling-ball background subtraction on a circular 1-D profile.

t    = t(:);              % ensure column vector
npts = numel(t) - 1;
rbr  = rbr * npts;        % "rolling ball radius"

t1    = [t(1:end-1); t; t(2:end)];
se    = strel('line', rbr, 90);
Itop  = imtophat(t1, se);
t_out = Itop(npts:2*npts);

end
