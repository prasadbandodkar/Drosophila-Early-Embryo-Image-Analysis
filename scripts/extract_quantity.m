% extract_quantity.m
% Extract a single quantity from all processed .mat files in a data folder.
%
% Usage: set the three parameters below, then run.
%
% Examples:
%   gene_name = 'sna';  field = 'Width';        % scalar per file
%   gene_name = 'sna';  field = 'LeftBorder';   % scalar per file
%   gene_name = 'sna';  field = 'sPeaks';       % array per file

% Usage: set the three parameters below, then run.
data_dir  = './data';   % root folder containing per-file subfolders
gene_name = 'sna';       % gene name (must match field in data.genes)
field     = 'Width';     % field inside the gene struct

out_dir   = './';                                       % directory for the output CSV
out_file  = sprintf('%s_%s.csv', gene_name, field);      % output CSV filename

% -------------------------------------------------------------------------

mat_files = dir(fullfile(data_dir, '**', '*.mat'));
fprintf('Found %d .mat files in %s\n', numel(mat_files), data_dir);

lsm_names  = {};
mat_paths  = {};
values     = {};

for i = 1:numel(mat_files)
    mat_path = fullfile(mat_files(i).folder, mat_files(i).name);
    try
        tmp = load(mat_path, 'data');
    catch
        fprintf('Could not load %s — skipping\n', mat_path);
        continue
    end

    if ~isfield(tmp, 'data')
        fprintf('No ''data'' variable in %s — skipping\n', mat_path);
        continue
    end
    d = tmp.data;

    if ~isfield(d, 'genes') || ~isfield(d.genes, gene_name)
        fprintf('Gene ''%s'' not found in %s — skipping\n', gene_name, mat_path);
        continue
    end

    gene = d.genes.(gene_name);

    if ~isfield(gene, field)
        fprintf('Field ''%s'' not found in gene ''%s'' in %s — skipping\n', field, gene_name, mat_path);
        continue
    end

    val = gene.(field);
    lsm_names{end+1} = char(d.filename); %#ok<SAGROW>
    mat_paths{end+1} = mat_path;         %#ok<SAGROW>
    values{end+1}    = val;              %#ok<SAGROW>
    fprintf('[%d/%d] %s → %s\n', i, numel(mat_files), mat_files(i).name, mat2str(val));
end

fprintf('\nExtracted ''%s.%s'' from %d/%d files.\n', gene_name, field, numel(lsm_names), numel(mat_files));

% Sort by mean value (increasing)
means_for_sort = cellfun(@(v) mean(v(:), 'omitnan'), values);
[~, sort_idx] = sort(means_for_sort, 'ascend');
lsm_names = lsm_names(sort_idx);
mat_paths = mat_paths(sort_idx);
values    = values(sort_idx);

% Write CSV
csv_path = fullfile(out_dir, out_file);
fid = fopen(csv_path, 'w');
fprintf(fid, 'filename,mat_path,%s,%s_mean\n', field, field);
for i = 1:numel(lsm_names)
    val      = values{i};
    val_mean = mean(val(:), 'omitnan');
    if isscalar(val)
        fprintf(fid, '%s,%s,%.6g,%.6g\n', lsm_names{i}, mat_paths{i}, val, val_mean);
    else
        % Array: write as quoted space-separated values
        val_str = strjoin(arrayfun(@(v) sprintf('%.6g', v), val(:)', 'UniformOutput', false), ' ');
        fprintf(fid, '%s,%s,"%s",%.6g\n', lsm_names{i}, mat_paths{i}, val_str, val_mean);
    end
end
fclose(fid);
fprintf('Saved to %s\n', csv_path);

% Summary
means = means_for_sort(sort_idx);
fprintf('Mean across files: %.4f  Std: %.4f  N: %d\n', mean(means), std(means), numel(means));
