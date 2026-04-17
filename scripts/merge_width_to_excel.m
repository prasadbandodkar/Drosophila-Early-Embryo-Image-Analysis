% merge_width_to_excel.m
% Match filenames in sna_Width.csv with rows in an Excel sheet and
% write the mean width as a new column.

csv_path       = './sna_Width.csv';
excel_path     = './files/staging_sna_fixed 2.xlsx';
excel_sheet    = 1;              % sheet index or name
filename_col   = 1;             % column index in Excel that contains filenames
out_col_header = 'Width_mean';  % header for the new column

% -------------------------------------------------------------------------

% Read the CSV (filename, Width, Width_mean)
csv = readtable(csv_path, 'TextType', 'string');

% Pre-compute stems from CSV filenames for fast lookup
csv_stems = arrayfun(@(f) get_stem(f), csv.filename, 'UniformOutput', false);
csv_stems = string(csv_stems);

% Read Excel as raw cell array (preserves all content)
raw = readcell(excel_path, 'Sheet', excel_sheet);

% Find or append the output column in the header row
header = raw(1, :);
out_col_idx = find(strcmpi(header, out_col_header));
if isempty(out_col_idx)
    out_col_idx = size(raw, 2) + 1;
end
raw{1, out_col_idx} = out_col_header;

% Match each Excel row to a CSV entry by base filename stem
n_matched = 0;
for r = 2:size(raw, 1)
    cell_val = raw{r, filename_col};

    if ismissing(cell_val) || (ischar(cell_val) && isempty(strtrim(cell_val)))
        raw{r, out_col_idx} = missing;
        continue
    end

    excel_stem = get_stem(string(cell_val));
    match_idx  = find(strcmpi(csv_stems, excel_stem), 1);

    if ~isempty(match_idx)
        raw{r, out_col_idx} = csv.Width_mean(match_idx);
        n_matched = n_matched + 1;
    else
        raw{r, out_col_idx} = missing;
    end
end

fprintf('Matched %d/%d Excel rows.\n', n_matched, size(raw, 1) - 1);

% Write updated cell array back to Excel
writecell(raw, excel_path, 'Sheet', excel_sheet);
fprintf('Saved to %s\n', excel_path);

% -------------------------------------------------------------------------
function stem = get_stem(fullpath)
    [~, stem, ~] = fileparts(char(fullpath));
end
