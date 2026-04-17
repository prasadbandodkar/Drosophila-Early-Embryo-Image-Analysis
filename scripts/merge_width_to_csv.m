% merge_width_to_csv.m
% Join staging_sna_fixed 2.csv with sna_Width.csv and write a new CSV
% with an added "new_width" column (Width_mean from sna_Width.csv).
% Matching is done on the full file path.

staging_csv = './files/staging_sna_fixed 2.csv';
width_csv   = './sna_Width.csv';
out_csv     = './files/staging_sna_fixed_with_new_width.csv';

% -------------------------------------------------------------------------

staging = readtable(staging_csv, 'TextType', 'string', 'Delimiter', ',', 'VariableNamingRule', 'preserve');
width   = readtable(width_csv,   'TextType', 'string', 'Delimiter', ',', 'VariableNamingRule', 'preserve');

staging_files = staging{:, 1};   % filelocation column by index
width_files   = width{:, 1};     % filename column by index

% Left-join: width_csv is the base. For each width row, find the matching
% staging row (if any) and append its columns. Unmatched rows are kept with
% NaN / "" for the staging columns.
staging_varnames = staging.Properties.VariableNames;

% Pre-build an empty staging row (NaN for numerics, "" for strings)
empty_staging = staging(1, :);
for c = 1:size(staging, 2)
    if isnumeric(staging{1, c}) || islogical(staging{1, c})
        empty_staging{1, c} = NaN;
    else
        empty_staging{1, c} = missing;
    end
end

n_width   = height(width);
n_matched = 0;
staging_rows = repmat(empty_staging, n_width, 1);   % pre-allocate

for i = 1:n_width
    idx = find(strcmp(staging_files, width_files(i)), 1);
    if ~isempty(idx)
        staging_rows(i, :) = staging(idx, :);
        n_matched = n_matched + 1;
    end
end

% The output is width_csv rows with staging columns appended
out_table = [width, staging_rows];
fprintf('Matched %d/%d rows.\n', n_matched, n_width);

writetable(out_table, out_csv);
fprintf('Saved to %s\n', out_csv);

% Also save as excel
[~, name, ~] = fileparts(out_csv);
out_excel = fullfile('./files', [name, '.xlsx']);
writetable(out_table, out_excel);
fprintf('Saved to %s\n', out_excel);

% -------------------------------------------------------------------------
% Comparison plot: sna width (staging) vs Width_mean (width_csv)
% -------------------------------------------------------------------------

x = out_table.("sna width");   % staging measurement
y = out_table.Width_mean;      % new width measurement
x = x(:);  y = y(:);   % corr() requires column vectors

% Only use rows where both values are valid (matched rows)
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);
n = sum(valid);

% --- Stats ---
r       = corr(x, y);
rmse    = sqrt(mean((x - y).^2));
bias    = mean(y - x);          % Bland-Altman bias
sd_diff = std(y - x);
loa_lo  = bias - 1.96 * sd_diff;
loa_hi  = bias + 1.96 * sd_diff;

fprintf('\n--- Comparison stats (n=%d matched rows) ---\n', n);
fprintf('  Pearson r          : %.4f\n', r);
fprintf('  RMSE               : %.4f\n', rmse);
fprintf('  Bias (mean diff)   : %.4f\n', bias);
fprintf('  95%% LoA            : [%.4f, %.4f]\n', loa_lo, loa_hi);

% --- Figure ---
fig = figure('Name', 'sna width vs Width_mean', 'Position', [100 100 900 420]);

% -- Subplot 1: scatter with identity line --
scatter(x, y, 40, 'filled', 'MarkerFaceAlpha', 0.6); hold on;
lims = [min([x; y]), max([x; y])];
plot(lims, lims, 'k--', 'LineWidth', 1.2);
xlabel('sna width (staging)');
ylabel('Width\_mean (width csv)');
title(sprintf('Scatter  |  r = %.3f,  RMSE = %.3f', r, rmse));
axis('equal'); grid('on'); box('on');

sgtitle(fig, 'sna width vs Width\_mean comparison');

% Save figure
fig_path = './files/width_comparison.png';
exportgraphics(fig, fig_path, 'Resolution', 150);
fprintf('Figure saved to %s\n', fig_path);

