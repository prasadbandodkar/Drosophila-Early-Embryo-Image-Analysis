% script_main
% This is the main script to analyze the fixed image dataset for staging.
% Cmd+Shift+B — runs the task directly (it's set as the default build task)
% Cmd+Shift+P → "Tasks: Run Task" → select "Run MATLAB (batch)"

clear
clc
close all

addpath(fullfile(fileparts(mfilename('fullpath')), 'src'))
addpath(fullfile(fileparts(mfilename('fullpath')), 'bfmatlab'))

% load configuration
cfg = stagingconfig2();

% Create timestamped output folder and start diary there
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
if ~cfg.subfolders
    cfg.path_data = fullfile(cfg.path_data, timestamp);
end
if ~isfolder(cfg.path_data)
    mkdir(cfg.path_data);
end
log_file = fullfile(cfg.path_data, sprintf('run_%s.log', timestamp));
diary(log_file)
fprintf('Log file: %s\n', log_file)

disp(cfg)

% First get all the files defined by the config file
files = [];
for i = 1:length(cfg.path_raw)
    files = [files; os.get_files(cfg.path_raw(i), cfg.filetype, cfg.subfolders)];
end
fprintf('Total number of files found: %d\n', length(files))
files = files(contains(files, cfg.substring));
fprintf('Total number of files after filtering for %s: %d\n', cfg.substring, length(files))


% Use parallel processing if requested and PCT is available
has_pct = license('test', 'Distrib_Computing_Toolbox');
if cfg.run_parallel && ~has_pct
    fprintf('Parallel processing requested but Parallel Computing Toolbox is not available. Running serially.\n')
end
if cfg.run_parallel && has_pct
    fprintf('Using parallel processing (%d workers)\n', cfg.num_workers)
    pool = gcp('nocreate');
    if isempty(pool)
        parpool(cfg.num_workers);
    end
    parfor i = 1:length(files)
        file = files(i);
        fprintf('\n\n\nProcessing file %d/%d: %s\n', i, length(files), file)
        t0 = tic;
        try
            run_analysis(file, cfg);
        catch ME
            fprintf('Error processing file %s: %s\n  Location: %s (line %d)\n', ...
                file, ME.message, ME.stack(1).file, ME.stack(1).line)
        end
        fprintf('  Time elapsed: %.1f s\n', toc(t0))
    end
else
    fprintf('Running serially\n')
    for i = 1:length(files)
        file = files(i);
        fprintf('\n\n\nProcessing file %d/%d: %s\n', i, length(files), file)
        tic
        try
            run_analysis(file, cfg);
        catch ME
            fprintf('Error processing file %s: %s\n  Location: %s (line %d)\n', ...
                file, ME.message, ME.stack(1).file, ME.stack(1).line)
        end
        fprintf('  Time elapsed: %.1f s\n', toc)
    end
end

diary off