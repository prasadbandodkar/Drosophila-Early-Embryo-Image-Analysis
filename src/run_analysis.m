function data = run_analysis(file, cfg)

    arguments
        file (1,1) string
        cfg  struct
    end
    % filename (for .mat naming) from fileparts
    [~, filename, ~] = fileparts(file);
    % path_data folder built from last cfg.num_file_parts path segments
    parts = strsplit(file, filesep);
    parts = parts(max(1, end - cfg.num_file_parts + 1) : end);
    [~, stem, ~] = fileparts(parts{end});
    parts{end}    = stem;
    cfg.path_data = fullfile(cfg.path_data, strjoin(parts, '_'));
    if ~isfolder(cfg.path_data)
        mkdir(cfg.path_data);
    elseif cfg.overwrite_data
        rmdir(cfg.path_data, 's');
        mkdir(cfg.path_data);
    end

    % custom data based on cfg
    data = create_data(file, cfg);

    % mat file path (written by analyze, read back if available)
    mat_path = fullfile(cfg.path_data, filename + ".mat");

    % run image analysis – skip if a cached .mat already exists
    if isfile(mat_path)
        if cfg.info
            fprintf('\nLoading cached analysis from %s\n', mat_path);
        end
        tmp  = load(mat_path, 'data');
        data = tmp.data;
    else
        if cfg.info
            fprintf('\nRunning image analysis for %s\n', file);
        end
        data = analyze(data, cfg);
        if cfg.save_analysis
            if cfg.info
                fprintf('Saving analysis data to: %s\n', mat_path);
            end
            save(mat_path, 'data');
        end
    end

    % run fitting
    if cfg.info
        fprintf('\nRunning fitting for %s\n', file);
    end
    data = fitting(data, cfg);

    % save results
    if cfg.info
        fprintf('Saving fitting data to: %s\n', mat_path);
    end
    save(mat_path, 'data');
end
