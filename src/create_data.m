function data = create_data(file, cfg)

    arguments
        file (1,1) string
        cfg  struct
    end

    % Use cfg.channel_order to create the data struct
    switch cfg.data_type
        case 'dgrp'
            g = struct();
            for i = 1:numel(cfg.channel_names)
                if cfg.skip_channels(i)
                    continue
                end
                gname = char(cfg.channel_names(i));
                g.(gname) = structs.gene( ...
                    name=cfg.channel_names(i),         ...
                    channel_type=cfg.channel_types(i), ...
                    channel_id=int32(i),               ...
                    t=[], ...
                    raw=[], ...
                    s=[] ...
                );
            end

            data = structs.data(filename=file, genes=g);
            if cfg.info
                fprintf('Data struct created for %s\n', file)
            end
        
        otherwise
            data = structs.data(filename=file);
            if cfg.info
                fprintf('Empty data struct created for %s\n', file)
            end
    end

end

