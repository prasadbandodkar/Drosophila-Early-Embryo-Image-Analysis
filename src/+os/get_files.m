function Files = get_files(folderpath, filetype, subfolders, substring)
% GET_FILES Extracts files with specific extensions from a directory
%
% Input:
%   folderpath:  Directory to search in
%   filetype:    String or string array of extensions (e.g., "lsm" or ".lsm")
%   subfolders:  (Optional) Search subfolders recursively. Default: true.
%
% Output:
%   Files: String column vector of full file paths

arguments
    folderpath  (1,1) string
    filetype    (1,:) string = "lsm"
    subfolders  (1,1) logical = true
    substring   (1,1) string = ""
end

midstr = "**";
if ~subfolders, midstr = ""; end

Files = [];
for ext = strip(filetype, 'left', '.')
    info = dir(fullfile(folderpath, midstr, "*." + ext));
    info = info(~[info.isdir]);
    if ~isempty(substring) && substring ~= ""
        info = info(arrayfun(@(f) contains(f.name, substring), info));
    end
    if ~isempty(info)
        Files = [Files; fullfile(string({info.folder}'), string({info.name}'))];
    end
end

end