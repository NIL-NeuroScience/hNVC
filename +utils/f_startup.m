
path = fileparts(fileparts(mfilename('fullpath')));
addPaths(path);

function addPaths(mainDir)
    addDirs = genpath(mainDir);
    addDirs = strsplit(addDirs, pathsep);
    addDirs = addDirs(~cellfun(@(x) any(startsWith(strsplit(x, ...
        filesep), '.')), addDirs));
    addDirs = strjoin(addDirs, pathsep);
    addpath(addDirs);
end