function path = f_getPackagePath()

    path = fileparts(fileparts(mfilename('fullpath')));

end