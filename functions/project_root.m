function root = project_root()
%PROJECT_ROOT  Absolute path of the project folder (one level above functions/).
root = fileparts(fileparts(mfilename('fullpath')));
end
