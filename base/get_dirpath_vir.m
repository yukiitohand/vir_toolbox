function [dirpath] = get_dirpath_vir(volumeid,subdir)
    global vir_env_vars
    if ~isempty(subdir)
        dirpath = joinPath(vir_env_vars.local_VIRrootDir,vir_env_vars.url_local_root,volumeid,subdir);
    else
        dirpath = '';
    end
end