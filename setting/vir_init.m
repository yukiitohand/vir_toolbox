function vir_init()
    global vir_env_vars
    str = fileread('vir_toolbox.json');
    vir_env_vars = jsondecode(str);
    
    check_mkdir(vir_env_vars.local_VIRrootDir);
    
    vir_env_vars.url_local_root = vir_env_vars.remote_VIR_URL;
    vir_env_vars.url_remote_root = vir_env_vars.remote_VIR_URL;
    
    global LUT_SCTIME_LIST
    if exist('LUT_SCTIME_LIST.mat','file')
        tmp = load('LUT_SCTIME_LIST.mat');
        LUT_SCTIME_LIST = tmp.lut_sctime_list;
        clear tmp;
    else
        error('LUT_SCTIME_LIST.mat is missing');
    end

end

function [] = check_mkdir(dirpath)
    exist_flg = exist(dirpath,'dir');
    if exist_flg
        %yesno = 1;
    else
        flg = 1;
        while flg
            prompt = sprintf('%s does not exist. Do you want to create?(y/n)',dirpath);
            ow = input(prompt,'s');
            if any(strcmpi(ow,{'y','n'}))
                flg=0;
            else
                fprintf('Input %s is not valid.\n',ow);
            end
        end
        if strcmpi(ow,'n')
            fprintf('No local database will be created.\n');
        elseif strcmpi(ow,'y')
            fprintf('%s is created...\n',dirpath);
        end
    end
end