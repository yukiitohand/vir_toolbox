
[dirs,files] = vir_downloader('','dwld',1);

for i=1:length(dirs)
    dname = dirs{i};
    if ~isempty(regexpi(dname,'^DWN','ONCE'))
        subdir = joinPath(dname,'INDEX');
        vir_downloader(subdir,'dwld',2);
    end
end

dirs1 = [];
for i=1:length(dirs)
    dname = dirs{i};
    if ~isempty(regexpi(dname,'^DWN.*_(I|V)1(A|B).*','ONCE'))
        dirs1 = [dirs1,{dname}];
    end
end

dirs1_a = [];
for i=1:length(dirs1)
    dname = dirs1{i};
    t = ~isempties(regexpi(dirs1,dname,'ONCE'));
    t = find(t);
    if length(t) == 1 && t==i
        dirs1_a = [dirs1_a,{dname}];
    else
    end
end

%%
baseIdx = 'INDEX';
for i=1:length(dirs1_a)
    dname = dirs1_a{i};
    subdir = joinPath(dname,'INDEX');
    dpath = joinPath(vir_env_vars.local_VIRrootDir,vir_env_vars.url_local_root,subdir);
    idxlbl = crismlblread_v2(joinPath(dpath,[baseIdx '.LBL']));
    idxtab = crismTABread(joinPath(dpath,[baseIdx '.TAB']),idxlbl,'skip_line',1);
    if isfield(idxtab.data,'IMAGE_MID_TIME') && strcmpi(idxtab.data(1).IMAGE_MID_TIME(1),',')
    % if strcmpi(dname,'DWNXVIR_I1B')
        fprintf('%s\n',dname);
        idxtab.data(10)
    end
end

%%
baseIdx = 'INDEX';
lut_sctime = [];
for i=1:length(dirs1_a)
    dname = dirs1_a{i};
    subdir = joinPath(dname,'INDEX');
    dpath = joinPath(vir_env_vars.local_VIRrootDir,vir_env_vars.url_local_root,subdir);
    idxlbl = crismlblread_v2(joinPath(dpath,[baseIdx '.LBL']));
    idxtab = crismTABread(joinPath(dpath,[baseIdx '.TAB']),idxlbl,'skip_line',1);
    volumeid = dname;
    for k=1:length(idxtab.data)
        if ~isempty(idxtab.data(k).PRODUCT_ID)
            [dsubpath,basename,ext] = fileparts(idxtab.data(k).FILE_SPECIFICATION_NAME);
            prop = getProp_basenameVIR(basename);
            if ~isempty(prop)
                sctime = prop.sctime;
                if strcmpi(dsubpath(1),'/')
                    dsubpath = dsubpath(2:end);
                end
                % volumeid = idxtab.data(k).VOLUME_ID;
                % I found the entry VOLUME_ID may be wrong for _v2 volumes.
                %volumeid = dname;
%                 mission_id_struct = regexpi(volumeid,'DWN(?<mission_id>[A-Z]{1})(?<mission_phase_id>[0-9A-Z]*)VIR_(?<sensor>I|V*)(?<level>[0-9]{1}[A-Z]{1})_*(?<version>v[0-9]{1})*','names');
%                 mission_id = mission_id_struct.mission_id;
%                 mission_phase_id = mission_id_struct.mission_phase_id;
%                 sensor = mission_id_struct.sensor;
%                 level = mission_id_struct.level;
%                 volume_ver = mission_id_struct.version;
%                 switch upper(mission_id)
%                     case 'C'
%                         mission = 'Ceres';
%                     case 'V'
%                         mission = 'Vesta';
%                     case 'X'
%                         mission = 'Cruies';
%                     otherwise
%                         error('Mission is right??');
%                 end
                if isempty(prop.type)
%                     dsubpath_sep = split(dsubpath,'/');
%                     if strcmpi(dsubpath_sep{1},'DATA')
%                         mission_phase_str = dsubpath_sep{2};
%                         mission_phase     = regexpi(mission_phase_str,...
%                             '(?<date>[\d]{8})_(?<phase_name>.*)','names');
%                         if length(dsubpath_sep)>2
%                             mission_period_str = dsubpath_sep{3};
%                             mission_period = regexpi(mission_period_str,...
%                                 '(?<date>[\d]{8})_(?<period_name>.*)','names');
%                         else
%                             mission_period = [];
%                         end
%                     else
%                         error('Probably something wrong with dsubpath');
%                     end
                    fldnm = sprintf('SCTIME%09d',sctime);
                    if ~isfield(lut_sctime,fldnm)
                        lut_sctime.(fldnm) = [];
%                         lut_sctime.(fldnm).subdir_VIS_1A   = '';
%                         lut_sctime.(fldnm).volumeid_VIS_1A = '';
%                         lut_sctime.(fldnm).basename_VIS_1A = '';
%                         lut_sctime.(fldnm).subdir_VIS_1B   = '';
%                         lut_sctime.(fldnm).volumeid_VIS_1B = '';
%                         lut_sctime.(fldnm).basename_VIS_1B = '';
%                         lut_sctime.(fldnm).subdir_IR_1A    = '';
%                         lut_sctime.(fldnm).volumeid_IR_1A  = '';
%                         lut_sctime.(fldnm).basename_IR_1A  = '';
%                         lut_sctime.(fldnm).subdir_IR_1B    = '';
%                         lut_sctime.(fldnm).volumeid_IR_1B  = '';
%                         lut_sctime.(fldnm).basename_IR_1B  = '';
%                         lut_sctime.(fldnm).mission         = mission;
%                         lut_sctime.(fldnm).mission_id      = mission_id;
%                         if isempty(mission_phase)
%                             lut_sctime.(fldnm).mission_phase   = '';
%                             lut_sctime.(fldnm).mission_phase_date   = '';
%                         else
%                             lut_sctime.(fldnm).mission_phase   = mission_phase.phase_name;
%                             lut_sctime.(fldnm).mission_phase_date   = mission_phase.date;
%                         end
%                         if isempty(mission_period)
%                             lut_sctime.(fldnm).mission_period       = '';
%                             lut_sctime.(fldnm).mission_period_date  = '';
%                         else
%                             lut_sctime.(fldnm).mission_period       = mission_period.period_name;
%                             lut_sctime.(fldnm).mission_period_date  = mission_period.date;
%                         end
                    end
                    if strcmpi(prop.level,'1A') && strcmpi(prop.sensor,'VIS')
                        lut_sctime.(fldnm).subdir_VIS_1A = dsubpath;
                        lut_sctime.(fldnm).volumeid_VIS_1A = volumeid;
                        lut_sctime.(fldnm).basename_VIS_1A = basename;
                    elseif strcmpi(prop.level,'1B') && strcmpi(prop.sensor,'VIS')
                        lut_sctime.(fldnm).subdir_VIS_1B = dsubpath;
                        lut_sctime.(fldnm).volumeid_VIS_1B = volumeid;
                        lut_sctime.(fldnm).basename_VIS_1B = basename;
                    elseif strcmpi(prop.level,'1A') && strcmpi(prop.sensor,'IR')
                        lut_sctime.(fldnm).subdir_IR_1A = dsubpath;
                        lut_sctime.(fldnm).volumeid_IR_1A = volumeid;
                        lut_sctime.(fldnm).basename_IR_1A = basename;
                    elseif strcmpi(prop.level,'1B') && strcmpi(prop.sensor,'IR')
                        lut_sctime.(fldnm).subdir_IR_1B = dsubpath;
                        lut_sctime.(fldnm).volumeid_IR_1B = volumeid;
                        lut_sctime.(fldnm).basename_IR_1B = basename;
                    end
                end
            end
        end
    end
end
%save LUT_SCTIME.mat lut_sctime
%%
%load LUT_SCTIME.mat lut_sctime
flds = fieldnames(lut_sctime);
L = length(flds);
sctimes = cellfun(@(x) str2num(x(7:end)),flds);
[vsort,isort] = sort(sctimes);

lut_sctime_list = struct('sctime',cell(L,1),...
    'volumeid_VIS_1A',cell(L,1),'volumeid_VIS_1B',cell(L,1),...
    'volumeid_IR_1A',cell(L,1),'volumeid_IR_1B',cell(L,1),...
    'subdir_VIS_1A',cell(L,1),'subdir_VIS_1B',cell(L,1),...
    'subdir_IR_1A',cell(L,1),'subdir_IR_1B',cell(L,1),...
    'basename_VIS_1A',cell(L,1),'basename_VIS_1B',cell(L,1),...
    'basename_IR_1A',cell(L,1),'basename_IR_1B',cell(L,1));

for ii=1:length(isort)
    lut_sctime_list(ii).sctime = vsort(ii);
    scfield = flds{isort(ii)};
    fldstmp = fieldnames(lut_sctime.(scfield));
    for jj=1:length(fldstmp)
        fldtmp = fldstmp{jj};
        if ~isempty(lut_sctime.(scfield).(fldtmp))
            lut_sctime_list(ii).(fldtmp) = lut_sctime.(scfield).(fldtmp);
        end
    end
end

save LUT_SCTIME_LIST.mat lut_sctime_list
