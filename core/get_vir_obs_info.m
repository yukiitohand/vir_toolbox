function [obs_info] = get_vir_obs_info(sctime,varargin)

dwld = 0;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case {'DWLD','DOWNLOAD'}
                dwld = varargin{i+1};
            otherwise
                % Something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

global LUT_SCTIME_LIST

sctimes = [LUT_SCTIME_LIST.sctime];
ii = find(sctime==sctimes);

%
volumeid_VIS_1A = load_field_val(LUT_SCTIME_LIST(ii),'volumeid_VIS_1A');
volumeid_VIS_1B = load_field_val(LUT_SCTIME_LIST(ii),'volumeid_VIS_1B');
volumeid_IR_1A  = load_field_val(LUT_SCTIME_LIST(ii),'volumeid_IR_1A');
volumeid_IR_1B  = load_field_val(LUT_SCTIME_LIST(ii),'volumeid_IR_1B');
subdir_VIS_1A   = load_field_val(LUT_SCTIME_LIST(ii),'subdir_VIS_1A');
subdir_VIS_1B   = load_field_val(LUT_SCTIME_LIST(ii),'subdir_VIS_1B');
subdir_IR_1A    = load_field_val(LUT_SCTIME_LIST(ii),'subdir_IR_1A');
subdir_IR_1B    = load_field_val(LUT_SCTIME_LIST(ii),'subdir_IR_1B');
basename_VIS_1A = load_field_val(LUT_SCTIME_LIST(ii),'basename_VIS_1A');
basename_VIS_1B = load_field_val(LUT_SCTIME_LIST(ii),'basename_VIS_1B');
basename_IR_1A  = load_field_val(LUT_SCTIME_LIST(ii),'basename_IR_1A');
basename_IR_1B  = load_field_val(LUT_SCTIME_LIST(ii),'basename_IR_1B');

if ~isempty(volumeid_IR_1A)
    [volumeid_info] = decompose_volumeid(volumeid_IR_1A);
elseif ~isempty(volumeid_VIS_1A)
    [volumeid_info] = decompose_volumeid(volumeid_VIS_1A);
end

mission = get_mission(volumeid_info.mission_id);

if ~isempty(subdir_IR_1A)
    [mission_phase_period] = get_mission_phases(subdir_IR_1A);
elseif ~isempty(subdir_VIS_1A)
    [mission_phase_period] = get_mission_phases(subdir_VIS_1A);
end


if isempty(volumeid_VIS_1B) && ~isempty(volumeid_VIS_1A)
    vol_info1B = decompose_volumeid(volumeid_VIS_1A);
    vol_info1B.level = '1B';
    volumeid_VIS_1B = get_volumeid(vol_info1B);
end

if isempty(volumeid_IR_1B) && ~isempty(volumeid_IR_1A)
    vol_info1B = decompose_volumeid(volumeid_IR_1A);
    vol_info1B.level = '1B';
    volumeid_IR_1B = get_volumeid(vol_info1B);
end


[dirpath_vis1A] = get_dirpath_vir(volumeid_VIS_1A,subdir_VIS_1A);
[dirpath_vis1B] = get_dirpath_vir(volumeid_VIS_1B,subdir_VIS_1B);
[dirpath_ir1A]  = get_dirpath_vir(volumeid_IR_1A, subdir_IR_1A);
[dirpath_ir1B]  = get_dirpath_vir(volumeid_IR_1B, subdir_IR_1B);

% get HK basenames
[basename_VIS_1A_HK] = get_basenameHK(basename_VIS_1A);
[basename_VIS_1B_HK] = get_basenameHK(basename_VIS_1B);
[basename_IR_1A_HK]  = get_basenameHK(basename_IR_1A);
[basename_IR_1B_HK]  = get_basenameHK(basename_IR_1B);

% get QQ basenames
[basename_VIS_1B_QQ] = get_basenameQQ(basename_VIS_1B);
[basename_IR_1B_QQ]  = get_basenameQQ(basename_IR_1B);

%download
if dwld>0
    vsubdir_VIS_1A = joinPath(volumeid_VIS_1A,subdir_VIS_1A);
    vsubdir_VIS_1B = joinPath(volumeid_VIS_1B,subdir_VIS_1B);
    vir_downloader_wrapper(vsubdir_VIS_1A,basename_VIS_1A,   dwld);
    vir_downloader_wrapper(vsubdir_VIS_1A,basename_VIS_1A_HK,dwld);
    vir_downloader_wrapper(vsubdir_VIS_1B,basename_VIS_1B,   dwld);
    vir_downloader_wrapper(vsubdir_VIS_1B,basename_VIS_1B_HK,dwld);
    vir_downloader_wrapper(vsubdir_VIS_1B,basename_VIS_1B_QQ,dwld);
    %
    vsubdir_IR_1A = joinPath(volumeid_IR_1A,subdir_IR_1A);
    vsubdir_IR_1B = joinPath(volumeid_IR_1B,subdir_IR_1B);
    vir_downloader_wrapper(vsubdir_IR_1A,basename_IR_1A,    dwld);
    vir_downloader_wrapper(vsubdir_IR_1A,basename_IR_1A_HK, dwld);
    vir_downloader_wrapper(vsubdir_IR_1B,basename_IR_1B,    dwld);
    vir_downloader_wrapper(vsubdir_IR_1B,basename_IR_1B_HK, dwld);
    vir_downloader_wrapper(vsubdir_IR_1B,basename_IR_1B_QQ, dwld);
    if ~isempty(volumeid_VIS_1B)
        vir_downloader(joinPath(volumeid_VIS_1B,'CALIB'),'dwld', dwld);
    end
    if ~isempty(volumeid_IR_1B)
        vir_downloader(joinPath(volumeid_IR_1B, 'CALIB'),'dwld', dwld);
    end
end

% get calibration
if ~isempty(volumeid_VIS_1B)
    dirpath_cal_VIS_1B  = get_dirpath_vir(volumeid_VIS_1B,'CALIB');
else
    dirpath_cal_VIS_1B = '';
end
if ~isempty(volumeid_IR_1B)
    dirpath_cal_IR_1B   = get_dirpath_vir(volumeid_IR_1B, 'CALIB');
else
    dirpath_cal_IR_1B = '';
end
basename_cal_VIS_RESP    = 'DAWN_VIR_VIS_RESP_V2';
basename_cal_IR_RESP     = 'DAWN_VIR_IR_RESP_V2';
basename_cal_VIS_SOLSPEC = 'DAWN_VIR_VIS_SOLAR_SPECTRUM_V2';
basename_cal_IR_SOLSPEC  = 'DAWN_VIR_IR_SOLAR_SPECTRUM_V2';



obs_info = [];
obs_info.sctime              = sctime;
obs_info.mission             = mission;
obs_info.mission_id          = volumeid_info.mission_id;
obs_info.mission_phase       = mission_phase_period.phase_name;
obs_info.mission_phase_date  = mission_phase_period.phase_date;
obs_info.mission_period      = mission_phase_period.period_name;
obs_info.mission_period_date = mission_phase_period.period_date;
obs_info.volumeid_vis_1A     = volumeid_VIS_1A;
obs_info.volumeid_vis_1B     = volumeid_VIS_1B;
obs_info.volumeid_ir_1A      = volumeid_IR_1A;
obs_info.volumeid_ir_1B      = volumeid_IR_1B;
obs_info.subdir_vis_1A       = subdir_VIS_1A;
obs_info.subdir_vis_1B       = subdir_VIS_1B;
obs_info.subdir_ir_1A        = subdir_IR_1A;
obs_info.subdir_ir_1B        = subdir_IR_1B;
obs_info.dirpath_vis_1A      = dirpath_vis1A;
obs_info.dirpath_vis_1B      = dirpath_vis1B;
obs_info.dirpath_ir_1A       = dirpath_ir1A;
obs_info.dirpath_ir_1B       = dirpath_ir1B;
obs_info.basename_vis_1A     = basename_VIS_1A;
obs_info.basename_vis_1B     = basename_VIS_1B;
obs_info.basename_ir_1A      = basename_IR_1A;
obs_info.basename_ir_1B      = basename_IR_1B;
obs_info.basename_vis_1A_hk  = basename_VIS_1A_HK;
obs_info.basename_vis_1B_hk  = basename_VIS_1B_HK;
obs_info.basename_ir_1A_hk   = basename_IR_1A_HK;
obs_info.basename_ir_1B_hk   = basename_IR_1B_HK;
obs_info.basename_vis_1B_qq  = basename_VIS_1B_QQ;
obs_info.basename_ir_1B_qq   = basename_IR_1B_QQ;

obs_info.dirpath_cal_vis_1B       = dirpath_cal_VIS_1B;
obs_info.dirpath_cal_ir_1B        = dirpath_cal_IR_1B;
obs_info.basename_cal_vis_resp    = basename_cal_VIS_RESP;
obs_info.basename_cal_ir_resp     = basename_cal_IR_RESP;
obs_info.basename_cal_vis_solspec = basename_cal_VIS_SOLSPEC;
obs_info.basename_cal_ir_solspec  = basename_cal_IR_SOLSPEC;



end

function [fldval] = load_field_val(struct_obj,fldnm)
    if ~isempty(struct_obj.(fldnm)) 
        fldval = struct_obj.(fldnm);
    else
        fldval = '';
    end
end

function [] = vir_downloader_wrapper(subdir,basename,dwld)
    if ~isempty(basename)
        vir_downloader(subdir,'BASENAMEPTRN',basename,'dwld',dwld);
    end
end

