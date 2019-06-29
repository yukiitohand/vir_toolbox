function [sctime_info] = get_sctime_info(volumeid)

global LUT_SCTIME_LIST



vol_info = decompose_volumeid(volumeid);

switch upper(vol_info.level)
    case '1A'
        switch upper(vol_info.sensor_id)
            case 'I'
                volumeid_IR_1A = {LUT_SCTIME_LIST.volumeid_IR_1A};
                is_volumeid = strcmpi(volumeid_IR_1A,volumeid);
            case 'V'
                volumeid_VIS_1A = {LUT_SCTIME_LIST.volumeid_VIS_1A};
                is_volumeid = strcmpi(volumeid_VIS_1A,volumeid);
        end
    case '1B'
        switch upper(vol_info.sensor_id)
            case 'I'
                volumeid_IR_1B = {LUT_SCTIME_LIST.volumeid_IR_1B};
                is_volumeid = strcmpi(volumeid_IR_1B,volumeid);
            case 'V'
                volumeid_VIS_1B = {LUT_SCTIME_LIST.volumeid_VIS_1B};
                is_volumeid = strcmpi(volumeid_VIS_1B,volumeid);
        end
end

sctime_info = LUT_SCTIME_LIST(is_volumeid);

end