function [volumeid] = get_volumeid(volumeid_info)
% get volumeid from volumeid_info struct
% INPUTS:
%   volumeid_info: struct having fields
%     mission_id  : {'C','V','X'} 'C' for Ceres, 'V' for Vesta, and 'X' for
%                   cruise
%     mission_phase_id : string, can be empty.
%     sensor_id   : {'I','V'} 'I' for IR and 'V' for VIS
%     level       : {'1A','1B'} 
%     version     : 'v2' or '' (for version 1)
% OUTPUTS
%   volumeid: something like 'DWNVVIR_I1B'
%   

if isempty(volumeid_info.version)
    volumeid = sprintf('DWN%s%sVIR_%s%s',...
        volumeid_info.mission_id,volumeid_info.mission_phase_id,...
        volumeid_info.sensor_id,volumeid_info.level);
else
    volumeid = sprintf('DWN%s%sVIR_%s%s_%s',...
        volumeid_info.mission_id,volumeid_info.mission_phase_id,...
        volumeid_info.sensor_id,volumeid_info.level,volumeid_info.version);


end