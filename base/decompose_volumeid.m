function [volumeid_info] = decompose_volumeid(volumeid)
% Decompose volume id of VIR data archive.
% INPUTS:
%   volumeid: something like 'DWNVVIR_I1B'
% OUTPUTS
%   volumeid_info: struct having fields
%     mission_id  : {'C','V','X'} 'C' for Ceres, 'V' for Vesta, and 'X' for
%                   cruise
%     mission_phase_id : string, can be empty.
%     sensor_id   : {'I','V'} 'I' for IR and 'V' for VIS
%     level       : {'1A','1B'} 
%     version     : 'v2' or '' (for version 1)

volume_id_ptrn = 'DWN(?<mission_id>[A-Z]{1})(?<mission_phase_id>[0-9A-Z]*)VIR_(?<sensor_id>I|V*)(?<level>[0-9]{1}[A-Z]{1})_*(?<version>v[0-9]{1})*';
volumeid_info = regexpi(volumeid,volume_id_ptrn,'names');

end