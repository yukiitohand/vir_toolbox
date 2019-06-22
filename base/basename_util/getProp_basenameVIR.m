function [prop] = getProp_basenameVIR(basenameVIR,varargin)
% [prop] = getProp_basenameVIR(basenameVIR,varargin)
%   Get properties from the basename of VIR
%  Input Parameters
%    basenameVIR: 
%     string like  "VIR_sss_ll_v_sctime[_TT_]_z"
%     where,
%      VIR    : indicates the instrument, fixed
%      sss    : indicates the sensor, IR or VIS for the infrared or 
%               visible spectrometers respectively {IR,VIS}
%      ll     : indicates the processing level, either 1A or 1B
%      v      : indicates the clock reset number
%      sctime : is the acquisition SC_CLOCK_START_COUNT (integer part)
%      type   : additional information for data (empty for QUB data, HK for
%               house keeping table, QQ for wavelength related information)
%      z      : indicates the data file version (1-9)
%  Output Parameters
%    prop: struct storing properties
%      'sensor'               
%      'level'                
%      'clock_reset_number'   
%      'sctime'               
%      'TYPE'                 
%      'version'              


[ prop_ori ] = create_propVIRbasename();
[basenameptrn] = get_basenameVIR_fromProp(prop_ori);

prop = regexpi(basenameVIR,basenameptrn,'names');
if ~isempty(prop)
    prop.clock_reset_number = str2num(prop.clock_reset_number);
    prop.sctime = str2num(prop.sctime);
    prop.version = str2num(prop.version);
end

end