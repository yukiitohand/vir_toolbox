function [basenameVIR] = get_basenameVIR_fromProp(prop)
% [basenameVIR] = get_basenameVIR_fromProp(prop)
%   return basename for given VIR data property
%   INPUTS
%    prop: struct storing properties
%      'sensor'                
%      'level'                 
%      'clock_reset_number'    
%      'sctime'                
%      'TYPE'                  
%      'version'               
%
%   OUTPUTS: 
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

sensor             = prop.sensor;
level              = prop.level;
clock_reset_number = prop.clock_reset_number;
sctime             = prop.sctime;
tp                 = prop.type; 
vr                 = prop.version;

if any(strcmpi(sensor,{'VIS','IR'}))
    sensor = upper(sensor);
else
    %error('property has incorrect sensor information');
end

if isnumeric(clock_reset_number)
    clock_reset_number = sprintf('%1d',clock_reset_number);
end
if isnumeric(sctime)
    sctime = sprintf('%09d',sctime);
end
if isnumeric(vr)
    vr = sprintf('%1d',vr);
end



if isempty(tp)
    basenameVIR = sprintf('VIR_%s_%s_%s_%s_%s',sensor,level,...
        clock_reset_number,sctime,vr);
elseif any(strcmpi(tp,{'QQ','HK'}))
    basenameVIR = sprintf('VIR_%s_%s_%s_%s_%s_%s',sensor,level,...
        clock_reset_number,sctime,tp,vr);
else 
    basenameVIR = sprintf('VIR_%s_%s_%s_%s%s%s',sensor,level,...
        clock_reset_number,sctime,tp,vr);
end

end