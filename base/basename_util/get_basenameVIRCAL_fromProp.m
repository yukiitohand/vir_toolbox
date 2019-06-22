function [basenameVIRCAL] = get_basenameVIRCAL_fromProp(prop)
% [basenameVIRCAL] = get_basenameVIRCAL_fromProp(prop)
%   return basename for given VIR Calibration data property
%   INPUTS
%   prop: struct storing properties
%    'sensor'
%    'name'   
%    'version'
%
%   OUTPUTS: 
%    basenameVIRCAL: 
%    DAWN_sss_N-------N_Vz
%    VIR    : indicates the instrument, fixed
%    sss    : indicates the sensor, IR or VIS for the infrared or 
%             visible spectrometers respectively {IR,VIS}
%    N---N  : indicates name of the calibration file
%    z      : indicates the data file version (1-9)

sensor   = prop.sensor;
name     = prop.name;
vr       = prop.version;


if any(strcmpi(sensor,{'VIS','IR'}))
    sensor = upper(sensor);
else
    %error('property has incorrect sensor information');
end

basenameVIRCAL = sprintf('DAWN_VIR_%s_%s_%s',sensor,name,vr);

end