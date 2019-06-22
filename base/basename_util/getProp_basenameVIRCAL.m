function [prop] = getProp_basenameVIRCAL(basenameVIRCAL,varargin)
% [prop] = getProp_basenameVIRCAL(basenameVIRCAL,varargin)
%   Get properties from the basename of VIR calibration data
%  Input Parameters
%    basenameVIR: 
%    DAWN_sss_N-------N_Vz
%    VIR    : indicates the instrument, fixed
%    sss    : indicates the sensor, IR or VIS for the infrared or 
%             visible spectrometers respectively {IR,VIS}
%    N---N  : indicates name of the calibration file
%    z      : indicates the data file version (1-9)
%  Output Parameters
%   prop: struct storing properties
%    'sensor'
%    'name'   
%    'version'            

[ prop_ori ] = create_propVIRCALbasename();
[basenameptrn] = get_basenameVIRCAL_fromProp(prop_ori);

prop = regexpi(basenameVIRCAL,basenameptrn,'names');

end