function [ prop ] = create_propVIRCALbasename( varargin )
% [ prop ] = create_propVIRCALbasename( varargin )
%   return struct of VIR CALIBRATION data property
% 
%   Output
%   prop: struct storing properties
%    'sensor'         : (default) '(?<sensor>VIS|IR)'
%    'name'           : (default) '(?<name>[0-9a-zA-Z_]*[A-Z0-9]{1})'
%    'version'        : (default) '(?<version>V[0-9]{1}){1}'
%   Optional Parameters
%    'SENSOR','NAME','VERSION'
% 
%  * Reference *
%    DAWN_sss_N-------N_Vz
%    VIR    : indicates the instrument, fixed
%    sss    : indicates the sensor, IR or VIS for the infrared or 
%             visible spectrometers respectively {IR,VIS}
%    N---N  : indicates name of the calibration file
%    z      : indicates the data file version (1-9)

sensor  = '(?<sensor>VIS|IR)';
name    = '(?<name>[0-9a-zA-Z_]*[A-Z0-9]{1})';
vr      = '(?<version>V[0-9]{1}){1}';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'SENSOR'
                sensor = varargin{i+1};
            case 'NAME'
                name = varargin{i+1};
            case 'VERSION'
                vr = varargin{i+1};               
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);   
        end
    end
end

prop = [];
prop.sensor   = sensor;
prop.name     = name;
prop.version  = vr;

end