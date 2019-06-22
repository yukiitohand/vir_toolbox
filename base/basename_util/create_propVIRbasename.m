function [ prop ] = create_propVIRbasename( varargin )
% [ prop ] = create_propVIRbasename( varargin )
%   return struct of VIR data property
% 
%   Output
%   prop: struct storing properties
%    'sensor'                : (default) '(?<sensor>VIS|IR)'
%    'level'                 : (default) '(?<level>[\d]{1}[a-zA-Z]{1})'
%    'clock_reset_number'    : (default) '(?<clock_reset_number>[\d]{1})'
%    'sctime'                : (default) '(?<sctime>[\d]{9})'
%    'type'                  : (default) '_*(?<type>HK|QQ){0,1}_*'
%    'version'               : (default) '(?<version>[0-9]{1}){0,1}'
%   Optional Parameters
%    'SENSOR','LEVEL','CLOCK_RESET_NUMBER','SCTIME','VERSION'
% 
%  * Reference *
%    VIR_sss_ll_v_sctime[_TT_]_z
%    VIR    : indicates the instrument, fixed
%    sss    : indicates the sensor, IR or VIS for the infrared or 
%             visible spectrometers respectively {IR,VIS}
%    ll     : indicates the processing level, either 1A or 1B
%    v      : indicates the clock reset number
%    sctime : is the acquisition SC_CLOCK_START_COUNT (integer part)
%    type   : additional information for data (empty for QUB data, HK for
%             house keeping table, QQ for wavelength related information)
%    z      : indicates the data file version (1-9)

sensor             = '(?<sensor>VIS|IR)';
level              = '(?<level>[\d]{1}[a-zA-Z]{1})';
clock_reset_number = '(?<clock_reset_number>[\d]{1})';
sctime             = '(?<sctime>[\d]{9})';
tp                 = '_*(?<type>HK|QQ){0,1}_*';
vr                 ='(?<version>[0-9]{1}){0,1}';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'SENSOR'
                sensor = varargin{i+1};
            case 'LEVEL'
                level = varargin{i+1};
            case 'CLOCK_RESET_NUMBER'
                clock_reset_number = varargin{i+1};
            case 'SCTIME'
                sctime = varargin{i+1};
            case 'TYPE'
                tp = varargin{i+1};
            case 'VERSION'
                vr = varargin{i+1};               
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);   
        end
    end
end

prop = [];
prop.sensor             = sensor;
prop.level              = level;
prop.clock_reset_number = clock_reset_number;
prop.sctime             = sctime;
prop.type               = tp; 
prop.version            = vr;

end

