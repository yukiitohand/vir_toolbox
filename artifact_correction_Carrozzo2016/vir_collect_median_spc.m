function [img_spcs_meds,idxes_valid,isinc_lt_theta_list] = vir_collect_median_spc(sctime_info,varargin)
% [img_spcs_meds,idxes_valid,isinc_lt_50_list] = vir_collect_median_spc(sctime_info,varargin)
% INPUTS
%  sctime_info: subset of LUT_SCTIME_LIST
% OUTPUTS
%  img_spc_meds: median spectra for each column (1 x S x B)
%  idxes_valid : valid indices for sctime_info the images associated with
%                which are used for computing median spectra. I noticed that
%                the number of lines of some images is inconsistent with the
%                number of lines in their house keeping table, 
%                not straightforward to get sctime for each line. 
%  ininc_lt_theta_list
%              : [ L' x S ], boolean matrix, storing whether or not the
%                incident angle is lower than theta deg or not. L' refers to
%                the number of all the lines stucked for computing median
%                spectra.
%  
%  OPTIONAL PARAMETERS
%   'Theta'   : threhsold value [degree] for incident angles. 
%               The pixels greater than this value will be marked with false,
%               otherwise true.
%               (default) : 50
%   'SENSOR_ID'  : {'V','I'} (default) : 'I' 

%%
theta = 50;
sensor_id = 'I';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'THETA'
                theta = varargin{i+1};
            case 'SENSOR_ID'
                sensor_id = varargin{i+1};
            otherwise
                error('Undefined option:%s.',upper(varargin{i}));
        end
    end
end

switch upper(sensor_id)
    case 'I'
        bands = 432;
        samples = 256;
    case 'V'
        error('not implemented yet');
    otherwise
        error('Invalid sensor_id %s',sensor_id);
end

%% evaluate 
[idxes_valid,isinc_lt_theta_list] = vir_isinc_lt_theta(sctime_info,'Theta',theta,'SENSOR_ID',sensor_id);

%% compute median spectra for each column
img_spcs_meds = nan(1,samples,bands);
for c = 1:samples
    tic;
    [I1Bimc_list_med,~] = vir_get_median_spc_column(sctime_info,...
                           idxes_valid,isinc_lt_theta_list,c);
    img_spcs_meds(1,c,:) = I1Bimc_list_med;
    toc;
end

end