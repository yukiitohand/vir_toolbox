function [idxes_valid,isinc_lt_theta_list,img_sctime_list] = vir_isinc_lt_theta(sctime_info,varargin)
% [idxes_valid,isinc_lt_50_list] = vir_isinc_lt_theta(sctime_info,varargin)
%  Evaluate incident angles at each pixels for given sets of images
%  specified as (sctim_info). The resulting images are row stucked.
%  Some invalid images are not used. The pixels with their incident angles
%  lower than theta will be marked as true, otherwise false.
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
%                incident angle is lower than theta or not. L' refers to
%                the number of all the lines stucked for computing median
%                spectra.
%  Optional Parameters
%  'Theta'  : threhsold value [degree] for incident angles. 
%             The pixels greater than this value will be marked with false,
%             otherwise true.
%             (default) : 50 
%   'SENSOR_ID'  : {'V','I'} (default) : 'I' 

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
        basename_1B_fld = 'basename_ir_1B';
        dirpath_1B_fld  = 'dirpath_ir_1B';
    case 'V'
        error('not implemented yet');
    otherwise
        error('Invalid sensor_id %s',sensor_id);
end

idxBool_valid = true(length(sctime_info),1);
for i=1:length(sctime_info)
    sctime_info_idx = sctime_info(i);
    sctime = sctime_info_idx.sctime;
    vir_obs = get_vir_obs_info(sctime,'dwld',0);
    vir1Bdata = VIRdata(vir_obs.(basename_1B_fld),vir_obs.(dirpath_1B_fld));
    
    [geom_info] = get_geom_info(vir1Bdata);
    
    isinc_lt_theta = geom_info.incident_angle < theta;
    
    % Ni = nansum(isinc_lt_theta,1);
    if size(isinc_lt_theta,1) ~= vir1Bdata.hdr.lines
        fprintf('%d: image size and HKT size are inconsistent\n',sctime);
        idxBool_valid(i) = false;
        continue;
    end

    if i==1
        isinc_lt_theta_list = isinc_lt_theta;
        img_sctime_list = repmat(sctime,[vir1Bdata.hdr.lines,1]);
    else
        isinc_lt_theta_list = cat(1,isinc_lt_theta_list,isinc_lt_theta);
        img_sctime_list = cat(1,img_sctime_list,repmat(sctime,[vir1Bdata.hdr.lines,1]));
    end
    
end

idxes_valid = find(idxBool_valid);

end