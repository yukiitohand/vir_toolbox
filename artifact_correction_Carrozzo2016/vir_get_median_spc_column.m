function [vir1Bimc_list_med,vir1Bimc_list] = vir_get_median_spc_column(sctime_info,idxes_valid,isinc_lt_theta_list,c,varargin)
% [vir1Bimc_list_med,vir1Bimc_list] = vir_get_median_spc_column(sctime_info,idxes_valid,isinc_lt_theta_list,c,varargin)
%  Compute median spectrum for a given column c over the valid pixels 
%  specified by isinc_lt_theta_list of the images sctime_info(idxes_valid)
%
% INPUTS
%  sctime_info: subset of LUT_SCTIME_LIST
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
%  c  : column for which median is computed.
% OUTPUTS
%  vir1Bimc_list_med: [1 x B] median spectrum
%  vir1Bimc_list: [L' x B] concatenation of all the spectra from 
%               sctime_info(idxes_valid). The actual spectra used for
%               median is obtained by I1Bimc_list(isinc_lt_theta_list(:,c),:)
% Optional Parameters
%   'SENSOR_ID'  : {'V','I'} (default) : 'I' 

sensor_id = 'I';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
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
        basename_cal_solspec_fld = 'basename_cal_ir_solspec';
        dirpath_cal_1B_fld = 'dirpath_cal_ir_1B';
    case 'V'
        error('not implemented yet');
    otherwise
        error('Invalid sensor_id %s',sensor_id);
end


for i=1:length(idxes_valid)
    idx = idxes_valid(i);
    sctime_info_idx = sctime_info(idx);
    sctime = sctime_info_idx.sctime;
    vir_obs = get_vir_obs_info(sctime,'dwld',0);
    vir1Bdata = VIRdata(vir_obs.(basename_1B_fld),vir_obs.(dirpath_1B_fld));

    vir1Bimc = vir1Bdata.lazyEnviReadc(c);
    SolSpcdata= VIRdata(vir_obs.(basename_cal_solspec_fld),vir_obs.(dirpath_cal_1B_fld));
    % convert radiance to I/F
    d_km = vir1Bdata.lbl.SPACECRAFT_SOLAR_DISTANCE{1};
    [ d_au ] = km2au( d_km );
    SolSpcdata.readTAB();
    SS= [SolSpcdata.tab.data.SPECTRAL_IRRADIANCE];

    vir1Bimc_if = vir1Bimc  .* ( (pi .* (d_au.^2)) ./ SS);

    if i==1
        vir1Bimc_list = vir1Bimc_if;
    else
        vir1Bimc_list = cat(1,vir1Bimc_list,vir1Bimc_if);
    end
end
vir1Bimc_list_valid = vir1Bimc_list(isinc_lt_theta_list(:,c),:);
vir1Bimc_list_med = nanmedian(vir1Bimc_list_valid,1);

end