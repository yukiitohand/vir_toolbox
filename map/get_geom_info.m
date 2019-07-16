function [geom_info] = get_geom_info(virdata_obj)
% [geom_info] = get_geom_info(virdata_obj)
%  Get geometric information for the non-map projected VIR data.
%   Spice kernels are used.
%  INPUTS
%    virdata_obj: VIRdata class object
%  OUTPUTS
%    geom_info: structure storing geometrical information
%      longitude: deg
%      latitude : deg
%      incidence_angle: deg
%      emission_angle : deg
%      phase_angle    : deg


%% SET UP ENVIRONMENTAL VARIABLES
% spice_kernel_TM_name = I1Bdata.lbl.SPICE_FILE_NAME;
[~,spice_kernel_TM_base,ext] = fileparts(virdata_obj.lbl.SPICE_FILE_NAME);
spice_kernel_TM_name = [spice_kernel_TM_base '_edited' ext];
virdata_obj.loadHKTAB();
hkt = virdata_obj.hkt;
hkt.readTAB();
isopen = isempties(cellfun(@(x) regexpi(x,'closed','Once'),{hkt.tab.data.SHUTTER_STATUS},'UniformOutput',false));
hkt_tab_data_open = hkt.tab.data(isopen);

cspice_furnsh( joinPath('spice',spice_kernel_TM_name) );


abcorr  = 'CN+S';
method  = 'Ellipsoid';
switch upper(virdata_obj.lbl.TARGET_NAME)
    case '4 VESTA'
        abnm = 'vesta'; % astronomical body name
    case '1 CERES'
        abnm = 'ceres'; % astronomical body name
end
fixref = virdata_obj.lbl.COORDINATE_SYSTEM_NAME; % not sure the meaning... {'J2000','IAF'}??
scnm = 'dawn'; % spacecraft name
switch upper(virdata_obj.lbl.CHANNEL_ID)
    case 'IR'
        instnm = 'dawn_vir_ir';
    case 'VIS'
        instnm = 'dawn_vir_vis';
end
NCORNR = 4;

%%
% get spacecraft id
[scid, foundsc] = cspice_bodn2c( scnm );
if ( ~foundsc )
    txt = sprintf( [ 'SPICE(NOTRANSLATION) ' ...
                     'Could not find ID code for spacecraft %s.' ], ...
                  scnm);
    error( txt )
end

% get instrument id
[instid, foundinst] = cspice_bodn2c( instnm );
if ( ~foundinst )
    txt = sprintf( [ 'SPICE(NOTRANSLATION) ' ...
                     'Could not find ID code for instrument %s.' ], ...
                   instnm);
    error( txt )
end

% Convert computer clock time to ET (seconds past % J2000, TDB).
% et = cspice_scs2e( scid, num2str(sclk) );

% Convert UTC time to ET (seconds past % J2000, TDB).
% Around same as et, derived above, we use the above. This is just for
% confirmation.
% et1 = cspice_str2et( time_strt );

% get the radius of an astronomical body
radii = cspice_bodvrd( abnm, 'RADII', 3 );

% cspice_getfov will return the name of the camera-fixed frame in the 
% string 'dref', the camera boresight vector in the array 'bsight', and the
% FOV corner vectors in the array 'bounds'.
[shape, dref, bsight, bounds] = cspice_getfov( instid, 4);


%%
obsvr = scnm;
target = abnm;
n_imgfr = length(hkt_tab_data_open);
lons_IFOV12 = zeros([n_imgfr,1]);
lats_IFOV12 = zeros([n_imgfr,1]);
lons_IFOV34 = zeros([n_imgfr,1]);
lats_IFOV34 = zeros([n_imgfr,1]);
radius_IFOV12 = zeros([n_imgfr,1]);
radius_IFOV34 = zeros([n_imgfr,1]);

% figure;
% hold on;

for ll=1:n_imgfr
    sclk = hkt_tab_data_open(ll).SCET_TIME__CLOCK_;
    % Convert computer clock time to ET (seconds past % J2000, TDB).
    et = cspice_scs2e( scid, num2str(sclk) );
    lats = zeros([1,4]);lons = zeros([1,4]);
    radiis = zeros([1,4]);
    for i=1:NCORNR
         dvec = bounds(:,i);
         % Compute the surface intercept point using
         % the specified aberration corrections.
         [ spoint, trgepc, srfvec, found ] =                   ...
                        cspice_sincpt( method, target,         ...
                                       et,     fixref, abcorr, ...
                                       obsvr,   dref,   dvec );
         if( found )
            % Convert rectangular coordinates to planetocentric
            % latitude and longitude. Convert radians to degrees.
            [ radius, lon, lat ] = cspice_reclat( spoint );

            lon = lon * cspice_dpr;
            lat = lat * cspice_dpr;

            lons(i) = lon;
            lats(i) = lat;
            radiis(i) = radius;
            % plot(lon,lat,'.');
         end
    end
    lons_IFOV12(ll,1) = (lons(1)+lons(2))/2;
    lats_IFOV12(ll,1) = (lats(1)+lats(2))/2;
    lons_IFOV34(ll,1) = (lons(3)+lons(4))/2;
    lats_IFOV34(ll,1) = (lats(3)+lats(4))/2;
    radius_IFOV12(ll,1) = (radiis(1)+radiis(2))/2;
    radius_IFOV34(ll,1) = (radiis(3)+radiis(4))/2;
end

% figure;
% plot(lons_IFOV12,lats_IFOV12);
% hold on;
% plot(lons_IFOV34,lats_IFOV34);
% plot([lons_IFOV12(1) lons_IFOV34(1)],[lats_IFOV12(1) lats_IFOV34(1)]);
% plot([lons_IFOV12(end) lons_IFOV34(end)],[lats_IFOV12(end) lats_IFOV34(end)]);

%
% It's always good form to unload kernels after use,
% particularly in MATLAB due to data persistence.
%


%%
% create DDR
smpls = virdata_obj.hdr.samples;
de_lon = nan(n_imgfr,smpls);
de_lat = nan(n_imgfr,smpls);
de_radius = nan(n_imgfr,smpls);

for ll=1:n_imgfr
    lons_ll = interp1([0.5 smpls+0.5],[lons_IFOV12(ll) lons_IFOV34(ll)],1:smpls);
    lats_ll = interp1([0.5 smpls+0.5],[lats_IFOV12(ll) lats_IFOV34(ll)],1:smpls);
    radius_ll = interp1([0.5 smpls+0.5],[radius_IFOV12(ll) radius_IFOV34(ll)],1:smpls);
    de_lon(ll,:) = flip(lons_ll);
    de_lat(ll,:) = flip(lats_ll);
    de_radius(ll,:) = flip(radius_ll);
end

%%
% compute incidence angles
% refer https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/MATLAB/info/mostused.html#L
de_inc = nan(n_imgfr,smpls);
de_emi = nan(n_imgfr,smpls);
de_phase = nan(n_imgfr,smpls);
for ll=1:n_imgfr
    sclk = hkt_tab_data_open(ll).SCET_TIME__CLOCK_;
    et = cspice_scs2e( scid, num2str(sclk) );
    for cc = 1:smpls
        % convert planetocentric r/lon/lat to Cartesian vector
        lon = de_lon(ll,cc);
        lat = de_lat(ll,cc);
        radius = de_radius(ll,cc);
        point = cspice_latrec( radius, lon * cspice_rpd, lat * cspice_rpd );
        % compute illumination angles
        [ trgepc, srfvec, phase, solar, emissn ] = ...
         cspice_ilumin( method, target, et, fixref, abcorr, obsvr, point );
        de_inc(ll,cc) = solar;
        de_emi(ll,cc) = emissn;
    end
end

cspice_kclear

geom_info = [];
geom_info.longitude = de_lon;
geom_info.latitude  = de_lat;
geom_info.incident_angle = rad2deg(de_inc);
geom_info.emission_angle = rad2deg(de_emi);
geom_info.phase_angle    = rad2deg(de_phase);

end