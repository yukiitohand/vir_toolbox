
dir_ceres_L1A = '/Volumes/LaCie/data/vir/DWNCHVIR_I1A/DATA/20150816_HAMO/20151012_CYCLE6';
basenameA_ceres = 'VIR_IR_1A_1_498259058_1';
lblA_ceres = crismlblread_v2(joinPath(dir_ceres_L1A,[basenameA_ceres '.LBL']));
[hdrA_ceres] = extract_imghdr_from_lbl_vir(lblA_ceres);
basenameA_ceres_TAB = 'VIR_IR_1A_1_498259058_HK_1';
lblA_ceresTAB = crismlblread_v2(joinPath(dir_ceres_L1A,[basenameA_ceres_TAB '.LBL']));
HKTABA_ceres = crismTABread(joinPath(dir_ceres_L1A,[basenameA_ceres_TAB '.TAB']),lblA_ceresTAB);

% sclk = HKTABA_ceres.data(1).SCET_TIME__CLOCK_;
% time_strt = lblA_ceres.START_TIME;
% time_strt_dtobj = datetime(time_strt,'format','yyyy-MM-dd''T''hh:mm:ss.SSS');

%%
cspice_furnsh( 'spice/ker' )

abcorr  = 'CN+S';
method  = 'Ellipsoid';
abnm = 'ceres'; % astronomical body name
fixref = 'CERES_FIXED'; % not sure the meaning... {'J2000','IAF'}??
scnm = 'dawn'; % spacecraft name
instnm = 'dawn_vir_ir';
NCORNR = 4;

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
n_imgfr = length(HKTABA_ceres.data);
lons_IFOV12 = zeros([n_imgfr,1]);
lats_IFOV12 = zeros([n_imgfr,1]);
lons_IFOV34 = zeros([n_imgfr,1]);
lats_IFOV34 = zeros([n_imgfr,1]);

% figure;
% hold on;

for ll=1:n_imgfr
    sclk = HKTABA_ceres.data(ll).SCET_TIME__CLOCK_;
    % Convert computer clock time to ET (seconds past % J2000, TDB).
    et = cspice_scs2e( scid, num2str(sclk) );
    lats = zeros([1,4]);lons = zeros([1,4]);
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
            plot(lon,lat,'.');
         end
    end
    lons_IFOV12(ll,1) = (lons(1)+lons(2))/2;
    lats_IFOV12(ll,1) = (lats(1)+lats(2))/2;
    lons_IFOV34(ll,1) = (lons(3)+lons(4))/2;
    lats_IFOV34(ll,1) = (lats(3)+lats(4))/2;
end

figure;
plot(lons_IFOV12,lats_IFOV12);
hold on;
plot(lons_IFOV34,lats_IFOV34);
plot([lons_IFOV12(1) lons_IFOV34(1)],[lats_IFOV12(1) lats_IFOV34(1)]);
plot([lons_IFOV12(end) lons_IFOV34(end)],[lats_IFOV12(end) lats_IFOV34(end)]);

%
% It's always good form to unload kernels after use,
% particularly in MATLAB due to data persistence.
%
%cspice_kclear

%%
% create DDR
smpls = hdrA_ceres.samples;
de_lon = zeros(n_imgfr,smpls);
de_lat = zeros(n_imgfr,smpls);

for ll=1:n_imgfr
    lons_ll = interp1([0.5 smpls+0.5],[lons_IFOV12(ll) lons_IFOV34(ll)],1:smpls);
    lats_ll = interp1([0.5 smpls+0.5],[lats_IFOV12(ll) lats_IFOV34(ll)],1:smpls);
    de_lon(ll,:) = flip(lons_ll);
    de_lat(ll,:) = flip(lats_ll);
end

%%
dir_ceres_L1B = '/Volumes/LaCie/data/vir/DWNCHVIR_I1B/DATA/20150816_HAMO/20151012_CYCLE6';
basename_ceres = 'VIR_IR_1B_1_498259058_1';
lbl_ceres = crismlblread_v2(joinPath(dir_ceres_L1B,[basename_ceres '.LBL']));
[hdr_ceres] = extract_imghdr_from_lbl_vir(lbl_ceres);
img_ceres = envidataread(joinPath(dir_ceres_L1B,[basename_ceres '.QUB']),hdr_ceres);
img_ceres(img_ceres<-32766) = nan;

status = {HKTABA_ceres.data.SHUTTER_STATUS};
status = cellfun(@(x) strip(x), status,'UniformOutput',false);

isclosed = strcmpi(status,'closed');
isopen = strcmpi(status,'open');

z = zeros(n_imgfr,smpls);
rgb = scx_rgb(img_ceres(:,:,250),0.00);
figure;
surf(de_lon(isopen,:),de_lat(isopen,:),z(isopen,:),rgb(:,:,:),'EdgeColor','none');

%%
% create GLT
range_latd = [35 60];
range_lond = [25 60];
latd0 = 50;
pixel_size = 450; % mete
place  = 'ceres_ernutet';

coslatd0 = cosd(latd0);
lat_dstep = pixel_size/(radii(3)*pi) * 180;
lon_dstep = pixel_size/(radii(1)*pi) * 180 / coslatd0;

% create GRID IMAGE
[latMAP,lonMAP,latNS,lonEW] =...
    create_grid_equirectangular(radii*10^3,latd0,range_latd,range_lond,pixel_size);

% create GLT image
dst_lmt_param = 1;
[x_glt,y_glt] = create_glt_equirectangular(de_lat(isopen,:),de_lon(isopen,:),latNS,lonEW,latd0,...
    'Dst_Lmt_Param',dst_lmt_param);

[img_proj] = img_proj_w_glt_quick(img_ceres,x_glt,y_glt);

    