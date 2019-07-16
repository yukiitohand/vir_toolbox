
dir_ceres_L1A = '/Volumes/LaCie/data/vir/DWNCHVIR_I1A/DATA/20150816_HAMO/20151012_CYCLE6';
basenameA_ceres = 'VIR_IR_1A_1_498259058_1';
lblA_ceres = crismlblread_v2(joinPath(dir_ceres_L1A,[basenameA_ceres '.LBL']));
basenameA_ceres_TAB = 'VIR_IR_1A_1_498259058_HK_1';
lblA_ceresTAB = crismlblread_v2(joinPath(dir_ceres_L1A,[basenameA_ceres_TAB '.LBL']));
HKTABA_ceres = crismTABread(joinPath(dir_ceres_L1A,[basenameA_ceres_TAB '.TAB']),lblA_ceresTAB);

sclk = HKTABA_ceres.data(1).SCET_TIME__CLOCK_;
time_strt = lblA_ceres.START_TIME;
time_strt_dtobj = datetime(time_strt,'format','yyyy-MM-dd''T''hh:mm:ss.SSS');

%%
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
et = cspice_scs2e( scid, num2str(sclk) );

% Convert UTC time to ET (seconds past % J2000, TDB).
% Around same as et, derived above, we use the above. This is just for
% confirmation.
et1 = cspice_str2et( time_strt );

% get the radius of an astronomical body
radii = cspice_bodvrd( target, 'RADII', 3 );

% cspice_getfov will return the name of the camera-fixed frame in the 
% string 'dref', the camera boresight vector in the array 'bsight', and the
% FOV corner vectors in the array 'bounds'.
[shape, dref, bsight, bounds] = cspice_getfov( instid, 4);


%%
obsvr = scnm;
target = abnm;
lons = [];
lats = [];
for i=1:NCORNR+1
     if( i <= NCORNR )
        fprintf( 'Corner vector %d\n\n', i)
        dvec = bounds(:,i);
     end

     if ( i == (NCORNR + 1) )
        fprintf( 'Boresight vector\n\n' )
        dvec = bsight;
     end

     %
     % Compute the surface intercept point using
     % the specified aberration corrections.
     %
     [ spoint, trgepc, srfvec, found ] =                   ...
                    cspice_sincpt( method, target,         ...
                                   et,     fixref, abcorr, ...
                                   obsvr,   dref,   dvec );
     if( found )
        %
        % Compute range from observer to apparent intercept.
        %
        dist = norm( srfvec );

        %
        % Convert rectangular coordinates to planetocentric
        % latitude and longitude. Convert radians to degrees.
        %
        [ radius, lon, lat ] = cspice_reclat( spoint );

        lon = lon * cspice_dpr;
        lat = lat * cspice_dpr;
            %
            % Display the results.
            %
            fprintf( '  Vector in %s frame = \n', dref )
            fprintf( '   %18.10e %18.10e %18.10e\n', dvec );

            fprintf( [ '\n'                                              ...
                       '  Intercept:\n'                                  ...
                       '\n'                                              ...
                       '     Radius                   (km)  = %18.10e\n' ...
                       '     Planetocentric Latitude  (deg) = %18.10e\n' ...
                       '     Planetocentric Longitude (deg) = %18.10e\n' ...
                       '     Range                    (km)  = %18.10e\n' ...
                       '\n' ],                                           ...
                        radius,  lat,  lon,  dist                          )
         else
            disp( 'Intercept not found.' )
     end
     if( i <= NCORNR )
        lons(i) = lon;
        lats(i) = lat;
     end
     

end

%
% It's always good form to unload kernels after use,
% particularly in MATLAB due to data persistence.
%
cspice_kclear

