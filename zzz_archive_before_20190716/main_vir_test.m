dir_ceres_L1B = '/Volumes/LaCie/data/vir/DWNCHVIR_I1B/DATA/20150816_HAMO/20151012_CYCLE6';
basename_ceres = 'VIR_IR_1B_1_498259058_1';
lbl_ceres = crismlblread_v2(joinPath(dir_ceres_L1B,[basename_ceres '.LBL']));
[hdr_ceres] = extract_imghdr_from_lbl_vir(lbl_ceres);
img_ceres = envidataread(joinPath(dir_ceres_L1B,[basename_ceres '.QUB']),hdr_ceres);
img_ceres(img_ceres<-32766) = nan;

basename_ceres_cal_SS = 'DAWN_VIR_IR_SOLAR_SPECTRUM_V2';
dir_ceres_cal = '/Volumes/LaCie/data/vir/DWNCHVIR_I1B/CALIB';
lbl_ceres_cal_SS = crismlblread_v2(joinPath(dir_ceres_cal,[basename_ceres_cal_SS '.LBL']));
SStab_ceres = crismTABread(joinPath(dir_ceres_cal,[basename_ceres_cal_SS '.TAB']),lbl_ceres_cal_SS);
SS_ceres = [SStab_ceres.data.SPECTRAL_IRRADIANCE];

d_km = lbl_ceres.SPACECRAFT_SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );

img_if_ceres = img_ceres ./ reshape(SS_ceres,[1,1,length(SS_ceres)]) .*pi .* (d_au.^2);

basenameB_ceres_TAB = 'VIR_IR_1B_1_498259058_HK_1';
lblB_ceresTAB = crismlblread_v2(joinPath(dir_ceres_L1B,[basenameB_ceres_TAB '.LBL']));
HKTABB_ceres = crismTABread(joinPath(dir_ceres_L1B,[basenameB_ceres_TAB '.TAB']),lblB_ceresTAB);

dir_ceres_L1A = '/Volumes/LaCie/data/vir/DWNCHVIR_I1A/DATA/20150816_HAMO/20151012_CYCLE6';
basenameA_ceres = 'VIR_IR_1A_1_498259058_1';
lblA_ceres = crismlblread_v2(joinPath(dir_ceres_L1A,[basenameA_ceres '.LBL']));
[hdrA_ceres] = extract_imghdr_from_lbl_vir(lblA_ceres);
imgA_ceres = envidataread(joinPath(dir_ceres_L1A,[basenameA_ceres '.QUB']),hdrA_ceres);
imgA_ceres(img_ceres<-32766) = nan;

basenameA_ceres_TAB = 'VIR_IR_1A_1_498259058_HK_1';
lblA_ceresTAB = crismlblread_v2(joinPath(dir_ceres_L1A,[basenameA_ceres_TAB '.LBL']));
HKTABA_ceres = crismTABread(joinPath(dir_ceres_L1A,[basenameA_ceres_TAB '.TAB']),lblA_ceresTAB);

%%
% updated version.
sctime = 516144471; % 1sec exposure time
vir_obs_ceres = get_vir_obs_info(sctime,'dwld',0);
% vir_data_1B = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);
% vir_data_1B.readimg();
vir_data_1A_ceres = VIRdata(vir_obs_ceres.basename_ir_1A,vir_obs_ceres.dirpath_ir_1A);
vir_data_1A_ceres.readimg();
vir_data_1B_ceres = VIRdata(vir_obs_ceres.basename_ir_1B,vir_obs_ceres.dirpath_ir_1B);
vir_data_1B_ceres.readimg();
respdata_ceres = VIRdata(vir_obs_ceres.basename_cal_ir_resp,vir_obs_ceres.dirpath_cal_ir_1B);
respdata_ceres.readimg();
solspec_ceres = VIRdata(vir_obs_ceres.basename_cal_ir_solspec,vir_obs_ceres.dirpath_cal_ir_1B);
solspec_ceres.readTAB();
SS_ceres = [solspec_ceres.tab.data.SPECTRAL_IRRADIANCE];

imgB_ceres_man = (vir_data_1A_ceres.img(2:end-1,:,:) - nanmean(vir_data_1A_ceres.img([1 end],:,:),1)) ./ (respdata_ceres.img*(2.2));

d_km = vir_data_1B_ceres.lbl.SPACECRAFT_SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );

img_if_ceres_man = imgB_ceres_man ./ reshape(SS_ceres,[1,1,length(SS_ceres)]) .*pi .* (d_au.^2);
img_if_ceres = vir_data_1B_ceres.img ./ reshape(SS_ceres,[1,1,length(SS_ceres)]) .*pi .* (d_au.^2);

%%
sctime = 371373500;
vir_obs_vesta = get_vir_obs_info(sctime,'dwld',2);
respdata_vesta = VIRdata(vir_obs_vesta.basename_cal_ir_resp,vir_obs_vesta.dirpath_cal_ir_1B);


%% compare calibration images
% updated version.
sctime = 370613575; % band 200 is dead
sctime = 371586110;
sctime = 371589712;
sctime = 366390345;
% sctime = 371191309;
sctime = 366390345;
sctime = 366388060;
sctime = 360823166;
sctime = 367193249; % 1sec exposure time
vir_obs_vesta = get_vir_obs_info(sctime,'dwld',0);
vir_data_1A_vesta = VIRdata(vir_obs_vesta.basename_ir_1A,vir_obs_vesta.dirpath_ir_1A);
vir_data_1A_vesta.readimg();

vir_data_1B_vesta = VIRdata(vir_obs_vesta.basename_ir_1B,vir_obs_vesta.dirpath_ir_1B);
vir_data_1B_vesta.readimg();
d_km = vir_data_1B_vesta.lbl.SPACECRAFT_SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );

solspec_vesta = VIRdata(vir_obs_vesta.basename_cal_ir_solspec,vir_obs_vesta.dirpath_cal_ir_1B);
solspec_vesta.readTAB();
SS_vesta = [solspec_vesta.tab.data.SPECTRAL_IRRADIANCE];
img_if_vesta = vir_data_1B_vesta.img ./ reshape(SS_vesta,[1,1,length(SS_vesta)]) .*pi .* (d_au.^2);


sctime = 493156585;
vir_obs_ceres = get_vir_obs_info(sctime,'dwld',2);
vir_data_1A_ceres = VIRdata(vir_obs_ceres.basename_ir_1A,vir_obs_ceres.dirpath_ir_1A);
vir_data_1A_ceres.readimg();

