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
