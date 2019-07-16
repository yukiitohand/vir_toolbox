dir_ceres_L1B = '/Volumes/LaCie/data/vir/DWNCHVIR_l1B/DATA/20150816_HAMO/20150818_CYCLE1/';
basename_ceres = 'VIR_IR_1B_1_493158996_1';
lbl_ceres = crismlblread_v2(joinPath(dir_ceres_L1B,[basename_ceres '.LBL']));
[hdr_ceres] = extract_imghdr_from_lbl_vir(lbl_ceres);
img_ceres = envidataread(joinPath(dir_ceres_L1B,[basename_ceres '.QUB']),hdr_ceres);
img_ceres(img_ceres<-32766) = nan;

basename_ceres_cal_SS = 'DAWN_VIR_IR_SOLAR_SPECTRUM_V2';
dir_ceres_cal = '/Volumes/LaCie/data/vir/DWNCHVIR_l1B/CALIB';
lbl_ceres_cal_SS = crismlblread_v2(joinPath(dir_ceres_cal,[basename_ceres_cal_SS '.LBL']));
SStab_ceres = crismTABread(joinPath(dir_ceres_cal,[basename_ceres_cal_SS '.TAB']),lbl_ceres_cal_SS);
SS_ceres = [SStab_ceres.data.SPECTRAL_IRRADIANCE];

d_km = lbl_ceres.SPACECRAFT_SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );

img_if_ceres = img_ceres ./ reshape(SS_ceres,[1,1,length(SS_ceres)]) .*pi .* (d_au.^2);


dir_ceres_L1A = '/Volumes/LaCie/data/vir/DWNCHVIR_l1A/DATA/20150816_HAMO/20150818_CYCLE1/';
basenameA_ceres = 'VIR_IR_1A_1_493158996_1';
lblA_ceres = crismlblread_v2(joinPath(dir_ceres_L1A,[basenameA_ceres '.LBL']));
[hdrA_ceres] = extract_imghdr_from_lbl_vir(lblA_ceres);
imgA_ceres = envidataread(joinPath(dir_ceres_L1A,[basenameA_ceres '.QUB']),hdrA_ceres);
imgA_ceres(img_ceres<-32766) = nan;


%%
dir_vesta_L1B = '/Volumes/LaCie/data/vir/DWNVVIR_l1B_v2/DATA/20110811_SURVEY/20110811_CYCLE1';
basename_vesta = 'VIR_IR_1B_1_366390345_3';
lbl_vesta = crismlblread_v2(joinPath(dir_vesta_L1B,[basename_vesta '.LBL']));
[hdr_vesta] = extract_imghdr_from_lbl_vir(lbl_vesta);
img_vesta = envidataread(joinPath(dir_vesta_L1B,[basename_vesta '.QUB']),hdr_vesta);

basename_vesta_cal_SS = 'DAWN_VIR_IR_SOLAR_SPECTRUM_V2';
dir_vesta_cal = '/Volumes/LaCie/data/vir/DWNVVIR_l1B_v2/CALIB';
lbl_vesta_cal_SS = crismlblread_v2(joinPath(dir_vesta_cal,[basename_vesta_cal_SS '.LBL']));
SStab_vesta = crismTABread(joinPath(dir_vesta_cal,[basename_vesta_cal_SS '.TAB']),lbl_vesta_cal_SS);
SS_vesta = [SStab_vesta.data.SPECTRAL_IRRADIANCE];

basename_vesta_cal_ITF = 'DAWN_VIR_IR_RESP_V2';
lbl_vesta_cal_ITF = crismlblread_v2(joinPath(dir_vesta_cal,[basename_vesta_cal_ITF '.LBL']));
[hdr_vesta_cal_ITF] = extract_dathdr_from_lbl_vir(lbl_vesta_cal_ITF);
img_vesta_cal_ITF = envidataread(joinPath(dir_vesta_cal,[basename_vesta_cal_ITF '.DAT']),hdr_vesta_cal_ITF);

d_km = lbl_vesta.SPACECRAFT_SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );

img_vesta(img_vesta<-32766) = nan;

img_if_vesta = img_vesta ./ reshape(SS_vesta,[1,1,length(SS_vesta)]) .*pi .* (d_au.^2);


% img_if_vesta = img_vesta ./ reshape(SS_vesta,[1,1,length(SS_vesta)]);

%%
dir_vesta_L1A = '/Volumes/LaCie/data/vir/DWNVVIR_l1A_v2/DATA/20110811_SURVEY/20110811_CYCLE1';
basenameA_vesta = 'VIR_IR_1A_1_366390345_3';
lblA_vesta = crismlblread_v2(joinPath(dir_vesta_L1A,[basenameA_vesta '.LBL']));
[hdrA_vesta] = extract_imghdr_from_lbl_vir(lblA_vesta);
imgA_vesta = envidataread(joinPath(dir_vesta_L1A,[basenameA_vesta '.QUB']),hdrA_vesta);

basenameA_vesta_TAB = 'VIR_IR_1A_1_366390345_HK_3';
lblA_vestaTAB = crismlblread_v2(joinPath(dir_vesta_L1A,[basenameA_vesta_TAB '.LBL']));
HKTABA_vesta = crismTABread(joinPath(dir_vesta_L1A,[basenameA_vesta_TAB '.TAB']),lblA_vestaTAB);

%%
imgB_vesta_man = (imgA_vesta(2:end-1,:,:) - nanmean(imgA_vesta([1 end],:,:),1)) ./ img_vesta_cal_ITF*2;

%%
dir_vesta_L1A = '/Volumes/LaCie/data/vir/DWNVVIR_l1A_v2/DATA/20110811_SURVEY/20110811_CYCLE1';
basenameA_vesta_cal = 'VIR_IR_1A_1_366388060_3';
lblA_vesta_cal = crismlblread_v2(joinPath(dir_vesta_L1A,[basenameA_vesta_cal '.LBL']));
[hdrA_vesta_cal] = extract_imghdr_from_lbl_vir(lblA_vesta_cal);
imgA_vesta_cal = envidataread(joinPath(dir_vesta_L1A,[basenameA_vesta_cal '.QUB']),hdrA_vesta_cal);

basenameA_vesta_cal_TAB = 'VIR_IR_1A_1_366388060_HK_3';
lblA_vesta_calTAB = crismlblread_v2(joinPath(dir_vesta_L1A,[basenameA_vesta_cal_TAB '.LBL']));
HKTABA_vesta_cal = crismTABread(joinPath(dir_vesta_L1A,[basenameA_vesta_cal_TAB '.TAB']),lblA_vesta_calTAB);

%%
imgA_bias_rmvd = imgA_vesta - nanmean(imgA_vesta_cal(1:5,:,:),1);
