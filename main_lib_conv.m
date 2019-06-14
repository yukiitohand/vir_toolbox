% 'mineral'
% 'biological' : AS-EAC,AS-JFM
% 'organic' : MS-CMP*
% 

load('spclib_relab2016Dec.mat')
specCode_list = {'^BC-FTIR','^BCF-FTIR'};
generalType1_list = {'mineral','Rock','Rocks','RockCoating','MineralPwdr','Coating','biological',...
    'organic','Chemical Compound','Regolith'};

spclib_relab_ftir = searchby('specCode',specCode_list,spclib_relab);


spclib_relab_ftir_a = searchby('generalType1',generalType1_list,spclib_relab_ftir);



%%
% convolve spectra
dir_ceres_L1B = '/Volumes/LaCie/data/vir/DWNCHVIR_I1B/DATA/20150816_HAMO/20151012_CYCLE6';
basenameQQ_ceres = 'VIR_IR_1B_1_498259058_1';
lbl_ceres = crismlblread_v2(joinPath(dir_ceres_L1B,[basenameQQ_ceres '.LBL']));
[hdr_ceres] = extract_imghdr_from_lbl_vir(lbl_ceres);
% img_ceres = envidataread(joinPath(dir_ceres_L1B,[basenameQQ_ceres '.QUB']),hdr_ceres);
% img_ceres(img_ceres<-32766) = nan;


wv = hdr_ceres.wavelength;
fwhm = hdr_ceres.fwhm;
[A,option] = libstruct_convoluter(spclib_relab_ftir_a,wv,2,'fwhm',fwhm,'xmult',0.001,'retainRatio',0.1);


%%
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

%%
c = 100;
bands = 180:294;
imgc = squeeze(img_if_ceres(:,100,bands))';
wvb = hdr_ceres.wavelength(bands);
wvb = wvb(:);

Ab = A(bands,:);
valid_idx = ~any(isnan(Ab),1);
Ab_vld = Ab(:,valid_idx);

[X,d,c] = vir_unmixing_denoiser(Ab_vld,imgc,wvb);