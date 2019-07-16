vir_obs = get_vir_obs_info(366886513,'dwld',0);
I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);
SolSpcdata= VIRdata(vir_obs.basename_cal_ir_solspec,vir_obs.dirpath_cal_ir_1B);
I1Bimg_if = vir_ra2if(I1Bdata,SolSpcdata);
I1Bdata.img = I1Bimg_if;
% hsiview_v3(I1Bdata.img(:,:,200),{I1Bdata},{'i1b vesta'});

c=178;
I1Bim_ifc = squeeze(I1Bimg_if(:,c,:))';

wv = I1Bdata.hdr.wavelength;
[I1Bim_ifc_smthd_vesta] = vir_oddeven_rmvl_wInterp(I1Bim_ifc,wv);
figure; plot(wv,I1Bim_ifc_smthd_vesta(:,22));

sctime = 486824235;
sctime = 498191098;
% vir_obs = get_vir_obs_info(498191098,'dwld',0);
vir_obs = get_vir_obs_info(sctime,'dwld',0);
I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);
SolSpcdata = VIRdata(vir_obs.basename_cal_ir_solspec,vir_obs.dirpath_cal_ir_1B);
I1Bimg_if = vir_ra2if(I1Bdata,SolSpcdata);
I1Bdata.img = I1Bimg_if;
hsiview_v3(I1Bdata.img(:,:,200),{I1Bdata},{'i1b ceres'});

c=100;
I1Bim_ifc = squeeze(I1Bimg_if(:,c,:))';

wv = I1Bdata.hdr.wavelength;
[I1Bim_ifc_smthd_ceres] = vir_oddeven_rmvl_wInterp(I1Bim_ifc,wv);
figure; plot(wv,I1Bim_ifc_smthd_ceres(:,47));

figure; plot(wv,I1Bim_ifc_smthd_ceres(:,47)./I1Bim_ifc_smthd_vesta(:,47));



%% check calibration
sctime = 366388060;
vir_obs = get_vir_obs_info(sctime,'dwld',0);
I1Adata = VIRdata(vir_obs.basename_ir_1A,vir_obs.dirpath_ir_1A);
I1Adata.readimg();
respdata = VIRdata(vir_obs.basename_cal_ir_resp,vir_obs.dirpath_cal_ir_1B);
respdata.readimg();

I1Aim_irlamp_d = I1Adata.img(16:20,:,:) - nanmean(I1Adata.img(6:10,:,:),1);
I1Aim_irlamp_d_ra = I1Aim_irlamp_d ./ respdata.img;

% sctime_c = 497915283;
%sctime_c = 493156585;
sctime_c = 486792974;
%sctime_c = 503528899;
vir_obs_c = get_vir_obs_info(sctime_c,'dwld',0);
I1Adata_c = VIRdata(vir_obs_c.basename_ir_1A,vir_obs_c.dirpath_ir_1A);
I1Adata_c.readimg();
respdata = VIRdata(vir_obs.basename_cal_ir_resp,vir_obs.dirpath_cal_ir_1B);
respdata.readimg();

I1Aim_c_irlamp_d = I1Adata_c.img(16:20,:,:) - nanmean(I1Adata_c.img(6:10,:,:),1);
I1Aim_c_irlamp_d_ra = I1Aim_c_irlamp_d ./ respdata.img;

% figure; plot(I1Adata.hdr.wavelength,squeeze(I1Aim_irlamp_d(3,178,:)));
% hold on;
% plot(I1Adata_c.hdr.wavelength,squeeze(I1Aim_c_irlamp_d(3,178,:)));

figure; plot(I1Adata.hdr.wavelength,squeeze(I1Aim_irlamp_d_ra(3,178,:)));
hold on;
plot(I1Adata_c.hdr.wavelength,squeeze(I1Aim_c_irlamp_d_ra(3,178,:)));


%%
% searching calibration
volumeid_ch = 'DWNCHVIR_I1A';
[sctime_info_ch] = get_sctime_info(volumeid_ch);

for i=1:length(sctime_info_ch)
    sctime = sctime_info_ch(i).sctime;
    vir_obs = get_vir_obs_info(sctime,'dwld',0);
    I1Adata = VIRdata(vir_obs.basename_ir_1A,vir_obs.dirpath_ir_1A);
    if strcmpi(I1Adata.lbl.TARGET_TYPE,'calibration')
        fprintf('%s: %d \n',vir_obs.subdir_ir_1A,sctime);
    end
end
% DATA/20150816_HAMO/20150818_CYCLE1: 493156585 
% DATA/20150816_HAMO/20150829_CYCLE2: 494111609 
% DATA/20150816_HAMO/20150909_CYCLE3: 495061440 
% DATA/20150816_HAMO/20150920_CYCLE4: 496015533 
% DATA/20150816_HAMO/20151001_CYCLE5: 496965137 
% DATA/20150816_HAMO/20151012_CYCLE6: 497915283 
% it seems that only one calibration is performed per cycle.

%% black body radiation
sctime_c = 498259058;
sctime_c = 493158996;
vir_obs_c = get_vir_obs_info(sctime_c,'dwld',2);
I1Bdata = VIRdata(vir_obs_c.basename_ir_1B,vir_obs_c.dirpath_ir_1B);
I1Bdata.readimg();
V1Bdata = VIRdata(vir_obs_c.basename_vis_1B,vir_obs_c.dirpath_vis_1B);
V1Bdata.readimg();

SolSpcdata = VIRdata(vir_obs_c.basename_cal_ir_solspec,vir_obs_c.dirpath_cal_ir_1B);
I1Bimg_if = vir_ra2if(I1Bdata,SolSpcdata);

% I1Bdata.loadHKTAB();
% I1Bdata.hkt.readTAB();
geom_I1Bdata = get_geom_info(I1Bdata);
% status = {I1Bdata.hkt.tab.data.SHUTTER_STATUS};
% status = cellfun(@(x) strip(x), status,'UniformOutput',false);
% isclosed = strcmpi(status,'closed');
% isopen = strcmpi(status,'open');

z = zeros(I1Bdata.hdr.lines,I1Bdata.hdr.samples);
rgb = scx_rgb(I1Bimg_if(:,:,250),0.01);
figure;
surf(geom_I1Bdata.longitude,geom_I1Bdata.latitude,z,rgb,'EdgeColor','none');

c = 178;
I1Bdata_c = squeeze(I1Bdata.img(:,c,:))';
figure; plot(I1Bdata.hdr.wavelength,I1Bdata_c);

TList = 90:5:250
figure; hold on;
for i=1:length(TList)
    T = TList(i);
    blck_spc = get_black_body_radiation(I1Bdata.hdr.wavelength,T);
    plot(I1Bdata.hdr.wavelength,blck_spc*10^-6,'DisplayName',['T=' num2str(T) 'K']);
end

c = 100;
I1Bimg_ifc = squeeze(I1Bimg_if(:,c,:))';

figure; plot(I1Bimg_ifc-blck_spc.*10^(-6));