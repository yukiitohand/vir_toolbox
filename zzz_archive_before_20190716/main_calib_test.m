sctime = 366391692;
vir_obs_vesta = get_vir_obs_info(sctime,'dwld',0);
I1Adata_vesta = VIRdata(vir_obs_vesta.basename_ir_1A,vir_obs_vesta.dirpath_ir_1A);
I1Adata_vesta.readimg();

figure;
hold on;
for i=0:6
    s = 1+0.1*i;
    I1Ad = I1Adata_vesta.img(:,:,:) - nanmean(I1Adata_vesta.img(1,:,:),1)*s;

    plot(I1Adata_vesta.hdr.wavelength,squeeze(I1Ad([50],100,:)));
end



sctime = 498191098;
vir_obs_ceres = get_vir_obs_info(sctime,'dwld',2);
I1Adata_ceres = VIRdata(vir_obs_ceres.basename_ir_1A,vir_obs_ceres.dirpath_ir_1A);
I1Adata_ceres.readimg();

respdata_ceres = VIRdata(vir_obs_ceres.basename_cal_ir_resp,vir_obs_ceres.dirpath_cal_ir_1B);
respdata_ceres.readimg();

figure;
hold on;
for i=0:4
    s = 1+0.1*i;
    I1Ad = I1Adata_ceres.img(:,:,:) - nanmean(I1Adata_ceres.img(1,:,:),1)*s;

    plot(I1Adata_ceres.hdr.wavelength,squeeze(I1Ad([50],100,:)));
end

I1Adr = I1Ad ./ (respdata_ceres.img*0.5);





d_km = I1Adata_ceres.lbl.SPACECRAFT_SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );