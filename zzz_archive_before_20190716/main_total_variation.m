% band region
% ?1.42 µm< ? < 1.56 µm
% 2.41 µm< ? < 2.61 µm
% 3.74 µm< ? < 3.83 µm
% 4.35 µm< ? < 4.54 µm


% sctime_cal = 497915283; % for ceres
sctime_cal = 366388060;
vir_obs_cal = get_vir_obs_info(sctime_cal,'dwld',0);
I1Adata_cal = VIRdata(vir_obs_cal.basename_ir_1A,vir_obs_cal.dirpath_ir_1A);
I1Adata_cal.readimg();
bias     = nanmean(I1Adata_cal.img(1:5,:,:),1);
irlamp   = nanmean(I1Adata_cal.img(16:20,:,:),1);

irlamp_d = irlamp - nanmean(I1Adata_cal.img(11:15,:,:),1);
irlamp_d_r = irlamp_d ./ (respdata.img*0.5);

respdata_v1 = VIRdata('DAWN_VIR_IR_RESP_V1','./');
respdata_v1.readimg();
irlamp_d_r_1 = irlamp_d ./ (respdata_v1.img*0.5);

% sctime = 498191098; % for ceres
sctime = 366391692;
vir_obs = get_vir_obs_info(sctime,'dwld',0);
I1Adata = VIRdata(vir_obs.basename_ir_1A,vir_obs.dirpath_ir_1A);
I1Adata.readimg();
I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);
I1Bdata.readimg();

I1A_b = I1Adata.img - bias;


wv = I1Adata.hdr.wavelength;
% wv1 = I1Adata_vesta.hdr.wavelength;

[~,b1420] = bisection(wv-1.42);
[~,b1560] = bisection(wv-1.56);
[~,b2410] = bisection(wv-2.41);
[~,b2610] = bisection(wv-2.61);
[~,b3740] = bisection(wv-3.74);
[~,b3830] = bisection(wv-3.83);
[~,b4350] = bisection(wv-4.35);
[~,b4540] = bisection(wv-4.54);

I1Adata.loadHKTAB();
I1Adata.hkt.readTAB();
isopen = isempties(cellfun(@(x) regexpi(x,'closed','Once'),{I1Adata.hkt.tab.data.SHUTTER_STATUS},'UniformOutput',false));

c=100;
imc = squeeze(I1Adata.img(isopen,c,:))';
dark = squeeze(I1Adata.img(1,c,:));

br = b2610:b3740;
s_opt = select_stretch_scaling_param4dark(imc(br,:),dark(br));
dark_s = stretch_dark(dark(br),s_opt);

s_opt_1 = select_stretch_scaling_param4dark(imc(br,:),dark(br),'den_idx','mean');

I1A_b_d = imc(br,:)-dark(br)*s_opt;
I1A_b_d_2 = imc(br,:) - dark_s;
I1A_b_d_1 = imc(br,:)-dark(br)*s_opt_1;

figure; plot(br,I1A_b_d(:,55)); hold on; plot(br,I1A_b_d_1(:,55));
figure; plot(br,I1A_b_d(:,55)./I1A_b_d(:,114)); hold on; plot(br,I1A_b_d_1(:,55)./I1A_b_d_1(:,114));


respdata = VIRdata(vir_obs.basename_cal_ir_resp,vir_obs.dirpath_cal_ir_1B);
respdata.readimg();

I1A_b_d_r = I1A_b_d./(squeeze(respdata.img(1,c,br))*2.2);
I1A_b_d_1_r =I1A_b_d_1./(squeeze(respdata.img(1,c,br))*2.2);

d_km = I1Adata.lbl.SPACECRAFT_SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );

solspec= VIRdata(vir_obs.basename_cal_ir_solspec,vir_obs.dirpath_cal_ir_1B);
solspec.readTAB();
SS_ceres= [solspec.tab.data.SPECTRAL_IRRADIANCE]';

I1A_b_d_r_if = I1A_b_d_r ./ SS_ceres(br) .*pi .* (d_au.^2);
I1A_b_d_1_r_if = I1A_b_d_1_r ./ SS_ceres(br) .*pi .* (d_au.^2);



figure; plot(imc-dark*s_opt);

dark_s = stretch_dark(dark(br),s_opt);
I1A_b_d = imc(br)-dark_s;


tvals = nan(100,1);
sList = linspace(0.8,2.5,100);
br = b2610:b3740;
for i=1:100
    yd = imc(br,100)-dark(br)*sList(i);
    [tvals(i)] = total_variation(yd,1);
end

figure; plot(sList,tvals);

s_opts = nan(I1Adata.hdr.samples,1);
br = b2610:b3740;
for c=1:I1Adata.hdr.samples
    imc = squeeze(I1Adata.img(isopen,c,:))';
    dark = squeeze(I1Adata.img(1,c,:));
    
    [s_opts(c)] = select_stretch_scaling_param4dark(imc(br,:),dark(br));
end
figure; plot(s_opts);

d_km = I1Adata.lbl.SPACECRAFT_SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );

solspec= VIRdata(vir_obs.basename_cal_ir_solspec,vir_obs.dirpath_cal_ir_1B);
solspec.readTAB();
SS_ceres= [solspec.tab.data.SPECTRAL_IRRADIANCE]';
img_if_ceres = imc_cor ./ SS_ceres(br) .*pi .* (d_au.^2);
