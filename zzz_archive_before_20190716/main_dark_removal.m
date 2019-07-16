sctime = 366391692;
vir_obs = get_vir_obs_info(sctime,'dwld',0);
I1Adata = VIRdata(vir_obs.basename_ir_1A,vir_obs.dirpath_ir_1A);
I1Adata.readimg();
I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);
I1Bdata.readimg();

% I1A_b = I1Adata.img - bias;


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



respdata = VIRdata(vir_obs.basename_cal_ir_resp,vir_obs.dirpath_cal_ir_1B);
respdata.readimg();


SolSpcdata = VIRdata(vir_obs.basename_cal_ir_solspec,vir_obs.dirpath_cal_ir_1B);
SolSpcdata.readTAB();
SS = [SolSpcdata.tab.data.SPECTRAL_IRRADIANCE]';

I1Bim_if = vir_ra2if(I1Bdata,SolSpcdata);

%%
I1Adata.loadHKTAB();
I1Adata.hkt.readTAB();
isopen = isempties(cellfun(@(x) regexpi(x,'closed','Once'),{I1Adata.hkt.tab.data.SHUTTER_STATUS},'UniformOutput',false));

c=49;
imc = squeeze(I1Adata.img(isopen,c,:))';
dark = squeeze(I1Adata.img(1,c,:));

br = b2610:b3740;
s_opt = select_stretch_scaling_param4dark(imc(br,:),dark(br));
dark_s = stretch_dark(dark(br),s_opt);

s_opt = select_stretch_scaling_param4dark(imc(br,:),dark(br));

%%
frparam = I1Adata.lbl.FRAME_PARAMETER{1};

I1A_d = imc(br,:)-dark_s;
I1A_d_r = I1A_d./(squeeze(respdata.img(1,c,br))*frparam(1));

d_km = I1Adata.lbl.SPACECRAFT_SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );
SS= [SolSpcdata.tab.data.SPECTRAL_IRRADIANCE]';

I1A_d_r_if = I1A_d_r  .* (pi .* (d_au.^2) ./ SS(br));

%%
figure; plot(wv(br),squeeze(I1Bim_if(:,c,br))');
figure; plot(wv(br),squeeze(I1A_d_r_if));


