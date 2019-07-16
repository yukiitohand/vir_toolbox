
opt = 'CI';

switch opt
    case 'VI'
        volumeid = 'DWNVVIR_I1B_v2'; % for vesta

        % get obs_info for Vesta measurement
        [sctime_info_vesta] = get_sctime_info(volumeid);

        % get HAMO measurement for Survey
        subdir_field = 'subdir_IR_1B';
        mission_phase_ptr = '^Survey';
        [sctime_info_vesta_mp] = get_sctimes_mission_phase(sctime_info_vesta,subdir_field,mission_phase_ptr);

        [img_spcs_meds,idxes_valid,isinc_lt_50_list] = vir_collect_median_spc(sctime_info_vesta_mp);

        save median_spectrum_vesta_Survey_I.mat sctime_info_vesta_mp idxes_valid isinc_lt_50_list img_spcs_meds

    case 'CI'
        volumeid = 'DWNCSVIR_I1B';
        % get obs_info for Ceres measurement
        [sctime_info_ceres] = get_sctime_info(volumeid);
        % get HAMO measurement for Survey
        subdir_field = 'subdir_IR_1B';
        mission_phase_ptr = '^Survey';
        [sctime_info_ceres_mp] = get_sctimes_mission_phase(sctime_info_ceres,subdir_field,mission_phase_ptr);
        [img_spcs_meds,idxes_valid,isinc_lt_50_list] = vir_collect_median_spc(sctime_info_ceres_mp);
        save median_spectrum_ceres_Survey_I.mat sctime_info_ceres_mp idxes_valid isinc_lt_50_list img_spcs_meds
        
    otherwise
        error('Undefined opt %s',opt);

end

load median_spectrum_vesta_Survey_I.mat sctime_info_vesta_mp idxes_valid isinc_lt_50_list img_spcs_meds

%%
% collect spectra
% idxBool_valid = true(length(sctime_info_vesta_mp),1);
% for i=1:length(sctime_info_vesta_mp)
%     sctime_info_idx = sctime_info_vesta_mp(i);
%     sctime = sctime_info_idx.sctime;
%     vir_obs = get_vir_obs_info(sctime,'dwld',0);
%     I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);
%     
%     [geom_info] = get_geom_info(I1Bdata);
%     
%     isinc_lt_50 = geom_info.incident_angle < 50;
%     
%     Ni = nansum(isinc_lt_50,1);
%     if size(isinc_lt_50,1) ~= I1Bdata.hdr.lines
%         fprintf('%d: image size and HKT size are inconsistent\n',sctime);
%         idxBool_valid(i) = false;
%         continue
%     end
% 
%     if i==1
%         isinc_lt_50_list = isinc_lt_50;
%     else
%         isinc_lt_50_list = cat(1,isinc_lt_50_list,isinc_lt_50);
%     end
%     
% end

% %% compute median spectra for each column
% idxes_valid = find(idxBool_valid);
% img_spcs_meds = nan(1,256,432);
% for c = 1:256
%     tic;
%     for i=1:length(idxes_valid)
%         idx = idxes_valid(i);
%         sctime_info_idx = sctime_info_vesta_mp(idx);
%         sctime = sctime_info_idx.sctime;
%         vir_obs = get_vir_obs_info(sctime,'dwld',0);
%         I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);
% 
%         I1Bimc = I1Bdata.lazyEnviReadc(c);
%         SolSpcdata= VIRdata(vir_obs.basename_cal_ir_solspec,vir_obs.dirpath_cal_ir_1B);
%         % convert radiance to I/F
%         d_km = I1Bdata.lbl.SPACECRAFT_SOLAR_DISTANCE{1};
%         [ d_au ] = km2au( d_km );
%         SolSpcdata.readTAB();
%         SS= [SolSpcdata.tab.data.SPECTRAL_IRRADIANCE];
% 
%         I1Bimc_if = I1Bimc  .* ( (pi .* (d_au.^2)) ./ SS);
% 
%         if i==1
%             I1Bimc_list = I1Bimc_if;
%         else
%             I1Bimc_list = cat(1,I1Bimc_list,I1Bimc_if);
%         end
%     end
%     I1Bimc_list_valid = I1Bimc_list(isinc_lt_50_list(:,c),:);
%     I1Bimc_list_med = nanmedian(I1Bimc_list_valid,1);
%     img_spcs_meds(1,c,:) = I1Bimc_list_med;
%     toc;
% end

%%
load median_spectrum_vesta_Survey_I.mat sctime_info_vesta_mp idxes_valid isinc_lt_50_list img_spcs_meds

%% Remove odd and even effects
% may be better to be performed after de-spiking procedure
wv = I1Bdata.hdr.wavelength;
im_med_t = squeeze(img_spcs_meds)';
[im_med_t_smthd] = vir_oddeven_rmvl_wInterp(im_med_t,wv);

%% Standard de-spiking procedure
[im_med_t_smthd_outrmvd,outliers] = standard_despiking_Carrozzo2016(im_med_t_smthd,wv);

% h = ones(3,1)/3;
% im_med_t_smthd_fil = imfilter(im_med_t_smthd,h);
% 
% r = im_med_t_smthd./im_med_t_smthd_fil;
% r_sgm = nanmean(((r-1).^2),'all');
% r_std = sqrt(r_sgm);
% 
% outliers = abs(r-1) > 3*r_std;
% 
% % interpolate outliers by a polynomial fit with 20 neighboring points.
% im_med_t_smthd_outrmvd = im_med_t_smthd;
% for c=1:256
%     for b=1:432
%         if outliers(b,c)
%             bdx_tested = max(1,(b-10)):min(432,(b+10));
%             x = wv(bdx_tested);
%             y = im_med_t_smthd(bdx_tested,c);
%             good_tmp = ~outliers(bdx_tested);
%             p = polyfit(x(good_tmp),y(good_tmp),2);
%             im_med_t_smthd_outrmvd(b,c) = polyval(p,wv(b));
%             
%         end
%     end
% end

%% Unique median spectrum is obtained
Umed = nanmedian(im_med_t_smthd_outrmvd,2);

% Maybe it is better to fit polynomials for different wavelength filter
% regions.
[~,b1420] = bisection(wv-1.42);
[~,b1560] = bisection(wv-1.56);
[~,b2410] = bisection(wv-2.41);
[~,b2610] = bisection(wv-2.61);
[~,b2590] = bisection(wv-2.59);
[~,b3773] = bisection(wv-3.773);
[~,b3740] = bisection(wv-3.74);
[~,b3830] = bisection(wv-3.83);
[~,b4350] = bisection(wv-4.35);
[~,b4540] = bisection(wv-4.54);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vesta correction over 2.61-4.07um
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.61-4.0um for Ceres
[~,b4200] = bisection(wv-4.2);
[~,b2200] = bisection(wv-2.20);
[~,b4070] = bisection(wv-4.07);
% In order to remove boundary effects, 
bdx_in = [b2200:b2410 b2610:b3740 b3830:b4200];
% bdx_in = [b2610:b3740 b3830:b4070];
x = wv(bdx_in);
y = Umed(bdx_in);
% p_Umed_coeff = polyfit(x,y,5);
p_Umed_coeff = polyfit(x,y,6);

bdx_out = [b2610:b3740 b3830:b4070];
bdx_out = [b2590:b3773 b3830:b4070];
P_Umed = polyval(p_Umed_coeff,wv(bdx_out));

figure; plot(wv,Umed); hold on; plot(wv(bdx_out),P_Umed);

%% for print
figure; plot(wv,Umed,'LineWidth',1,'Displayname','Median spectrum'); 
hold on; 
plot(wv(bdx_out),P_Umed,'LineWidth',1,'Displayname','Polynomial Fit');
xlim([2.3 4.1]);
% ylim([-0.02 0.3]);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F');
set(gca,'FontSize',14);
legend
export_fig('Figure_3.png','-transparent','-r150');


%%
%%ARTIFACT MATRIX
% question: im_med_t_smthd_outrmvd or im_med_t_smthd ? 
im_med_t_smthd_outrmvd_bdx_out = im_med_t_smthd_outrmvd(bdx_out,:);

Aout = (im_med_t_smthd_outrmvd_bdx_out-P_Umed)./P_Umed;

artifact_matrix = nan(I1Bdata.hdr.bands,I1Bdata.hdr.samples);
artifact_matrix(bdx_out,:) = Aout;

% save artifact_matrix_v2.mat artifact_matrix

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ceres correction over 4.07-5.0um
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% TEST!!
sctime = 366391692;
sctime = 498259058;
vir_obs = get_vir_obs_info(sctime,'dwld',0);
I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);
I1Bdata.readimg();

SolSpcdata= VIRdata(vir_obs.basename_cal_ir_solspec,vir_obs.dirpath_cal_ir_1B);

% convert radiance to I/F
[I1Bim_if] = vir_ra2if(I1Bdata,SolSpcdata);

c=100;
I1Bim_ifc = squeeze(I1Bim_if(:,c,:))';

wv = I1Bdata.hdr.wavelength;
[I1Bim_ifc_smthd] = vir_oddeven_rmvl_wInterp(I1Bim_ifc,wv);

I1Bim_ifc_smthd_cor = I1Bim_ifc_smthd(bdx_out,:) ./ (1+artifact_matrix(:,c));

figure; plot(wv(bdx_out),I1Bim_ifc(bdx_out,:));
figure; plot(wv(bdx_out),I1Bim_ifc_smthd(bdx_out,:));

figure; plot(wv(bdx_out),I1Bim_ifc_smthd_cor);
hold on; plot(wv(bdx_out),P_Umed,'r-','LineWidth',2);


