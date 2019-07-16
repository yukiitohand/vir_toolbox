sctime = 498259058;
% sctime = 494253260;
vir_obs = get_vir_obs_info(sctime,'dwld',0);
I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);
I1Adata = VIRdata(vir_obs.basename_ir_1A,vir_obs.dirpath_ir_1A);
I1Adata.readimg();
frame_parms = I1Adata.lbl.FRAME_PARAMETER{1};
ir_expo = frame_parms(1);

d_km = I1Adata.lbl.SPACECRAFT_SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );

I1Adata.loadHKTAB();
I1Adata.hkt.readTAB();
isopen = isempties(cellfun(@(x) regexpi(x,'closed','Once'),{I1Adata.hkt.tab.data.SHUTTER_STATUS},'UniformOutput',false));

respdata = VIRdata(vir_obs.basename_cal_ir_resp,vir_obs.dirpath_cal_ir_1B);
respdata.readimg();

SolSpcdata = VIRdata(vir_obs.basename_cal_ir_solspec,vir_obs.dirpath_cal_ir_1B);
SolSpcdata.readTAB();
SS = [SolSpcdata.tab.data.SPECTRAL_IRRADIANCE]';

I1Bim_if = vir_ra2if(I1Bdata,SolSpcdata);

%%
load Library_2018_ftir_bdbnir A Ahat spclib_relab_ftir_bdvnir_a

[howardite,h_i] = searchby_multfield({'type1','type2','subType'},'howardite',spclib_relab_ftir_bdvnir_a);
[Eucrite,e_i] = searchby_multfield({'type1','type2','subType'},'eucrite',spclib_relab_ftir_bdvnir_a);
[Diogenite,d_i] = searchby_multfield({'type1','type2','subType'},'diogenite',spclib_relab_ftir_bdvnir_a);
[achondrite,a_i] = searchby_multfield({'type1','type2','subType'},'achondrite',spclib_relab_ftir_bdvnir_a);
[chondrite,c_i] = searchby_multfield({'type1','type2','subType'},'chondrite',spclib_relab_ftir_bdvnir_a);
[silicate,s_i] = searchby_multfield({'type1','type2','subType'},'silicate',spclib_relab_ftir_bdvnir_a);

[carbonate,carb_i] = searchby_multfield({'type1','type2','subType'},'carbonate',spclib_relab_ftir_bdvnir_a);
[sulfate,sulf_i] = searchby_multfield({'type1','type2','subType'},'sulfate',spclib_relab_ftir_bdvnir_a);
[oxide,oxi_i] = searchby_multfield({'type1','type2','subType'},'oxide',spclib_relab_ftir_bdvnir_a);

[kerite,k_i] = searchby('sampleID','MA-ATB-043',spclib_relab_ftir_bdvnir_a);
[asphaltite1,a1_i] = searchby('sampleID','MS-CMP-',spclib_relab_ftir_bdvnir_a);
[asphaltite2,a2_i] = searchby('sampleID','AS-LXM-',spclib_relab_ftir_bdvnir_a);

%% create black body emission spectra
TList = 130:5:255;
for i=1:length(TList)
    T = TList(i);
    blck_spc = get_black_body_radiation(I1Bdata.hdr.wavelength,T);
    if i==1
        BB = blck_spc;
    else
        BB = cat(2,BB,blck_spc);
    end
end

% multiply with sun flux (After I/F conversion, the contribution of
% emission is also multiplied.
d_km = I1Bdata.lbl.SPACECRAFT_SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );
BBif = BB .* ( (pi .* (d_au.^2)) ./ SS) * 10^(-6);

wv = I1Bdata.hdr.wavelength;
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
bL = length(wv);

bands = [b2610:b3740];
bands = [b2590:b3773];
wvb = wv(bands);

Alib_b = Ahat(bands,:);
Alib_b_nrmed = normalizevec(Alib_b,1);
% valid_idx = ~any(isnan(Alib_b),1);
% Alib_b_vld = Alib_b(:,valid_idx);
BBif_b = BBif(bands,:);
BBif_b_nrmed = normalizevec(BBif_b,1);

idx_lib = [h_i e_i d_i a_i c_i s_i a1_i a2_i k_i carb_i sulf_i oxi_i];
idx_lib = unique(idx_lib);
Alib_b_nrmed_sub = Alib_b_nrmed(:,idx_lib);
valid_idx_b = ~any(isnan(Alib_b_nrmed_sub),1);
Alib_b_nrmed_sub_valid = Alib_b_nrmed_sub(:,valid_idx_b);

infoAlib = spclib_relab_ftir_bdvnir_a(idx_lib(valid_idx_b));

% some spectra seems to be same
% 'BKR1CB015A' and 'BKR2CB015A'
% 'BKR1CB018A' and 'BKR2CB018A'
% 'BKR1CB043A' and 'BKR2CB043A'
% search duplicate entries
XX = Alib_b_nrmed_sub_valid'*Alib_b_nrmed_sub_valid;
issame_XX = (XX>(1-1e-8));
dup_idxes = false(1,size(XX,2));
for i=1:(size(XX,1)-1)
    idx = (i+1):size(XX,1);
    dup_idxes(idx) = issame_XX(i,idx);
end

no_duplicates = ~dup_idxes;

infoAlib_nodup = infoAlib(no_duplicates);
Alib_b_nrmed_sub_valid_nodup = Alib_b_nrmed_sub_valid(:,no_duplicates);


% Alib_b_vld_nrmed = normalizevec(Alib_b_vld,1);
% [ continua,bases ] = learnContinuum(wvb,Alib_b_vld_nrmed,'convhull');
% 
% [ Acntrmvd ] = CntRmvl(Alib_b_vld_nrmed,'additive','CONTINUA',...
%                           continua );
% Acntrmvd = -Acntrmvd;

%%
vesta_parameters = load('VIR_IR_1B_1_366393039_3_cor_b2590t3773_v2.mat');
biasList = nan(1,256,432);
multList = nan(1,256,432);

%% version 1
for c = 5:256
    %%
%     c = 49;
    % bands = 180:294;
    imgc = squeeze(I1Bim_if(:,c,:))';
    imgc_b = imgc(bands,:);
    
    % initialization parameters
    %load artifact_matrix.mat artifact_matrix
    %d0 = 1+artifact_matrix(bands,c);
    imgc_mean = nanmean(imgc,2);
    [imgc_mean_smthd] = vir_oddeven_rmvl_wInterp(imgc_mean,wv);
    imgc_mean_b = imgc_mean(bands);
    c00 = imgc_mean - imgc_mean_smthd;
    c0 = c00(bands);
    
    d0 = squeeze(vesta_parameters.multList(1,c,bands));
    
    if any(imgc_mean_b<0.01) % exception when we have bad detectors. should be deal with more smartly.
        imgc = squeeze(I1Bim_if(:,c,:))';
        imgc_b = imgc(bands,:);

        % initialization parameters
        %load artifact_matrix.mat artifact_matrix
        %d0 = 1+artifact_matrix(bands,c);
        imgc_mean = nanmean(imgc,2);
        imgc_mean_b = imgc_mean(bands);
        valid_bands = imgc_mean_b>0.01;
        [imgc_mean_smthd] = vir_oddeven_rmvl_wInterp(imgc_mean_b(valid_bands),wvb(valid_bands));
        c0 = imgc_mean_b(valid_bands) - imgc_mean_smthd;
        c0(1) = 0;
        c0(end) = 0;
        
        Alib_b_nrmed_sub = Alib_b_nrmed(:,[h_i e_i d_i a_i c_i s_i a1_i a2_i k_i carb_i sulf_i oxi_i]);
        valid_idx_b = ~any(isnan(Alib_b_nrmed_sub),1);
        [X,d1,c1] = vir_unmixing_denoiser_e(Alib_b_nrmed_sub(:,valid_idx_b),...
            BBif_b_nrmed(valid_bands,:),imgc_b(valid_bands,:),wvb(valid_bands),'C0',c0,'D0',d0(valid_bands),'Lambda_a_lib',0.01);
        biasList(1,c,bands(valid_bands)) = c1;
        multList(1,c,bands(valid_bands)) = d1;
        biasList(1,c,bands(~valid_bands)) = nan;
        multList(1,c,bands(~valid_bands)) = nan;
    else
        idx_lib = [h_i e_i d_i a_i c_i s_i a1_i a2_i k_i carb_i sulf_i oxi_i];
        idx_lib = unique(idx_lib);
        Alib_b_nrmed_sub = Alib_b_nrmed(:,idx_lib);
        valid_idx_b = ~any(isnan(Alib_b_nrmed_sub),1);
        % d0 = squeeze(vesta_parameters.multList(1,47,bands));
        [X,d1,c1,C,Z,Xlib,XBB] = vir_unmixing_denoiser_e(Alib_b_nrmed_sub(:,valid_idx_b),BBif_b_nrmed,imgc_b,wvb,'C0',c0,'D0',d0,'Lambda_a_lib',0.1);
        biasList(1,c,bands) = c1;
        multList(1,c,bands) = d1;
    end
        
    
   
end

I1Bim_if_cor = (I1Bim_if-biasList)./multList;

%% Improved dark removal on I1A data.
img_b1_if = nan(sum(isopen),256,432);
img_b1_dark = nan(1,256,432);
img_b1_d  = nan(sum(isopen),256,432);
[defective_pixels] = vir_defective_pixels();
good_pixels = (defective_pixels==0);
good_pixels_1nan = convertBoolTo1nan(good_pixels);
img_A = I1Adata.img.*good_pixels_1nan; % only consider good pixels.
for c = 5:256
    %%
    %c = 66;
    % bands = 180:294;
    %imgc = squeeze(I1Bim_if(:,c,:))';
    %imgc_b = imgc(bands,:);
    % covert I1A to I/F
    imgc_A = squeeze(img_A(isopen,c,:))';
    dark = squeeze(img_A(1,c,:));
    valid_bands = squeeze(good_pixels(:,c,bands));
    s_opt = select_stretch_scaling_param4dark(imgc_A(bands(valid_bands),:),dark(bands(valid_bands)));
    dark_s = stretch_dark(dark(bands),s_opt);
    
    img_b1_dark(1,c,bands) = reshape(dark_s',[1,1,length(bands)]);

    imgc_A_d = imgc_A(bands,:) - dark_s;
    imgc_A_d_r = imgc_A_d./(squeeze(respdata.img(1,c,bands))*ir_expo);
    imgc_b1 = imgc_A_d_r ./ SS(bands) .*pi .* (d_au.^2);
    
    img_b1_d(:,c,bands) = reshape(imgc_A_d',[sum(isopen),1,length(bands)]);
    img_b1_if(:,c,bands) = reshape(imgc_b1',[sum(isopen),1,length(bands)]);
    
%     figure(1); plot(bands,imgc_A_d);
%     title(num2str(c));
%     pause;

end
%I1Bim_if_cor_v2 = (img_b1_if-biasList)./multList;

%%
vesta_parameters = load('VIR_IR_1B_1_366393039_3_cor_b2590t3773_v4.mat');
biasList = nan(1,256,432);
multList = nan(1,256,432);

%% version 3 (version 2, you need to turn off defective pixel removal)
for c = 5:256
    %%
  %  c=50;
  % c = 49;
% bands = 180:294;

% covert I1A to I/F
%     imgc_A = squeeze(I1Adata.img(isopen,c,:))';
%     dark = squeeze(I1Adata.img(1,c,:));
%     s_opt = select_stretch_scaling_param4dark(imgc_A(bands,:),dark(bands));
%     dark_s = stretch_dark(dark(bands),s_opt);
% 
%     imgc_A_d = imgc_A(bands,:) - dark_s;
%     imgc_A_d_r = imgc_A_d./(squeeze(respdata.img(1,c,bands))*ir_expo);
%     imgc_b1 = imgc_A_d_r ./ SS(bands) .*pi .* (d_au.^2);
    
    imgc_b1 = squeeze(img_b1_if(:,c,bands))';
    valid_bands = squeeze(good_pixels(:,c,bands));
    d0 = squeeze(vesta_parameters.multList(1,c,bands(valid_bands)));

    [X,d1,c1,C,Z,Xlib,XBB] = vir_unmixing_denoiser_e(Alib_b_nrmed_sub_valid(valid_bands,:),...
        BBif_b_nrmed(valid_bands,:),imgc_b1(valid_bands,:),wvb(valid_bands),'D0',d0,'Lambda_a_lib',0.1);
    
    %%
    biasList(1,c,bands(valid_bands)) = c1;
    multList(1,c,bands(valid_bands)) = d1;
end

%% version 4
BPs = false(size(img_b1_if));
img_model = nan(size(img_b1_if));

for c = 5:256
    %%
  %  c=50;
  % c = 49;
% bands = 180:294;

% covert I1A to I/F
%     imgc_A = squeeze(I1Adata.img(isopen,c,:))';
%     dark = squeeze(I1Adata.img(1,c,:));
%     s_opt = select_stretch_scaling_param4dark(imgc_A(bands,:),dark(bands));
%     dark_s = stretch_dark(dark(bands),s_opt);
% 
%     imgc_A_d = imgc_A(bands,:) - dark_s;
%     imgc_A_d_r = imgc_A_d./(squeeze(respdata.img(1,c,bands))*ir_expo);
%     imgc_b1 = imgc_A_d_r ./ SS(bands) .*pi .* (d_au.^2);
    
    imgc_b1 = squeeze(img_b1_if(:,c,bands))';
    valid_bands = squeeze(good_pixels(:,c,bands));
    d0 = squeeze(vesta_parameters.multList(1,c,bands(valid_bands)));

    [X,d1,c1,C,Z,Xlib,XBB,bp] = vir_unmixing_denoiser_e_r_ceres(Alib_b_nrmed_sub_valid_nodup(valid_bands,:),...
        BBif_b_nrmed(valid_bands,:),imgc_b1(valid_bands,[1:16,21:end]),wvb(valid_bands),'D0',d0,'Lambda_a_lib',0.01);
    
    CZ = nan(size(imgc_b1));
    CZ(valid_bands,:) = C*Z;
    Bg = interp_nan_column(CZ,isnan(CZ),wvb);
    
    Ym = Alib_b_nrmed_sub_valid_nodup*Xlib + BBif_b_nrmed*XBB + Bg;
    
    %%
    biasList(1,c,bands(valid_bands)) = c1;
    multList(1,c,bands(valid_bands)) = d1;
    BPs(:,c,bands(valid_bands)) = reshape(bp',[size(imgc_b1,2),1,sum(valid_bands)]);
    img_model(:,c,bands) = reshape(Ym',[size(imgc_b1,2),1,length(bands)]);
end

%%
% load VIR_IR_1B_1_498259058_1_cor_b2590t3773_v4.mat

I1Bim_if_cor_v4 = (img_b1_if-biasList)./multList;

% postprocessing (filling bad pixels)
I1Bim_if_cor_v4_bands = I1Bim_if_cor_v4(:,:,bands);
img_model_bands = img_model(:,:,bands);
isnan_bands = isnan(I1Bim_if_cor_v4_bands);
I1Bim_if_cor_v4_bands(isnan_bands) = img_model_bands(isnan_bands);
I1Bim_if_cor_v4(:,:,bands) = I1Bim_if_cor_v4_bands;
I1Bim_if_cor_v4(BPs) = img_model(BPs);



%%
% a3 = load('VIR_IR_1B_1_366393039_3_cor_b2590t3773_v3.mat');
% I1Bim_if_cor_v3 = (img_b1_if-a3.biasList)./a3.multList;
        
    
%%
c=50;
figure;
set_figsize(gcf,367,420);
plot(wvb,squeeze(I1Bim_if_cor_v4(:,c,bands))');
xlim([2.58 3.8]);
ylim([-0.0164 0.0836]);
set(gca,'FontSize',12);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F');
fname = sprintf('vir_%s_proposed_c%03d_v2.png',I1Bdata.basename,c);
export_fig(fname,'-transparent','-r150');
