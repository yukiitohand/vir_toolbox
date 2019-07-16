% correct vesta image
% convert radiance to I/F
sctime = 366393039;
% sctime = 498259058;
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
% load('spclib_relab2018Dec31.mat')
% specCode_list = {'^BC-FTIR','^BCF-FTIR',};
% generalType1_list = {'mineral','Rock','Rocks','RockCoating','MineralPwdr','Coating','biological',...
%     'organic','Chemical Compound','Regolith'};
% % generalType1_list = {'mineral','biological'};
% % spclib_relab_ftir = searchby('specCode',specCode_list,spclib_relab);
% % spclib_relab_ftir_a = searchby('generalType1',generalType1_list,spclib_relab_ftir);
% 
% specCode_list1 = {'^BD-VNIR\+BC-FTIR','^BD-VNIR\+BCF-FTIR'};
% spclib_relab_ftir_bdvnir = searchby('specCode',specCode_list1,spclib_relab);
% spclib_relab_ftir_bdvnir_a = searchby('generalType1',generalType1_list,spclib_relab_ftir_bdvnir);

% wv = I1Bdata.hdr.wavelength;
% fwhm = I1Bdata.hdr.fwhm;
% [A,option] = libstruct_convoluter(spclib_relab_ftir_bdvnir_a,wv,4,'fwhm',fwhm,'xmult',0.001,'retainRatio',0.1);
% 
% % b_end = cell2mat([spclib_relab_ftir_a.wavelength_end]);
% 
% idxb_bad = isnan(A);
% idxb_good = ~idxb_bad;
% Ahat = A;
% for i=1:size(A,2)
%     a = A(:,i);
%     if (sum(idxb_good(:,i))>1) && any(idxb_bad(:,i))
%         a_hat = interp1(wv(idxb_good(:,i)),a(idxb_good(:,i)),wv(idxb_bad(:,i)),'pchip',nan);
%         Ahat(idxb_bad(:,i),i) = a_hat;
%     end
% end

load Library_2018_ftir_bdbnir A Ahat spclib_relab_ftir_bdvnir_a

[howardite,h_i] = searchby_multfield({'type1','type2','subType'},'howardite',spclib_relab_ftir_bdvnir_a);
[Eucrite,e_i] = searchby_multfield({'type1','type2','subType'},'eucrite',spclib_relab_ftir_bdvnir_a);
[Diogenite,d_i] = searchby_multfield({'type1','type2','subType'},'diogenite',spclib_relab_ftir_bdvnir_a);
[achondrite,a_i] = searchby_multfield({'type1','type2','subType'},'achondrite',spclib_relab_ftir_bdvnir_a);
[chondrite,c_i] = searchby_multfield({'type1','type2','subType'},'chondrite',spclib_relab_ftir_bdvnir_a);
[silicate,s_i] = searchby_multfield({'type1','type2','subType'},'silicate',spclib_relab_ftir_bdvnir_a);



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
[~,b2590] = bisection(wv-2.59);
[~,b3773] = bisection(wv-3.773);
[~,b2610] = bisection(wv-2.61);
[~,b3740] = bisection(wv-3.74);
[~,b3830] = bisection(wv-3.83);
[~,b4350] = bisection(wv-4.35);
[~,b4540] = bisection(wv-4.54);
bL = length(wv);

bands = [b2590:b3773];
%bands = [b2610:b3740];
wvb = wv(bands);

Alib_b = Ahat(bands,:);
Alib_b_nrmed = normalizevec(Alib_b,1);
%valid_idx = ~any(isnan(Alib_b),1);
%Alib_b_vld = Alib_b(:,valid_idx);
BBif_b = BBif(bands,:);
BBif_b_nrmed = normalizevec(BBif_b,1);

%Alib_b_vld_nrmed = normalizevec(Alib_b_vld,1);
%[ continua,bases ] = learnContinuum(wvb,Alib_b_vld_nrmed,'convhull');

%[ Acntrmvd ] = CntRmvl(Alib_b_vld_nrmed,'additive','CONTINUA',...
%                          continua );
%Acntrmvd = -Acntrmvd;

idx_lib = [h_i e_i d_i a_i c_i];
idx_lib = unique(idx_lib);
Alib_b_nrmed_sub = Alib_b_nrmed(:,idx_lib);
valid_idx_b = ~any(isnan(Alib_b_nrmed_sub),1);
Alib_b_nrmed_sub_valid = Alib_b_nrmed_sub(:,valid_idx_b);

infoAlib = spclib_relab_ftir_bdvnir_a(idx_lib(valid_idx_b));

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


%% version 1
biasList = nan(1,256,432);
multList = nan(1,256,432);
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
        
        %idx_lib = [h_i e_i d_i a_i c_i];
        %idx_lib = unique(idx_lib);

        [X,d1,c1] = vir_unmixing_denoiser_e(Alib_b_nrmed_sub_valid(valid_bands,:),...
            BBif_b_nrmed(valid_bands,:),imgc_b(valid_bands,:),wvb(valid_bands),'C0',c0);
        biasList(1,c,bands(valid_bands)) = c1;
        multList(1,c,bands(valid_bands)) = d1;
        biasList(1,c,bands(~valid_bands)) = nan;
        multList(1,c,bands(~valid_bands)) = nan;
    else
        %idx_lib = [h_i e_i d_i a_i c_i];
        %idx_lib = unique(idx_lib);
        [X,d1,c1] = vir_unmixing_denoiser_e(Alib_b_nrmed_sub_valid,BBif_b_nrmed,imgc_b,wvb,'C0',c0);
        % [X,d1,c1] = vir_unmixing_denoiser_e(Alib_b_nrmed(:,[h_i e_i d_i a_i]),BBif_b_nrmed,imgc_b,wvb,'C0',c0);
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

%% version 3 (version 2, you need to turn off defective pixel removal)
biasList = nan(1,256,432);
multList = nan(1,256,432);
% load artifact_matrix_v2.mat artifact_matrix
% artifact_matrix = reshape(artifact_matrix',[1,256,432]);
for c = 5:256
    %%
    % c = 49;
    % bands = 180:294;
    %imgc = squeeze(I1Bim_if(:,c,:))';
    %imgc_b = imgc(bands,:);

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
%     d0 = 1+artifact_matrix(bands,c);

    valid_bands = squeeze(good_pixels(:,c,bands));

    [X,d11,c11] = vir_unmixing_denoiser_e(Alib_b_nrmed_sub_valid(valid_bands,:),BBif_b_nrmed(valid_bands,:),imgc_b1(valid_bands,:),wvb(valid_bands),'lambda_a_lib',0.01);
    biasList(1,c,bands(valid_bands)) = c11;
    multList(1,c,bands(valid_bands)) = d11;

end
a3 = load('VIR_IR_1B_1_366393039_3_cor_b2590t3773_v3.mat');
I1Bim_if_cor_v3 = (img_b1_if-a3.biasList)./a3.multList;

%% version 4
BPs = false(size(img_b1_if));
img_model = nan(size(img_b1_if));
biasList = nan(1,256,432);
multList = nan(1,256,432);

for c = 5:256
    imgc_b1 = squeeze(img_b1_if(:,c,bands))';
    valid_bands = squeeze(good_pixels(:,c,bands));
    %d0 = squeeze(vesta_parameters.multList(1,c,bands(valid_bands)));
    
    [X,d1,c1,C,Z,Xlib,XBB,bp] = vir_unmixing_denoiser_e_r(Alib_b_nrmed_sub_valid_nodup(valid_bands,:),...
        BBif_b_nrmed(valid_bands,:),imgc_b1(valid_bands,:),wvb(valid_bands),'Lambda_a_lib',0.01);
    
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

% load VIR_IR_1B_1_366393039_3_cor_b2590t3773_v4.mat

I1Bim_if_cor_v4 = (img_b1_if-biasList)./multList;
I1Bim_if_cor_v4(BPs) = img_model(BPs);

%% destriping (post processing, if needed, just tested)
img_means = nanmean(img_model,1);
img_means_fil = ones(size(img_means));
n=5;
h = ones(1,n)/n;
for b=1:size(img_means,3)
    img_mean = img_means(:,:,b);
    img_mean_fil = imfilter(img_mean,h,'replicate');
    img_means_fil(:,:,b) = img_mean_fil;
end
coeff = img_means_fil./img_means;

bb = 200;
figure; scx_rgb(img_model(:,:,bb).*coeff(:,:,bb));

I1Bim_if_cor_v4_fil = I1Bim_if_cor_v4 .* coeff;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for printing
%
%%
a2 = load('VIR_IR_1B_1_366393039_3_cor_b2590t3773_v2.mat');
I1Bim_if_cor_v2 = (img_b1_if-a2.biasList)./a2.multList;
%%

a = load('VIR_IR_1B_1_366393039_3_cor_b2590t3773.mat');
I1Bim_if_cor = (I1Bim_if-a.biasList)./a.multList;

%%
figure;
plot(wvb,squeeze(I1Bim_if(10,100,bands)),'DisplayName','I1B I/F'); 
hold on;
plot(wvb,squeeze(img_b1_if(10,100,bands)-0.1),'DisplayName','I1B I/F dark optimized'); 
%%
hold on;
plot(wvb,squeeze(I1Bim_if_cor(10,100,bands)),'DisplayName','I1B I/F corrected');
plot(wvb,squeeze(I1Bim_if_cor_v2(10,100,bands)),'DisplayName','I1B I/F corrected v2');

%% for print
c=100;
figure;
set_figsize(gcf,300,300);
plot(wvb,squeeze(I1Bim_if_cor_v2(:,c,bands))','Color',[0.7 0.7 0.7]);
hold on;
% plot(wvb,squeeze(I1Bim_if_cor_v2([7 34 63 39],c,bands))','LineWidth',1);
plot(wvb,squeeze(I1Bim_if_cor_v2([1 19 54 61],c,bands))','LineWidth',2);
xlim([2.58 3.8]);
ylim([0.25 0.44]);
set(gca,'FontSize',16);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F');
fname = sprintf('vir_%s_proposed_demo_c%03d.png',I1Bdata.basename,c);
export_fig(fname,'-transparent','-r150');
