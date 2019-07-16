% 'mineral'
% 'biological' : AS-EAC,AS-JFM
% 'organic' : MS-CMP*
% 

load('spclib_relab2018Dec31.mat')
specCode_list = {'^BC-FTIR','^BCF-FTIR',};
generalType1_list = {'mineral','Rock','Rocks','RockCoating','MineralPwdr','Coating','biological',...
    'organic','Chemical Compound','Regolith'};
% generalType1_list = {'mineral','biological'};
spclib_relab_ftir = searchby('specCode',specCode_list,spclib_relab);
spclib_relab_ftir_a = searchby('generalType1',generalType1_list,spclib_relab_ftir);

specCode_list1 = {'^BD-VNIR\+BC-FTIR','^BD-VNIR\+BCF-FTIR'};
spclib_relab_ftir_bdvnir = searchby('specCode',specCode_list1,spclib_relab);
spclib_relab_ftir_bdvnir_a = searchby('generalType1',generalType1_list,spclib_relab_ftir_bdvnir);

% RELAB samples in "?Localized aliphatic organic material on the surface of
% Ceres", Science 2017
% MA-ATB-043  Kerite
% MS-CMP-030  Asphaltite-1
% AS-LXM-004  Asphaltite-2


[kerite,k_i] = searchby('sampleID','MA-ATB-043',spclib_relab_ftir_bdvnir);

[asphaltite1,a1_i] = searchby('sampleID','MS-CMP-',spclib_relab_ftir_bdvnir_a);
[asphaltite2,a2_i] = searchby('sampleID','AS-LXM-',spclib_relab_ftir_a);

[A,option] = libstruct_convoluter(asphaltite1,wv,4,'fwhm',fwhm,'xmult',0.001,'retainRatio',0.1);

% b = cell2mat([spclib_relab_ftir_a.wavelength_end]);




%%
% convolve spectra


wv = I1Bdata.hdr.wavelength;
fwhm = I1Bdata.hdr.fwhm;
[A,option] = libstruct_convoluter(spclib_relab_ftir_bdvnir_a,wv,4,'fwhm',fwhm,'xmult',0.001,'retainRatio',0.1);

% b_end = cell2mat([spclib_relab_ftir_a.wavelength_end]);

idxb_bad = isnan(A);
idxb_good = ~idxb_bad;
Ahat = A;
for i=1:size(A,2)
    a = A(:,i);
    if (sum(idxb_good(:,i))>1) && any(idxb_bad(:,i))
        a_hat = interp1(wv(idxb_good(:,i)),a(idxb_good(:,i)),wv(idxb_bad(:,i)),'pchip',nan);
        Ahat(idxb_bad(:,i),i) = a_hat;
    end
end


% plot spectra
figure(1); 
splib_view = pyroxene;
ii = p_i;
for i=1:length(splib_view)
    cla;
    jj = ii(i);
    plot(wv,Ahat(:,jj),'.-');
    title(sprintf('%04d:%s (%s:%s)',jj,splib_view(i).subType,...
        splib_view(i).sampleID,splib_view(i).spectrumID));
    pause;
end

[howardite,h_i] = searchby_multfield({'type1','type2','subType'},'howardite',spclib_relab_ftir_bdvnir_a);
[Eucrite,e_i] = searchby_multfield({'type1','type2','subType'},'eucrite',spclib_relab_ftir_bdvnir_a);
[Diogenite,d_i] = searchby_multfield({'type1','type2','subType'},'diogenite',spclib_relab_ftir_bdvnir_a);
[achondrite,a_i] = searchby_multfield({'type1','type2','subType'},'achondrite',spclib_relab_ftir_bdvnir_a);
[chondrite,c_i] = searchby_multfield({'type1','type2','subType'},'chondrite',spclib_relab_ftir_bdvnir_a);
[silicate,s_i] = searchby_multfield({'type1','type2','subType'},'silicate',spclib_relab_ftir_bdvnir_a);
[diopside,diop_i] = searchby_multfield({'type1','type2','subType'},'diopside',spclib_relab_ftir_bdvnir_a);

load Library_2018_ftir_bdbnir A Ahat spclib_relab_ftir_bdvnir_a

load library

%%
% convert radiance to I/F
sctime = 366393039;
sctime = 498259058;
vir_obs = get_vir_obs_info(sctime,'dwld',0);
I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);

SolSpcdata = VIRdata(vir_obs.basename_cal_ir_solspec,vir_obs.dirpath_cal_ir_1B);
SolSpcdata.readTAB();
SS = [SolSpcdata.tab.data.SPECTRAL_IRRADIANCE]';

I1Bim_if = vir_ra2if(I1Bdata,SolSpcdata);

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

%% second version considering emission
wv = I1Bdata.hdr.wavelength;
[~,b1420] = bisection(wv-1.42);
[~,b1560] = bisection(wv-1.56);
[~,b2410] = bisection(wv-2.41);
[~,b2610] = bisection(wv-2.61);
[~,b3740] = bisection(wv-3.74);
[~,b3830] = bisection(wv-3.83);
[~,b4350] = bisection(wv-4.35);
[~,b4540] = bisection(wv-4.54);
bL = length(wv);

bands = [b2610:b3740];

c = 49;
% bands = 180:294;
imgc = squeeze(I1Bim_if(:,c,:))';
imgc_b = imgc(bands,:);
wvb = wv(bands);

Alib_b = Ahat(bands,:);
valid_idx = ~any(isnan(Alib_b),1);
Alib_b_vld = Alib_b(:,valid_idx);
BBif_b = BBif(bands,:);
BBif_b_nrmed = normalizevec(BBif_b,1);

Alib_b_vld_nrmed = normalizevec(Alib_b_vld,1);
[ continua,bases ] = learnContinuum(wvb,Alib_b_vld_nrmed,'convhull');

[ Acntrmvd ] = CntRmvl(Alib_b_vld_nrmed,'additive','CONTINUA',...
                          continua );
Acntrmvd = -Acntrmvd;

% initialization parameters
load artifact_matrix.mat artifact_matrix
d0 = 1+artifact_matrix(bands,c);
imgc_mean = nanmean(imgc,2);
[imgc_mean_smthd] = vir_oddeven_rmvl_wInterp(imgc_mean,wv);
c00 = imgc_mean - imgc_mean_smthd;
c0 = c00(bands);

%%

load d00.mat d
d0 = d00(bands,c);

Alib_b_nrmed = normalizevec(Alib_b,1);
% [X,d1,c1] = vir_unmixing_denoiser_e(Alib_b_nrmed(:,[h_i e_i d_i a_i c_i]),BBif_b_nrmed,imgc_b,wvb,'C0',c0);

[X,d1,c1] = vir_unmixing_denoiser_e(Alib_b_nrmed(:,[h_i e_i d_i a_i c_i s_i]),BBif_b_nrmed,imgc_b,wvb,'C0',c0,'d0',d0,'lambda_a_lib',0.01);

% [X,d,c] = vir_unmixing_denoiser(Alib_b_vld_nrmed,imgc_b,wvb);




