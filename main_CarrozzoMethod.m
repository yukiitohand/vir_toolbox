

sctime = 366393039;
sctime = 498259058;
vir_obs = get_vir_obs_info(sctime,'dwld',0);
I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);

SolSpcdata = VIRdata(vir_obs.basename_cal_ir_solspec,vir_obs.dirpath_cal_ir_1B);
SolSpcdata.readTAB();
SS = [SolSpcdata.tab.data.SPECTRAL_IRRADIANCE]';

I1Bim_if = vir_ra2if(I1Bdata,SolSpcdata);


%%
load artifact_matrix_v2.mat artifact_matrix

artifact_matrix = reshape(artifact_matrix',[1,256,432]);

I1Bim_if_odd_even_rmvd = nan(size(I1Bim_if));
for c=1:256
    imgc = squeeze(I1Bim_if(:,c,:))';
%     imgc_b = imgc(bands,:);
    
    [imgc_oe] = vir_oddeven_rmvl_wInterp(imgc,wv);
    I1Bim_if_odd_even_rmvd(:,c,:) = reshape(imgc_oe',[size(I1Bim_if,1),1,432]); 
end
I1Bim_if_cor_carrozo = I1Bim_if_odd_even_rmvd./(1+artifact_matrix);


%%
b = 210;
r = [0.17 0.43];
r = [0.0265 0.0507];
% r = [0.0265 0.0507];
[ im ] = scx_rgb(I1Bim_if(:,:,b),r);
[ im_cor_carrozo ] = scx_rgb(I1Bim_if_cor_carrozo(:,:,b),r);
[ im_cor4 ] = scx_rgb(I1Bim_if_cor_v4(:,:,b),r);
% [ im_cor2 ] = scx_rgb(I1Bim_if_cor_v2(:,:,b),r);

figure; sc(im,'jet');
figure; sc(im_cor_carrozo,'jet');
figure; sc(im_cor4,'jet');
% figure; sc(im_cor2);

figure;
hold on;
c=50;
c=49;
plot(wvb,squeeze(I1Bim_if(5,50,bands)));

plot(wvb,squeeze(I1Bim_if_cor_carrozo(5,50,bands)));

plot(wvb,squeeze(I1Bim_if_cor(5,50,bands)));

figure; 
subplot(1,3,1);
plot(wvb,squeeze(I1Bim_if(:,c,bands))');
subplot(1,3,2);
plot(wvb,squeeze(I1Bim_if_cor_carrozo(:,c,bands))');
subplot(1,3,3);
plot(wvb,squeeze(I1Bim_if_cor(:,c,bands))');

%% for print
c=49;
figure;
set_figsize(gcf,367,420);
plot(wvb,squeeze(I1Bim_if_cor_carrozo(:,c,bands))');
xlim([2.58 3.8]);
set(gca,'FontSize',16);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F');
fname = sprintf('vir_%s_carrozo_c%03d.png',I1Bdata.basename,c);
export_fig(fname,'-transparent','-r150');


%% for print 2 ceres
figure;
set_figsize(gcf,361,243);
plot(wvb,squeeze(I1Bim_if_cor_carrozo(:,c,bands))','Color',[0.7 0.7 0.7]);
hold on;
plot(wvb,squeeze(I1Bim_if_cor_carrozo([15 44 108 118],c,bands))','LineWidth',1);
xlim([2.58 3.8]);
set(gca,'FontSize',16);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F');
fname = sprintf('vir_%s_carrozo_demo_c%03d.png',I1Bdata.basename,c);
export_fig(fname,'-transparent','-r150');

%% for print 2 ceres
figure;
set_figsize(gcf,361,243);
plot(wvb,squeeze(I1Bim_if_odd_even_rmvd(:,c,bands))','Color',[0.7 0.7 0.7]);
hold on;
plot(wvb,squeeze(I1Bim_if_odd_even_rmvd([15 44 108 118],c,bands))','LineWidth',1);
xlim([2.58 3.8]);
set(gca,'FontSize',16);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F');
fname = sprintf('vir_%s_odd_even_demo_c%03d.png',I1Bdata.basename,c);
export_fig(fname,'-transparent','-r150');

%% for print 3 ceres
figure;
set_figsize(gcf,361,243);
plot(wvb,squeeze(I1Bim_if(:,c,bands))','Color',[0.7 0.7 0.7]);
hold on;
plot(wvb,squeeze(I1Bim_if([15 44 108 118],c,bands))','LineWidth',1);
xlim([2.58 3.8]);
set(gca,'FontSize',16);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F');
fname = sprintf('vir_%s_c%03d_demo.png',I1Bdata.basename,c);
export_fig(fname,'-transparent','-r150');

%% for print 2 vesta
% c=49;
% c=100;
figure;
% set_figsize(gcf,300,300);
set_figsize(gcf,361,243);
plot(wvb,squeeze(I1Bim_if_cor_carrozo(:,c,bands))','Color',[0.7 0.7 0.7]);
hold on;
plot(wvb,squeeze(I1Bim_if_cor_carrozo([7 34 63 39],c,bands))','LineWidth',1);
% plot(wvb,squeeze(I1Bim_if_cor_carrozo([1 19 54 61],c,bands))','LineWidth',2);
xlim([2.58 3.8]);
% ylim([0.25 0.44]);
set(gca,'FontSize',16);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F');
fname = sprintf('vir_%s_carrozo_demo_c%03d.png',I1Bdata.basename,c);
export_fig(fname,'-transparent','-r150');

%% for print 2 vesta

figure;
set_figsize(gcf,361,243);
plot(wvb,squeeze(I1Bim_if_odd_even_rmvd(:,c,bands))','Color',[0.7 0.7 0.7]);
hold on;
plot(wvb,squeeze(I1Bim_if_odd_even_rmvd([7 34 63 39],c,bands))','LineWidth',1);
xlim([2.58 3.8]);
set(gca,'FontSize',16);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F');
fname = sprintf('vir_%s_odd_even_demo_c%03d.png',I1Bdata.basename,c);
export_fig(fname,'-transparent','-r150');

%% for print 3 vesta
figure;
set_figsize(gcf,300,300);
plot(wvb,squeeze(I1Bim_if(:,c,bands))','Color',[0.7 0.7 0.7]);
hold on;
% plot(wvb,squeeze(I1Bim_if([7 34 63 39],c,bands))','LineWidth',1);
plot(wvb,squeeze(I1Bim_if([1 19 54 61],c,bands))','LineWidth',1);
xlim([2.58 3.8]);
ylim([0.25 0.44]);
set(gca,'FontSize',16);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F');
fname = sprintf('vir_%s_c%03d_demo.png',I1Bdata.basename,c);
export_fig(fname,'-transparent','-r150');

%%
b = 210;
r = [0.0265 0.0507];
[ im ] = scx_rgb(I1Bim_if(:,:,b),r);
[ im_cor_carrozo ] = scx_rgb(I1Bim_if_cor_carrozo(:,:,b),r);
[ im_cor4 ] = scx_rgb(I1Bim_if_cor_v4(:,:,b),r);

fname = sprintf('vir_%s_if_b%03d.png',I1Bdata.basename,b);
[im_if] = sc(im,'jet');
imwrite(im_if,fname);

fname = sprintf('vir_%s_if_carrozzo_b%03d.png',I1Bdata.basename,b);
[im_carrozo] = sc(im_cor_carrozo,'jet');
imwrite(im_carrozo,fname);


fname = sprintf('vir_%s_if_proposed_v4_b%03d.png',I1Bdata.basename,b);
[im_4] = sc(im_cor4,'jet');
imwrite(im_4,fname);

figure; sc(im_cor_carrozo,'jet');
figure; sc(im_cor4,'jet');
