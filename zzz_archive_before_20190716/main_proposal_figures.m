figure;
hold on;
for i=0:4
    s = 1+0.1*i;
    I1Ad = I1Adata.img(:,:,:) - nanmean(I1Adata.img(1,:,:),1)*s;

    plot(I1Adata.hdr.wavelength,squeeze(I1Ad([50],100,:)),'DisplayName',sprintf('DN - Background x %3.1f',s));
end

set(gca,'FontSize',16);
legend();
set(gca,'YTickLabel',[]);
xlim([1 5]);
xlabel('Wavelength [{\mu}m]');
ylabel('DN - Background x const.');
fname = sprintf('vir_%s_background_subtract_demo_c%03d.png',I1Adata.basename,100);
export_fig(fname,'-transparent','-r150');



%%
figure; plot(sList,tvals,'LineWidth',2);
xlabel('Stretch Parameter')
ylabel('Total Variation');
set(gca,'FontSize',20);
export_fig('Total_variation.png','-transparent','-r150');

%%
sctime = 366390345; 
vir_obs = get_vir_obs_info(sctime,'dwld',0);
I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);
I1Adata = VIRdata(vir_obs.basename_ir_1A,vir_obs.dirpath_ir_1A);
I1Adata.readimg();
I1Bdata.readimg();
wv = I1Bdata.hdr.wavelength;

SolSpcdata = VIRdata(vir_obs.basename_cal_ir_solspec,vir_obs.dirpath_cal_ir_1B);
SolSpcdata.readTAB();
SS = [SolSpcdata.tab.data.SPECTRAL_IRRADIANCE]';

I1Bim_if = vir_ra2if(I1Bdata,SolSpcdata);
%%
c=217;
figure;
plot(wv,squeeze(I1Bim_if(:,c,:))','Color',[0.8 0.8 0.8]);
hold on;
plot(wv,squeeze(I1Bim_if([3 56 63],c,:))','LineWidth',1);
xlim([2.15 4.23]);
ylim([0.1 0.75]);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F');
set(gca,'FontSize',14);
export_fig('Figure_1.png','-transparent','-r150');
%%
c=217;
figure;
% plot(wv,squeeze(I1Bim_if(:,c,:))'-squeeze(I1Bim_if(61,c,:)),'Color',[0.8 0.8 0.8]);
hold on;
plot(wv,squeeze(I1Bim_if(:,c,:))'-squeeze(I1Bim_if(61,c,:)),'LineWidth',1);
xlim([2.15 4.23]);
ylim([-0.02 0.3]);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F - I/F_{darkest}');
set(gca,'FontSize',14);
export_fig('Figure_1_b.png','-transparent','-r150');

%%
figure;
hold on;
for i=0:10
    s = 1+0.05*i;
    I1Ad = I1Adata.img(:,:,:) - nanmean(I1Adata.img(1,:,:),1)*s;

    plot(I1Adata.hdr.wavelength,squeeze(I1Ad(37,c,:)),'DisplayName',sprintf('DN - Background x %3.2f',s));
end
legend
xlabel('Wavelength [{\mu}m]');
ylabel('DN - Background x const.');
set(gca,'YTick',[]);
set(gca,'FontSize',14);
export_fig('Figure_2.png','-transparent','-r150');
