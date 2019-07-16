figure;
set_figsize(gcf,361,243);
plot(wvb,squeeze(I1Bim_if(:,49,bands))','Color',[0.7 0.7 0.7]);
hold on;
plot(wvb,squeeze(I1Bim_if([15 44 108 118],49,bands))','LineWidth',1);
xlim([2.58 3.8]);
ylim([-0.0164 0.0836]);
set(gca,'FontSize',12);
xlabel('Wavelength [{\mu}m]');
ylabel('I/F');
fname = sprintf('vir_%s_c%03d_demo_v2.png',I1Bdata.basename,49);
export_fig(fname,'-transparent','-r150');

%%
figure;
set_figsize(gcf,361,243);
% set_figsize(gcf,367,420);
plot(wvb,squeeze(I1Bim_if_odd_even_rmvd(:,49,bands))','Color',[0.7 0.7 0.7]);
hold on;
plot(wvb,squeeze(I1Bim_if_odd_even_rmvd([15 44 108 118],49,bands))','LineWidth',1);
yyaxis left
ylim([-0.0164 0.0836]);
yyaxis right
plot(wvb,(1+squeeze(artifact_matrix(1,49,bands))),'LineWidth',1);
ylim([0.8349 1.5349]);
xlim([2.56,3.8]);
set(gca,'FontSize',12);

set(gca,'XTick',[2.6 wvb(48) wvb(76) wvb(93) 3.6]);
set(gca,'XGrid','on');
yyaxis left;
ylabel('I/F');
xlabel('Wavelength[{\mu}m]');

fname = sprintf('vir_%s_odd_even_intro_c%03d_v2.png',I1Bdata.basename,49);
export_fig(fname,'-transparent','-r150');

%%
figure;
set_figsize(gcf,361,243);
% set_figsize(gcf,367,420);
% set_figsize(gcf,361,243);
plot(wvb,squeeze(I1Bim_if_cor_carrozo(:,49,bands))','Color',[0.7 0.7 0.7]);
hold on;
plot(wvb,squeeze(I1Bim_if_cor_carrozo([15 44 108 118],49,bands))','LineWidth',1);
yyaxis left
ylim([-0.0164 0.0836]);
yyaxis right
plot(wvb,(1+squeeze(artifact_matrix(1,49,bands))),'LineWidth',1);
ylim([0.8349 1.5349]);
xlim([2.56,3.8]);
set(gca,'FontSize',12);

set(gca,'XTick',[2.6 wvb(48) wvb(76) wvb(93) 3.6]);
set(gca,'XGrid','on');
yyaxis left;
ylabel('I/F');
xlabel('Wavelength[{\mu}m]');

fname = sprintf('vir_%s_carrozo_intro_c%03d_v2.png',I1Bdata.basename,49);
export_fig(fname,'-transparent','-r150');

%%
figure;
% set_figsize(gcf,361,243);
set_figsize(gcf,367,420);
plot(wvb,squeeze(I1Bim_if_odd_even_rmvd(:,49,bands))');
yyaxis left
ylim([-0.0164 0.0836]);
yyaxis right
plot(wvb,(1+squeeze(artifact_matrix(1,49,bands))),'LineWidth',1);
ylim([0.8349 1.5349]);
xlim([2.56,3.8]);
set(gca,'FontSize',12);

set(gca,'XTick',[2.6 wvb(48) wvb(76) wvb(93) 3.6]);
set(gca,'XGrid','on');
yyaxis left;
ylabel('I/F');
xlabel('Wavelength[{\mu}m]');

fname = sprintf('vir_%s_odd_even_demo_c%03d_v2.png',I1Bdata.basename,49);
export_fig(fname,'-transparent','-r150');

%%
figure;
% set_figsize(gcf,361,243);
set_figsize(gcf,367,420);
% set_figsize(gcf,361,243);
plot(wvb,squeeze(I1Bim_if_cor_carrozo(:,49,bands))');
yyaxis left
ylim([-0.0164 0.0836]);
yyaxis right
plot(wvb,(1+squeeze(artifact_matrix(1,49,bands))),'LineWidth',1);
ylim([0.8349 1.5349]);
xlim([2.56,3.8]);
set(gca,'FontSize',12);

set(gca,'XTick',[2.6 wvb(48) wvb(76) wvb(93) 3.6]);
set(gca,'XGrid','on');
yyaxis left;
ylabel('I/F');
xlabel('Wavelength[{\mu}m]');

fname = sprintf('vir_%s_carrozo_demo_c%03d_v2.png',I1Bdata.basename,49);
export_fig(fname,'-transparent','-r150');

%%

%%
figure;
set_figsize(gcf,361,243);
set_figsize(gcf,367,420);
plot(wvb,squeeze(I1Bim_if_odd_even_rmvd(:,50,bands))');
yyaxis left
ylim([-0.0164 0.0836]);
yyaxis right
plot(wvb,(1+squeeze(artifact_matrix(1,50,bands))),'LineWidth',1);
ylim([0.8349 1.5349]);
xlim([2.56,3.8]);
set(gca,'FontSize',12);

set(gca,'XTick',[2.6 wvb(26) wvb(66) wvb(76) wvb(93) 3.6]);
set(gca,'XGrid','on');
yyaxis left;
ylabel('I/F');
xlabel('Wavelength[{\mu}m]');

fname = sprintf('vir_%s_odd_even_demo_c%03d_v2.png',I1Bdata.basename,50);
export_fig(fname,'-transparent','-r150');

%%
figure;
set_figsize(gcf,367,420);
% set_figsize(gcf,361,243);
plot(wvb,squeeze(I1Bim_if_cor_carrozo(:,50,bands))');
yyaxis left
ylim([-0.0164 0.0836]);
yyaxis right
plot(wvb,(1+squeeze(artifact_matrix(1,50,bands))),'LineWidth',1);
ylim([0.8349 1.5349]);
xlim([2.56,3.8]);
set(gca,'FontSize',12);

set(gca,'XTick',[2.6 wvb(26) wvb(66) wvb(83) wvb(93) 3.6]);
set(gca,'XGrid','on');
yyaxis left;
ylabel('I/F');
xlabel('Wavelength[{\mu}m]');

fname = sprintf('vir_%s_carrozo_demo_c%03d_v2.png',I1Bdata.basename,50);
export_fig(fname,'-transparent','-r150');
