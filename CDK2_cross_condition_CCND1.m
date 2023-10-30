clear all;
%% parameters to set
nuc_col     = 10;
cyto_col    = 15;
other_col   = 19;
framerate   = 5;
datadir     = 'Z:\20190912_2A7t_1AArepeat\Data\';
movie_leng  = 333;
sitemat     = [1:2];
plot_time   = [-8 50]; % time range to plot
plot_signal = [0 90]; % signal range of the protein of interest to plot 
drug_frame  = 94;
grab_time   = [4 5];

%% grabbing cell traces of interest
[cdk2_trace, other_trace] = grab_align('rowmat', [1:8], ...
    'colmat', [1], ...
    'sitemat', sitemat, ...
    'datapath', datadir, ...
    'cyto_col', cyto_col, ...
    'nuc_col', nuc_col, ...
    'other_col', other_col, ...
    'framerate', framerate, ...
    'movie_leng', movie_leng, ...
    'POI_name', 'mitosis', ...
    'grab_cdk2state','cdk2inc', ...
    'grab_generation','current',...
    'grab_time', grab_time,...
    'POI2_name', 'drug_addition',...
    'drug_frame', drug_frame);

[cdk2_trace2, other_trace2] = grab_align('rowmat', 1:8, ...
    'colmat', [10], ...
    'sitemat', sitemat, ...
    'datapath', datadir, ...
    'cyto_col', cyto_col, ...
    'nuc_col', nuc_col, ...
    'other_col', other_col, ...
    'framerate', framerate, ...
    'movie_leng', movie_leng, ...
    'POI_name', 'mitosis', ...
    'grab_cdk2state','cdk2inc', ...
    'grab_generation','current',...
    'grab_time', grab_time,...
    'POI2_name', 'drug_addition',...
    'drug_frame', drug_frame);

% [cdk2_trace3, other_trace3] = grab_align('rowmat', 7, ...
%     'colmat', [1:3 7:9], ...
%     'sitemat', sitemat, ...
%     'datapath', datadir, ...
%     'cyto_col', cyto_col, ...
%     'nuc_col', nuc_col, ...
%     'other_col', other_col, ...
%     'framerate', framerate, ...
%     'movie_leng', movie_leng, ...
%     'POI_name', 'mitosis', ...
%     'grab_cdk2state','cdk2inc', ...
%     'grab_generation','current',...
%     'grab_time', grab_time,...
%     'POI2_name', 'drug_addition',...
%     'drug_frame', drug_frame);


%% ploting
time   = (-movie_leng:movie_leng)/framerate;
% colors = [0 0 0;
%           0 0 1;
%           0.9100 0.4100 0.1700];

colors = [0 0 0;
         [5 147 72]/255;
         [232 22 140]/255];

% colors = [0 0 0;
%     [57 181 74]/255;
%     0 0 1];
% 0.9100 0.4100 0.1700 %Orange

figure, hold on

% subplot(2,2,1)
% patch([grab_time grab_time(2) grab_time(1)],[0 0 4 4],[0.85 0.85 0.85],'LineStyle','none');
% for i = 1:min(size(cdk2_trace,1),50)
%     line(time,cdk2_trace(i,:),'color','k');
%     hold on
% end
% for i = 1:min(size(cdk2_trace2,1),50)
%     line(time,cdk2_trace2(i,:),'color',colors(2,:));
%     hold on
% end
% % for i = 1:min(size(cdk2_trace3,1),50)
% %     line(time,cdk2_trace3(i,:),'color',colors(3,:));
% %     hold on
% % end
% ylim([0 4]);
% xlim(plot_time);
% xlabel('Time relative to anaphase (hr)');
% ylabel('CDK2 activity');
% % a = get(gca,'XTickLabel');set(gca,'fontsize',6);
% % a = get(gca,'YTickLabel');set(gca,'fontsize',6);
% 
% subplot(2,2,2)
% patch([grab_time grab_time(2) grab_time(1)],[plot_signal(1) plot_signal(1) plot_signal(2) plot_signal(2)],[0.85 0.85 0.85],'LineStyle','none');
% for i = 1:min(size(cdk2_trace,1),50)
%     line(time,other_trace(i,:),'color','k');
%     hold on
% end
% for i = 1:min(size(cdk2_trace2,1),50)
%     line(time,other_trace2(i,:),'color',colors(3,:));
%     hold on
% end
% xlim(plot_time);
% ylim(plot_signal);
% xlabel('Time relative to anaphase (hr)')
% ylabel('Cyclin D1 level')

subplot(2,2,1)
patch([grab_time grab_time(2) grab_time(1)],[0.2 0.2 2.5 2.5],[0.85 0.85 0.85],'LineStyle','none');
ciplot_pro(cdk2_trace,time,colors(1,:));
ciplot_pro(cdk2_trace2,time,colors(2,:));
% ciplot_pro(cdk2_trace3,time,colors(3,:));
xlim(plot_time);
ylim([0.2 2]);
xlabel('Time relative to anaphase (hr)');
ylabel('Mean CDK2 activity','color',colors(1,:));
% a = get(gca,'XTickLabel');set(gca,'fontsize',6);
% a = get(gca,'YTickLabel');set(gca,'fontsize',6);

subplot(2,2,3)
patch([grab_time grab_time(2) grab_time(1)],[plot_signal(1) plot_signal(1) plot_signal(2) plot_signal(2)],[0.85 0.85 0.85],'LineStyle','none');
ciplot_pro(other_trace,time,colors(1,:));
ciplot_pro(other_trace2,time,colors(3,:));
xlim(plot_time);
ylim([0 50]);
xlabel('Time relative to anaphase (hr)')
ylabel('Mean Cyclin D1 level','color',colors(1,:))
% 
% set(gcf, 'PaperPosition', [0 0 1.5 2.25]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [1.5 2.25]);
% 
% print(gcf, '-dpdf', '-r300', 'test.pdf');