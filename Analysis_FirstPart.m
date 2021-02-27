%%

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',12)
set(0,'defaultLegendFontSize',20)
cb = brewermap(30,'*RdYlBu');
cbFan_PCx = brewermap([],'*Blues');
cbFan_PRh = brewermap([],'*Reds');
cb_PCx = cb(6,:);
cb_PRh = cb(24,:);

cbz = brewermap(18,'*Spectral');
cb_setA = cbz(5,:);
cb_setB = cbz(13,:);

%% Analisi Tesi


%% Baseline Activity

load('database_PCx_THESIS');
load('cfg_PCx_THESIS');
tempVect= []; tempValues= []; PCx_bsl = [];
totSecs = 8 * 15 *10; i = 0;

for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        i = i + 1;
        tempVect = [];
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)
            for idxTrial = 1:length(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikesPerTrial)
                tempValues = [exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikesPerTrial(idxTrial).Trial];
                listSpikes = tempValues(tempValues<0);
                tempVect = [tempVect listSpikes'];
            end
        end
        PCx_bsl = [PCx_bsl numel(tempVect)/totSecs];
    end
end

load('database_PRh_THESIS');
load('cfg_PRh_THESIS');
tempVect= []; tempValues= []; PRh_bsl = [];
i = 0;

for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        i = i + 1;
        tempVect = [];
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)
            for idxTrial = 1:length(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikesPerTrial)
                tempValues = [exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikesPerTrial(idxTrial).Trial];
                listSpikes = tempValues(tempValues<0);
                tempVect = [tempVect listSpikes'];
            end
        end
        PRh_bsl = [PRh_bsl numel(tempVect)/totSecs];
    end
end

%[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
 figure('Position', [10 10 900 400])
d{1} =  PCx_bsl; d{2} = PRh_bsl;
h1 = raincloud_plot(d{1}, 'band_width',1,'box_on', 1, 'color', cb_PCx, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .10, 'dot_dodge_amount', .29,...
    'box_col_match', 1);
h2 = raincloud_plot(d{2},'band_width',1, 'box_on', 1, 'color', cb_PRh, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .48, 'dot_dodge_amount', .67, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';
h1{2}.SizeData = 4; h2{2}.SizeData = 4;
xlabel('Spikes/s');
ylabels = ylabel('PDF');
ylim([-0.195 0.30])
xlim([0 30]); yticks([0:0.1:0.3])
legend([h1{1},h2{1}],'PCx','PRh')
legend boxoff
set(gcf,'color','white'); set(gca, 'box', 'off', 'tickDir', 'out');
title('Spontaneous Activity during Baseline') 
[p,h]=ranksum(PCx_bsl,PRh_bsl)
ylabels.Position(2) = 0.15;
median(PCx_bsl)
median(PRh_bsl)
%% Clustering 
load('database_PCx_THESIS');
load('cfg_PCx_THESIS');
trough2peak = [];
asymmetry = [];
fr = [];
X = []; Xz = [];
for idxExp = 1:length(cfg.folders.folderList)
    for idxUnit = 1:size(exp.expID(idxExp).unit,2)
        trough2peak = [trough2peak exp.expID(idxExp).unit(idxUnit).trough2peak_latency];
        asymmetry = [asymmetry exp.expID(idxExp).unit(idxUnit).peak_asymmetry];
        fr = [fr exp.expID(idxExp).unit(idxUnit).meanFiringRate ];
    end
end

X(:,1) = trough2peak';
X(:,2) = asymmetry';
X(:,3) = fr';
Xz = zscore(X);
% Xz = zscore(X(:,1:2));

set(0,'defaultAxesFontSize',14)
set(0,'defaultLegendFontSize',22)
kidx = kmeans(Xz,2, 'Distance', 'sqeuclidean', 'MaxIter', 10000, 'Replicates', 1000);
figure('Renderer', 'painters', 'Position', [10 10 1400 900])
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
bro = brewermap(15,'*RdYlGn');
tiledlayout(2,2,'TileSpacing','compact')
nexttile
l = scatter3(X(kidx==1,1), X(kidx==1,2), X(kidx==1,3), 18,bro(3,:),'filled');
hold on
f =scatter3(X(kidx==2,1), X(kidx==2,2), X(kidx==2,3), 18, bro(12,:),'filled');
xl = xlabel('T2P'), yl = ylabel('Asym'), zlabel('FR');
set(gcf,'color','white', 'PaperPositionMode', 'auto');
title('PCx')
axis square
set(gca, 'CameraPosition', [-262.783861231323,-12.5607196171800,102.930423464843]);

xlim([0 50]); ylim([-1 1]); zlim([0 40]);
set(get(gca,'ylabel'),'rotation',341.9);
set(get(gca,'ylabel'),'Position',[4.3 .25 -4.5]);
set(get(gca,'xlabel'),'rotation',13);
set(get(gca,'xlabel'),'Position',[35 -0.6 -4.5]);
set(get(gca,'zlabel'),'rotation',360);
set(get(gca,'zlabel'),'Position',[-36 0 27 ]);

PCx_RS= length(find(kidx ==2));
PCx_FS =length(find(kidx ==1));

load('database_PRh_THESIS');
load('cfg_PRh_THESIS');
trough2peak = [];
asymmetry = [];
fr = [];
X =[];
Xz = [];
for idxExp = 1:length(cfg.folders.folderList)
    for idxUnit = 1:size(exp.expID(idxExp).unit,2)
        trough2peak = [trough2peak exp.expID(idxExp).unit(idxUnit).trough2peak_latency];
        asymmetry = [asymmetry exp.expID(idxExp).unit(idxUnit).peak_asymmetry];
        fr = [fr exp.expID(idxExp).unit(idxUnit).meanFiringRate ];
    end
end

X(:,1) = trough2peak';
X(:,2) = asymmetry';
X(:,3) = fr';
Xz = zscore(X);
% Xz = zscore(X(:,1:2));

kidx = kmeans(Xz,2, 'Distance', 'sqeuclidean', 'MaxIter', 10000, 'Replicates', 1000);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
bro = brewermap(15,'*RdYlGn');
nexttile(3)
l = scatter3(X(kidx==1,1), X(kidx==1,2), X(kidx==1,3), 18,bro(3,:),'filled');
hold on
f =scatter3(X(kidx==2,1), X(kidx==2,2), X(kidx==2,3), 18, bro(12,:),'filled');
xl = xlabel('T2P'), yl = ylabel('Asym'), zlabel('FR');
title('PRh')
axis square
set(gca, 'CameraPosition', [-262.783861231323,-12.5607196171800,102.930423464843]);
% lgd = legend('Regular Spiking', 'Fast Spiking','FontSize',12','Interpreter','latex');
% legend boxoff
%set(gca, 'XColor', 'w', 'box','off');
set(gcf,'color','white', 'PaperPositionMode', 'auto');
xlim([0 50]); ylim([-1 1]); zlim([0 40]);
set(get(gca,'ylabel'),'rotation',341.9);
set(get(gca,'ylabel'),'Position',[4.3 .25 -4.5]);
set(get(gca,'xlabel'),'rotation',13);
set(get(gca,'xlabel'),'Position',[35 -0.6 -4.5]);
set(get(gca,'zlabel'),'rotation',360);
set(get(gca,'zlabel'),'Position',[-36 0 27 ]);

nexttile(2,[2,1])
PRh_RS = length(find(kidx == 2));
PRh_FS = length(find(kidx == 1));
sbar = bar([PCx_FS/(PCx_FS+PCx_RS),PCx_RS/(PCx_FS+PCx_RS);PRh_FS/(PRh_FS+PRh_RS),PRh_RS/(PRh_FS+PRh_RS)], 'grouped');
sbar(1, 1).FaceColor =  bro(3,:); sbar(1, 2).FaceColor =  bro(12,:);
sbar(1, 1).LineStyle =  'none'; sbar(1, 2).LineStyle =  'none';
title('Putative Cell Type Fraction')
ylabel('Fraction of Neurons')
xticklabels({'PCx','PRh'})
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box','off');
ylim([0 1]);
nbars = size([PCx_FS/(PCx_FS+PCx_RS),PCx_RS/(PCx_FS+PCx_RS);PRh_FS/(PRh_FS+PRh_RS),PRh_RS/(PRh_FS+PRh_RS)], 2);
error_PCx_FS = sqrt((PCx_FS/(PCx_FS+PCx_RS))*(1-(PCx_FS/(PCx_FS+PCx_RS)))/(PCx_FS+PCx_RS))
error_PCx_RS = sqrt((PCx_RS/(PCx_FS+PCx_RS))*(1-(PCx_RS/(PCx_FS+PCx_RS)))/(PCx_FS+PCx_RS))
error_PRh_FS = sqrt((PRh_FS/(PRh_FS+PRh_RS))*(1-(PRh_FS/(PRh_FS+PRh_RS)))/(PRh_FS+PRh_RS))
error_PRh_RS = sqrt((PRh_RS/(PRh_FS+PRh_RS))*(1-(PRh_RS/(PRh_FS+PRh_RS)))/(PRh_FS+PRh_RS))
%Get the x coordinate of the bars
hold on
x = [];
for i = 1:nbars
    x = [x ; sbar(i).XEndPoints];
end
%Plot the errorbars
errb = errorbar(x',[PCx_FS/(PCx_FS+PCx_RS),PCx_RS/(PCx_FS+PCx_RS);PRh_FS/(PRh_FS+PRh_RS),PRh_RS/(PRh_FS+PRh_RS)],[error_PCx_FS error_PCx_RS; error_PRh_FS error_PRh_RS],'k','linestyle','none')'
errb(2,1).Color = bro(12,:)
errb(1,1).Color = bro(3,:)
% hold off
lg  = legend('Regular Spiking', 'Fast Spiking', 'Location', 'northwest'); 
legend boxoff
% lg.Layout.Tile = 'east';

%%
 % Observed data
 n1 = PCx_FS; N1 = PCx_FS+PCx_RS;
 n2 = PRh_FS; N2 = PRh_FS+PRh_RS;
 x1 = [repmat('a',N1,1); repmat('b',N2,1)];
 x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
 [tbl,chi2stat,pval] = crosstab(x1,x2)

n1 = PCx_FS; N1 = PCx_FS+PCx_RS;
n2 = PRh_FS; N2 = PRh_FS+PRh_RS;
% Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2)
% Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
chi2stat = sum((observed-expected).^2 ./ expected)


n1 = PCx_FS; N1 = PCx_FS+PCx_RS;
n2 = PRh_FS; N2 = PRh_FS+PRh_RS;
% Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2)
% Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
[h,p,stats] = chi2gof([1 2 3 4],'freq',observed,'expected',expected,'ctrs',[1 2 3 4],'nparams',2)

%% Spike Response Distribution
figure
load('database_PCx_THESIS');
load('cfg_PCx_THESIS');
PCx_Mat = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        tempMat = mean([exp.expID(idxExp).unit(idxUnit).event.spikeResponse]);
        PCx_Mat = [PCx_Mat ; tempMat];
    end
end

load('database_PRh_THESIS');
load('cfg_PRh_THESIS');
PRh_Mat = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        tempMat = mean([exp.expID(idxExp).unit(idxUnit).event.spikeResponse]);
        PRh_Mat = [PRh_Mat ; tempMat];
    end
end

A_PCx = nonzeros(PCx_Mat(:));
A_PRh = nonzeros(PRh_Mat(:));
pcx = histcounts(A_PCx, [-20:1:20]);
prh = histcounts(A_PRh, [-20:1:20]);
hPCx = plot(1:40, pcx/length(PCx_Mat(:))); hold on
hPRh = plot(1:40, prh/length(PRh_Mat(:)));
hPCx.LineWidth = 2; hPRh.LineWidth = 2;
hPCx.Color = cb_PCx; hPRh.Color = cb_PRh;
xlim([1 41])
xticks(1:5:41), xticklabels({[-20:5:20]})
legend('PCx','PRh')
legend boxoff

%% Raster Plots

 figure('Position', [10 10 900 600])
 
%% Fraction of Neurons with E/I responses

figure
load('database_PCx_THESIS');
load('cfg_PCx_THESIS');
PCx_Exc = 0; PCx_Inh = 0; PCx_counter = 0;
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        PCx_counter = PCx_counter + 1;
        Exc_C = 0; Inh_C = 0;
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse == 1
                Inh_C = Inh_C + 1;
            end
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse== 1
                Exc_C = Exc_C + 1;
            end
        end
        if Inh_C > 0
            PCx_Inh = PCx_Inh +1;
        elseif Exc_C>0
            PCx_Exc = PCx_Exc +1;
        end
    end
end


load('database_PRh_THESIS');
load('cfg_PRh_THESIS');
PRh_Exc = 0; PRh_Inh = 0; PRh_counter = 0;  a = 0;
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        PRh_counter = PRh_counter + 1;
        Exc_C = 0; Inh_C = 0; l = 0;
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse == 1
                Inh_C = Inh_C + 1;
            end
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse== 1
                Exc_C = Exc_C + 1;
            end
        end
        if Inh_C > 0
            PRh_Inh = PRh_Inh +1;
        elseif Exc_C>0
            PRh_Exc = PRh_Exc +1;
        end
    end
end
set(0,'defaultAxesFontSize',12)
set(0,'defaultLegendFontSize',20)
figure('Renderer', 'painters', 'Position', [10 10 900 300])
b = bar([PCx_Exc/PCx_counter,PRh_Exc/PRh_counter;PCx_Inh/PCx_counter,PRh_Inh/PRh_counter])
b(1,1).FaceColor = cb_PCx, b(1,1).LineStyle = 'none'
b(1,2).FaceColor = cb_PRh, b(1,2).LineStyle = 'none'
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
ylabel('Fraction of Neurons'), xticklabels({'Excitatory','Inhibitory'});
 text(0.96, 0.6, '***','FontSize',14);
 nbins = 2
 hold on
x = [];
for i = 1:nbins
    x = [x ; b(i).XEndPoints];
end
%Plot the errorbars
error_PCx_Exc = sqrt((PCx_Exc/(PCx_counter))*(1-(PCx_Exc/(PCx_counter)))/(PCx_counter))
error_PCx_Inh = sqrt((PCx_Inh/(PCx_counter))*(1-(PCx_Inh/(PCx_counter)))/(PCx_counter))
error_PRh_Exc = sqrt((PRh_Exc/(PRh_counter))*(1-(PRh_Exc/(PRh_counter)))/(PRh_counter))
error_PRh_Inh = sqrt((PRh_Inh/(PRh_counter))*(1-(PRh_Inh/(PRh_counter)))/(PRh_counter))
errb = errorbar(x',[PCx_Exc/PCx_counter,PRh_Exc/PRh_counter;PCx_Inh/PCx_counter,PRh_Inh/PRh_counter],[error_PCx_Exc error_PRh_Exc; error_PCx_Inh error_PRh_Inh],'k','linestyle','none')'
errb(1,1).Color = cb_PCx
errb(2,1).Color = cb_PRh
legend('PCx','PRh')
ylim([0 1]);
legend boxoff

 n1 = 712; N1 = 2996;
 n2 = 1741; N2 = 7350;
 x1 = [repmat('a',N1,1); repmat('b',N2,1)];
 x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
 [tbl,chi2stat,pval] = crosstab(x1,x2)

 n1 = PCx_Inh; N1 = PCx_counter;
 n2 = PRh_Inh; N2 = PRh_counter;
 x1 = [repmat('a',N1,1); repmat('b',N2,1)];
 x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
 [tbl,chi2stat,pval] = crosstab(x1,x2)

%% auROCs distributions



load('cfg_PCx_THESIS');
load('database_PCx_THESIS');
auROCs_PCx = [];
inH = [];
exC = [];

for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        %if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                auROCs_PCx = [auROCs_PCx exp.expID(idxExp).unit(idxUnit).event(idxEvent).auROC];
                inH = [ inH exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse];
                exC = [ exC exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse];
            end
        %end
    end
end
figure('Renderer', 'painters', 'Position', [10 10 900 300])
%subplot(1,2,1)
tiledlayout(1,2)
nexttile
pl = 0:0.05:1;
h1 = histogram(auROCs_PCx, pl, 'Normalization', 'probability');
h1.FaceColor = 'none';
h1.EdgeColor = cb_PCx;
h1.EdgeAlpha = 1;
h1.LineWidth = 1;

excNaN = length(auROCs_PCx) - length(auROCs_PCx(exC==1));
hold on
exc1 = histogram([auROCs_PCx(exC==1) nan(1,excNaN)], pl, 'Normalization', 'probability');
exc1.FaceColor = cb_PCx;
exc1.FaceAlpha = .85;
exc1.EdgeColor = [0 0 0];
exc1.EdgeAlpha = 0;

inhNaN = length(auROCs_PCx) - length(auROCs_PCx(inH==1));
hold on
inh = histogram([auROCs_PCx(inH==1) nan(1,inhNaN)], pl, 'Normalization', 'probability');
inh.FaceColor = cb_PCx;
inh.FaceAlpha = .85;
inh.EdgeColor = [0 0 0];
inh.EdgeAlpha = 0;
title('PCx');
xlabel('auROC'); ylabel('Fraction of Cell-Odor Pairs');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
% legend(exc1, 'Significant Responses')
% legend boxoff
ylim([0 0.2]);

load('database_PRh_THESIS');
load('cfg_PRh_THESIS');
auROCs_PRh = [];
inH = [];
exC = [];

for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        %if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                auROCs_PRh = [auROCs_PRh exp.expID(idxExp).unit(idxUnit).event(idxEvent).auROC];
                inH = [ inH exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse];
                exC = [ exC exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse];
            end
        %end
    end
end
%subplot(1,2,2)
nexttile
pl = 0:0.05:1;
h2 = histogram(auROCs_PRh, pl, 'Normalization', 'probability');
h2.FaceColor = 'none';
h2.EdgeColor = cb_PRh;
h2.EdgeAlpha = 1;
h2.LineWidth = 1;

excNaN = length(auROCs_PRh) - length(auROCs_PRh(exC==1));
hold on
exc2 = histogram([auROCs_PRh(exC==1) nan(1,excNaN)], pl, 'Normalization', 'probability');
exc2.FaceColor = cb_PRh;
exc2.FaceAlpha = .85;
exc2.EdgeColor = [0 0 0];
exc2.EdgeAlpha = 0;

inhNaN = length(auROCs_PRh) - length(auROCs_PRh(inH==1));
hold on
inh = histogram([auROCs_PRh(inH==1) nan(1,inhNaN)], pl, 'Normalization', 'probability');
inh.FaceColor = cb_PRh;
inh.FaceAlpha = .85;
inh.EdgeColor = [0 0 0];
inh.EdgeAlpha = 0;
title('PRh');
xlabel('auROC');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');set(gca, 'YTick', []);
ax1 = gca;            
ax1.YAxis.Visible = 'off';
% legend(exc2, 'Significant Responses')
% legn = legend([h1, h2],'PCx','PRh')
% legn.Position = [0.8 0.8 0.1 0.1]
ylim([0 0.2]);

[p,h]=ranksum(auROCs_PRh,auROCs_PCx)
[h,p]=kstest2(auROCs_PRh,auROCs_PCx)
[h,p]=ttest2(auROCs_PRh,auROCs_PCx)
mean(auROCs_PRh),mean(auROCs_PCx)
median(auROCs_PRh),median(auROCs_PCx)
sK_PRh = skewness(auROCs_PRh)
sK_PCx = skewness(auROCs_PCx)
%% NResponses

figure('Renderer', 'painters', 'Position', [10 10 900 300])
tiledlayout(1,2,'TileSpacing','compact')
nexttile
load('cfg_PCx_THESIS.mat');
load('database_PCx_THESIS.mat');
listNeuron = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        excCounter = 0;
        inhCounter = 0;
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                inhCounter = inhCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse;
                excCounter = excCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse;
        end
        listNeuron = [listNeuron ; excCounter inhCounter];
    end
end

listNeuron_PCx = listNeuron
xRange = 0:1:14;
[N_PCx, edges] = histcounts(listNeuron(:,1), -0.5:1:14.5);
N_PCx = N_PCx ./ numel(listNeuron(:,1));
t0 = plot(xRange, N_PCx, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', cb_PCx, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 7)
hold on
errExc_PCx = sqrt((N_PCx.*(1-N_PCx))/numel(listNeuron(:,1)));
errorbar(xRange,  N_PCx, errExc_PCx, 'color', cb_PCx,'LineStyle', 'none')
xlim([-1 15]);
xticks([0:1:14]);
title('Excitatory');
xlabel('Number of Odors');
ylabel('Fraction of Neurons');
ylim([-0.05 1]);
hold on

load('database_PRh_THESIS.mat');
load('cfg_PRh_THESIS.mat');
listNeuron = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        excCounter = 0;
        inhCounter = 0;
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                inhCounter = inhCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse;
                excCounter = excCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse;
        end
        listNeuron = [listNeuron ; excCounter inhCounter];
    end
end
listNeuron_PRh = listNeuron
[N_PRh, edges] = histcounts(listNeuron(:,1), -0.5:1:14.5);
N_PRh = N_PRh ./ numel(listNeuron(:,1));
t0 = plot(xRange, N_PRh, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', cb_PRh, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 7)
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
hold on
errExc_PRh = sqrt((N_PRh.*(1-N_PRh))/numel(listNeuron(:,1)));
errorbar(xRange,  N_PRh, errExc_PRh , 'color', cb_PRh,'LineStyle', 'none');



load('cfg_PCx_THESIS.mat');
load('database_PCx_THESIS.mat');
listNeuron = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        excCounter = 0;
        inhCounter = 0;
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                inhCounter = inhCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse;
                excCounter = excCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse;
        end
        listNeuron = [listNeuron ; excCounter inhCounter];
    end
end
nexttile
xRange = 0:1:14;
[N_PCx, edges] = histcounts(listNeuron(:,2), -0.5:1:14.5);
N_PCx = N_PCx ./ numel(listNeuron(:,2));
tz = plot(xRange, N_PCx, 'LineStyle', 'none', 'Marker', 'O', 'MarkerFaceColor', cb_PCx, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 7)
xlim([-1 15]);
xticks([0:1:14]);
title('Inhibitory');
xlabel('Number of Odors');
ylabel('Fraction of Neurons');

hold on
errInh_PCx = sqrt((N_PCx.*(1-N_PCx))/numel(listNeuron(:,1)));
errorbar(xRange,  N_PCx, errInh_PCx, 'color', cb_PCx,'LineStyle', 'none')
ylim([-0.05 1]);

load('database_PRh_THESIS.mat');
load('cfg_PRh_THESIS.mat');
listNeuron = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        excCounter = 0;
        inhCounter = 0;
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                inhCounter = inhCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse;
                excCounter = excCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse;
        end
        listNeuron = [listNeuron ; excCounter inhCounter];
    end
end

[N_PRh, edges] = histcounts(listNeuron(:,2), -0.5:1:14.5);
N_PRh = N_PRh ./ numel(listNeuron(:,2));
t1 = plot(xRange, N_PRh, 'LineStyle', 'none', 'Marker', 'O', 'MarkerFaceColor', cb_PRh, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 7)
hold on
errInh_PRh = sqrt((N_PRh.*(1-N_PRh))/numel(listNeuron(:,1)));
errorbar(xRange,  N_PRh, errInh_PRh , 'color', cb_PRh,'LineStyle', 'none');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');set(gca, 'YTick', []);
ax1 = gca;            
ax1.YAxis.Visible = 'off';
legend([tz, t1],{'PCx', 'PRh'});
legend boxoff

[h,p]=kstest2(listNeuron_PCx(:,1),listNeuron_PRh(:,1))
[h,p]=kstest2(listNeuron_PCx(:,2),listNeuron_PRh(:,2))

 n1 = N_PCx(1); N1 = sum(N_PCx);
 n2 = N_PRh(1); N2 = sum(N_PRh);
 x1 = [repmat('a',N1,1); repmat('b',N2,1)];
 x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
 [tbl,chi2stat,pval] = crosstab(x1,x2)


%% Reliability

load('database_PCx_THESIS');
load('cfg_PCx_THESIS');
sR = []; indices = []; exC = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
            sR = exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse;
            indices = [indices length(find(sR > 0))];
            exC = [exC exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse];
        end
    end
end

subplot(1,2,1)
extnan = nan(1,length(indices) - length(indices(exC==1)));
pl = -.5:1:10.5;
h = histogram(indices, pl, 'Normalization', 'probability');
h.FaceColor = 'none';
h.EdgeColor = cb_PCx;
h.EdgeAlpha = 1;
hold on
h = histogram([indices(exC ==1) extnan], pl, 'Normalization', 'probability');
h.FaceColor = cb_PCx;
h.EdgeColor = 'none';
h.FaceAlpha = .85;
ylim([0, 0.3]);
title('PCx'); ylabel('Fraction of Excitatory Responses'); xlabel('Response Probability');
xticks(0:1:10); xtickangle(45);
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');

load('cfg_PRh_THESIS');
load('database_PRh_THESIS');
sR = []; indices = []; exC = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
            sR = exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse;
            indices = [indices length(find(sR > 0))];
            exC = [exC exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse];
        end
    end
end
subplot(1,2,2)
extnan = nan(1,length(indices) - length(indices(exC==1)));
pl = -.5:1:10.5;
h = histogram(indices, pl, 'Normalization', 'probability');
h.FaceColor = 'none';
h.EdgeColor = cb_PRh;
h.EdgeAlpha = 1;
hold on
h = histogram([indices(exC ==1) extnan], pl, 'Normalization', 'probability');
h.FaceColor = cb_PRh;
h.EdgeColor = 'none';
h.FaceAlpha = .85;
ylim([0, 0.3]);
title('PRh'); ylabel('Fraction of Excitatory Responses'); xlabel('Response Probability');
xticks(0:1:10); xtickangle(45);
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');

% sistemato con significance trial

%% Fano Factor
load('database_PCx_THESIS');
load('cfg_PCx_THESIS');
l = 0;
vect_PCx = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
             if  exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                l = l + 1;
                vect_PCx(l,:) = exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse;
             end
        end
    end
end

fanoF_PCx = var(vect_PCx,[],2)./mean(vect_PCx,2);
fastTr_PCx = mean(vect_PCx(:,1:5),2);
lastTr_PCx = mean(vect_PCx(:,6:10),2);

load('database_PRh_THESIS');
load('cfg_PRh_THESIS');
l = 0; vect_PRh = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
            if  exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                l = l + 1;
                vect_PRh(l,:) = exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse;
            end
        end
    end
end

fanoF_PRh = var(vect_PRh,[],2)./mean(vect_PRh,2);
fastTr_PRh = mean(vect_PRh(:,1:5),2);
lastTr_PRh = mean(vect_PRh(:,6:10),2);

d{1} =  fanoF_PCx; d{2} = fanoF_PRh
h1 = raincloud_plot(d{1}, 'band_width',1,'box_on', 1, 'color', cb_PCx, 'alpha', 0.9,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .25,...
    'box_col_match', 1);
h2 = raincloud_plot(d{2},'band_width',1, 'box_on', 1, 'color', cb_PRh, 'alpha', 0.9,...
    'box_dodge', 1, 'box_dodge_amount', .45, 'dot_dodge_amount', .55, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';
ylim([-0.1 0.2]); xlim([0 70]);title('Fano Factor Distributions for Excitatory Neurons');
ylabel('Fraction of Neurons'); xlabel('Fano Factor');
legend([h1{1},h2{1}],'PCx','PRh'); yticks([0:0.05:0.2]);
legend boxoff 
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
yline(1)
% la girerei verticalmente
camroll(90)
[h,p,l]=ttest2(fanoF_PCx,fanoF_PRh)
[p,h] = ranksum(fanoF_PCx,fanoF_PRh)
%%

figure( 'Position', [10 10 1400 800])
for ii = 1:length(fastTr_PCx)
    colsel = randi(50)+90;
    HabitPCx = plot(.5:1.5, [fastTr_PCx(ii) lastTr_PCx(ii)], '-o', 'color', cbFan_PCx(colsel,:),...
        'MarkerFaceColor', cbFan_PCx(colsel,:),'MarkerSize',2,'LineWidth',.5);
    hold on
end
for ii = 1:length(fastTr_PRh)
    colsel = randi(70)+90;
     HabitPRh = plot(2:3, [fastTr_PRh(ii) lastTr_PRh(ii)], '-o', 'color', cbFan_PRh(colsel,:),...
      'MarkerFaceColor', cbFan_PRh(colsel,:),'MarkerSize',2,'LineWidth',.5);
    hold on
end
ylim([-10 110]), xlim([0 3.5]);
xticks([.5,1.5,2,3]); xticklabels({'First Trials','Last Trials','First Trials', 'Last Trials'})
scatter([.5,1.5,2,3],[mean(fastTr_PCx),mean(lastTr_PCx),mean(fastTr_PRh),mean(lastTr_PRh)],80,'s','k','filled')
delta_PCx = lastTr_PCx - fastTr_PCx;
delta_PRh = lastTr_PRh - fastTr_PRh;
[h,p]= ttest(delta_PCx)
[h,p]= ttest(delta_PRh)
legend([HabitPCx,HabitPRh],{'PCx','PRh'});
legend boxoff 
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
title('Spike Response Change across Trials')
% Test significance 
ranksum(fastTr_PRh,lastTr_PRh)
ranksum(fastTr_PCx,lastTr_PCx)
text(1, 100, '***','FontSize',14);
text(2.5, 100, '***','FontSize',14);
%% Tuning Curve

load('cfg_PCx_THESIS');
load('database_PCx_THESIS');
nResponsiveEventsFS_PCx = [];
nResponsiveEventsRS_PCx = [];
normTuningCurveFS_PCx = [];
normTuningCurveRS_PCx = [];
zs_TuningCurve_PCx = [];
grad_PCx = [];

for idxExp = 1:size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        eventResponse = 0;
        tuningCurve = nan(1, size(exp.expID(idxExp).unit(idxUnit).event, 2)-1);
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2) -1
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
            tuningCurve(idxEvent) = mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse);    
        end
        if eventResponse > 0
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'FS')
                nResponsiveEventsFS_PCx = [nResponsiveEventsFS_PCx; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveFS_PCx = [normTuningCurveFS_PCx; app];
            else
                nResponsiveEventsRS_PCx = [nResponsiveEventsRS_PCx; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveRS_PCx = [normTuningCurveRS_PCx; app];
                app_grad_PCx = polyfit(1:10,fliplr(app(5:end)),1);
                grad_PCx = [grad_PCx app_grad_PCx(1)];
                zs_TuningCurve_PCx = [zs_TuningCurve_PCx; zscore(tuningCurve)];
            end
        end
    end
end

load('database_PRh_THESIS');
load('cfg_PRh_THESIS');

nResponsiveEventsFS_PRh = [];
nResponsiveEventsRS_PRh = [];
normTuningCurveFS_PRh = [];
normTuningCurveRS_PRh = [];
zs_TuningCurve_PRh = [];
grad_PRh = [];

for idxExp = 1:size(exp.expID,2) %size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        eventResponse = 0;
        tuningCurve = nan(1, size(exp.expID(idxExp).unit(idxUnit).event, 2)-1);
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2)-1
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
            tuningCurve(idxEvent) = mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse);    
        end
        if eventResponse > 0
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'FS')
                nResponsiveEventsFS_PRh = [nResponsiveEventsFS_PRh; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveFS_PRh = [normTuningCurveFS_PRh; app];
            else
                nResponsiveEventsRS_PRh = [nResponsiveEventsRS_PRh; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveRS_PRh = [normTuningCurveRS_PRh; app];
                app_grad_PRh = polyfit(1:10,fliplr(app(5:end)),1);
                grad_PRh = [grad_PRh app_grad_PRh(1)];
                zs_TuningCurve_PRh = [zs_TuningCurve_PRh; zscore(tuningCurve)];
            end
        end
    end
end

figure( 'Position', [10 10 800 500])
KK_PCx = normTuningCurveRS_PCx;
meanKK_PCx = fliplr(mean(KK_PCx));
stdKK_PCx = fliplr(std(KK_PCx));
% [maxKK, imaxKK] = max(meanKK);
errKK_PCx = stdKK_PCx./sqrt(size(KK_PCx,1));
k2 = shadedErrorBar(1:14, meanKK_PCx, errKK_PCx,'-k',1);
set(k2.edge,'LineWidth',.5,'LineStyle','none');
k2.mainLine.LineWidth = 2; k2.mainLine.Color= cb_PCx;
k2.patch.FaceAlpha = .3; k2.patch.FaceColor= cb_PCx;
hold on
KK_PCx_FS = normTuningCurveFS_PCx;
meanKK_PCx_FS = fliplr(mean(KK_PCx_FS));
stdKK_PCx_FS = fliplr(std(KK_PCx_FS));
% [maxKK, imaxKK] = max(meanKK);
errKK_PCx_FS = stdKK_PCx_FS./sqrt(size(KK_PCx_FS,1));
k2_FS = shadedErrorBar(1:14, meanKK_PCx_FS, errKK_PCx_FS,'-k',1);
set(k2_FS.edge,'LineWidth',.5,'LineStyle','none');
k2_FS.mainLine.LineWidth = 2; k2_FS.mainLine.Color= cb_PCx; k2_FS.mainLine.LineStyle= '--';
k2_FS.patch.FaceAlpha = .1;  k2_FS.patch.FaceColor= cb_PCx;

KK_PRh = normTuningCurveRS_PRh;
meanKK_PRh = fliplr(mean(KK_PRh));
stdKK_PRh = fliplr(std(KK_PRh));
% [maxKK, imaxKK] = max(meanKK);
errKK_PRh = stdKK_PRh./sqrt(size(KK_PRh,1));
k1 = shadedErrorBar(1:14, meanKK_PRh, errKK_PRh,'-k',1);
set(k1.edge,'LineWidth',.5,'LineStyle','none');
k1.mainLine.LineWidth = 2; k1.mainLine.Color= cb_PRh;
k1.patch.FaceAlpha = .3;  k1.patch.FaceColor= cb_PRh;

KK_PRh_FS = normTuningCurveFS_PRh;
meanKK_PRh_FS = fliplr(mean(KK_PRh_FS));
stdKK_PRh_FS = fliplr(std(KK_PRh_FS));
% [maxKK, imaxKK] = max(meanKK);
errKK_PRh_FS = stdKK_PRh_FS./sqrt(size(KK_PRh_FS,1));
k1_FS = shadedErrorBar(1:14, meanKK_PRh_FS, errKK_PRh_FS,'-k',1);
set(k1_FS.edge,'LineWidth',.5,'LineStyle','none');
k1_FS.mainLine.LineWidth = 2; k1_FS.mainLine.Color= cb_PRh; k1_FS.mainLine.LineStyle= '--';
k1_FS.patch.FaceAlpha = .1;  k1_FS.patch.FaceColor= cb_PRh;

hold on
plot(2:11, ones(10,1)*-0.2, 'k-', 'LineWidth', 1)
 text(6, -0.3, '***','FontSize',14);
xlabel('Ordered Stimuli based on Preference')
ylabel('Scaled Response')
legend([k2.mainLine, k2_FS.mainLine, k1.mainLine, k1_FS.mainLine],{'PCx RS','PCx FS','PRh RS','PRh FS'}, 'NumColumns',2,'Location','northeast')
legend('boxoff')
title('Scaled Average Tuning Curve')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
xticks(1:14)
xlim([0 15])

clear pp
clear hh
for ii = 1:14
[h,p]=ranksum(fliplr(normTuningCurveRS_PRh(:,ii)),fliplr(normTuningCurveRS_PCx(:,ii)))
hh(ii)=h
pp(ii)=p
end

clear pp
clear hh
for ii = 1:14
[h,p]=ranksum(fliplr(normTuningCurveFS_PRh(:,ii)),fliplr(normTuningCurveFS_PCx(:,ii)))
hh(ii)=h
pp(ii)=p
end

% figure('Position', [10 10 1800 600])
% plot(fliplr(mean(normTuningCurveFS_PRh)), '--', 'color', cb_PRh, 'LineWidth', 2)
% hold on
% plot(fliplr(mean(normTuningCurveRS_PRh)), 'color', cb_PRh, 'LineWidth', 2)
% xlabel('Number of Stimuli')
% ylabel('Scaled response')
% legend('FS','RS')
% legend('boxoff')
% title('Scaled average tuning curve')
% set(gcf,'color','white', 'PaperPositionMode', 'auto');
% xticks(1:14)
% xlim([0 15])
% set(gca, 'box', 'off', 'tickDir', 'out')
% hold on
% 
% plot(fliplr(mean(normTuningCurveFS_PCx)), '--', 'color',cb_PCx, 'LineWidth', 2)
% hold on
% plot(fliplr(mean(normTuningCurveRS_PCx)), 'color',cb_PCx, 'LineWidth', 2)
% xlabel('Ordered Stimulus Number')
% ylabel('Scaled response')
% legend('PRh FS','PRh RS', 'PCx FS','PCx RS', 'NumColumns',2,'Location','north')
% legend('boxoff')
% title('Scaled average tuning curve')
% set(gcf,'color','white', 'PaperPositionMode', 'auto');
% set(gca, 'box', 'off', 'tickDir', 'out')


%% 
nexttile
aPCx=area([3.5,6],[14 14]);
aPCx.FaceColor = cb_PRh;
aPCx.FaceAlpha = 0.3;
aPCx.LineStyle = 'none'
hold on
aPRh=area([1,3.5],[14 14]);
aPRh.FaceColor = cb_PCx;
aPRh.FaceAlpha = 0.3;
aPRh.LineStyle = 'none'

zRS =plot(4, mean(nResponsiveEventsRS_PRh), 'o', 'markersize', 8, 'markeredgecolor', 'k', 'markerfacecolor', 'k')
hold on
zFS = plot(5, mean(nResponsiveEventsFS_PRh), 'o', 'markersize', 8, 'markeredgecolor', 'k')
hold on
plot(2, mean(nResponsiveEventsRS_PCx), 'o', 'markersize', 8, 'markeredgecolor', 'k', 'markerfacecolor', 'k')
hold on
plot(3, mean(nResponsiveEventsFS_PCx), 'o', 'markersize', 8, 'markeredgecolor','k')
hold on
errbar(5, mean(nResponsiveEventsFS_PRh), std(nResponsiveEventsFS_PRh)./sqrt(length(nResponsiveEventsFS_PRh)), 'color', 'k', 'linewidth', 1); %
hold on
errbar(4, mean(nResponsiveEventsRS_PRh), std(nResponsiveEventsRS_PRh)./sqrt(length(nResponsiveEventsRS_PRh)), 'color', 'k', 'linewidth', 1); %
hold on
errbar(3, mean(nResponsiveEventsFS_PCx), std(nResponsiveEventsFS_PCx)./sqrt(length(nResponsiveEventsFS_PCx)), 'color', 'k', 'linewidth', 1); %
hold on
errbar(2, mean(nResponsiveEventsRS_PCx), std(nResponsiveEventsRS_PCx)./sqrt(length(nResponsiveEventsRS_PCx)), 'color', 'k', 'linewidth', 1); %

xlim([1 6])
xticks([2.25 4.75]);
xticklabels({'PCx', 'PRh'});
ylim([0 14])
legend([zRS,zFS],'RS','FS')
legend('boxoff')
set(gca, 'box','off')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
title('Mean Number of Responses per Neuron')
ylabel('Number of Stimuli')

ranksum(nResponsiveEventsRS_PRh,nResponsiveEventsRS_PCx)
ranksum(nResponsiveEventsFS_PRh,nResponsiveEventsFS_PCx)
%% Lifetime Sparseness
load('database_PCx_THESIS');
load('cfg_PCx_THESIS');

LS_PCx = [];
LS_PCx_RS = [];
LS_PCx_FS = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        i = 0;
         for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)
             if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse==1
                i = i + 1;
             end
         end
         if i > 0
            LS_PCx = [LS_PCx exp.expID(idxExp).unit(idxUnit).lifetimeSparseness];
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
               LS_PCx_RS = [LS_PCx_RS exp.expID(idxExp).unit(idxUnit).lifetimeSparseness];
            else
               LS_PCx_FS = [LS_PCx_FS exp.expID(idxExp).unit(idxUnit).lifetimeSparseness];
            end
         end
    end
end

load('database_PRh_THESIS');
load('cfg_PRh_THESIS');
LS_PRh = [];
LS_PRh_RS = [];
LS_PRh_FS = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        i = 0;
         for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)
             if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse==1
                i = i + 1;
             end
         end
         if i > 0
            LS_PRh = [LS_PRh exp.expID(idxExp).unit(idxUnit).lifetimeSparseness];
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
               LS_PRh_RS = [LS_PRh_RS exp.expID(idxExp).unit(idxUnit).lifetimeSparseness];
            else
               LS_PRh_FS = [LS_PRh_FS exp.expID(idxExp).unit(idxUnit).lifetimeSparseness];
            end
         end
    end
end
figure('Position', [10 10 900 300])
tiledlayout(1,2)
nexttile()
%subplot(1,2,1)
hLS_PCxRS = histogram(LS_PCx_RS, 0:0.2:1, 'Normalization','probability')
hLS_PCxRS.FaceColor = cb_PCx;
hLS_PCxRS.FaceAlpha = .8;
hLS_PCxRS.EdgeColor = [0 0 0];
hLS_PCxRS.EdgeAlpha = 0;
ylim([0 0.65])
yticks(0:0.1:0.5)
hold on
hLS_PRhRS = histogram(LS_PRh_RS, 0:0.2:1, 'Normalization','probability')
hLS_PRhRS.FaceColor = cb_PRh;
hLS_PRhRS.FaceAlpha = .8;
hLS_PRhRS.EdgeColor = [0 0 0];
hLS_PRhRS.EdgeAlpha = 0;
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
title('Regular Spiking')
scatter(LS_PCx_RS, 0.57 +0.02*rand(1,length(LS_PCx_RS)), 'filled', 'MarkerFaceColor', cb_PCx,'MarkerEdgeColor','none',...
    'SizeData',5)
scatter(LS_PRh_RS, 0.60 +0.02*rand(1,length(LS_PRh_RS)), 'filled', 'MarkerFaceColor', cb_PRh,'MarkerEdgeColor','none',...
    'SizeData',5)
zzzzlabel = ylabel('Fraction of Neurons');
xlabel('Lifetime Sparseness') 
zzzzlabel.Position(2) = 0.25;
zzzzlabel.Position(1) = -0.2;
hold on
% sb = boxplot(LS_PRh_RS, 'Orientation','horizontal','positions', 0)
% axes('Position',[0.145 0.12 .3 .10])
% set(gcf,'color','white', 'PaperPositionMode', 'auto');
% set(gca, 'box', 'off', 'tickDir', 'out')
% xlim([0 1])
 text(.45, 0.4, '***','FontSize',14);

nexttile
hLS_PCxFS = histogram(LS_PCx_FS, 0:0.2:1, 'Normalization','probability')
hLS_PCxFS.FaceColor = cb_PCx;
hLS_PCxFS.FaceAlpha = .8;
hLS_PCxFS.EdgeColor = [0 0 0];
hLS_PCxFS.EdgeAlpha = 0;
ylim([0 0.65])
yticks(0:0.1:0.5)
hold on
hLS_PRhFS = histogram(LS_PRh_FS, 0:0.2:1, 'Normalization','probability')
hLS_PRhFS.FaceColor = cb_PRh;
hLS_PRhFS.FaceAlpha = .8;
hLS_PRhFS.EdgeColor = [0 0 0];
hLS_PRhFS.EdgeAlpha = 0;
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
title('Fast Spiking')
scatter(LS_PCx_FS, 0.57 +0.02*rand(1,length(LS_PCx_FS)), 'filled', 'MarkerFaceColor', cb_PCx,'MarkerEdgeColor','none',...
    'SizeData',5)
scatter(LS_PRh_FS, 0.60 +0.02*rand(1,length(LS_PRh_FS)), 'filled', 'MarkerFaceColor', cb_PRh,'MarkerEdgeColor','none',...
    'SizeData',5)
zzzzlabel = ylabel('Fraction of Neurons');
xlabel('Lifetime Sparseness') 
zzzzlabel.Position(2) = 0.25;
zzzzlabel.Position(1) = -0.2;
hold off
% 
% nexttile
% z1= cdfplot(LS_PCx_RS)
% hold on
% z2=cdfplot(LS_PRh_RS)
% z1.Color = cb_PCx
% z2.Color = cb_PRh
% grid off
% xlabel('Lifetime Sparseness') 
% set(gcf,'color','white', 'PaperPositionMode', 'auto');
% set(gca, 'box', 'off', 'tickDir', 'out')

[p_RS,h_RS]= ranksum(LS_PRh_RS,LS_PCx_RS)
[p_RS,h_RS]= ranksum(LS_PRh_FS,LS_PCx_FS)
[p_RS,h_RS]= kstest2(LS_PRh_RS,LS_PCx_RS)
[p_FS,h_FS]= kstest2(LS_PRh_FS,LS_PCx_FS)
mean(LS_PRh_RS)
mean(LS_PCx_RS)
%%
subplot(1,2,1)
d{1} =  LS_PCx_RS; d{2} = LS_PRh_RS
h1 = raincloud_plot(d{1}, 'band_width',0.1,'box_on', 1, 'color', cb_PCx, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .25,...
    'box_col_match', 1);
h2 = raincloud_plot(d{2},'band_width',0.1, 'box_on', 1, 'color', cb_PRh, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .45, 'dot_dodge_amount', .55, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';

subplot(1,2,2)
d{1} =  LS_PCx_FS; d{2} = LS_PRh_FS
h1 = raincloud_plot(d{1}, 'band_width',0.1,'box_on', 1, 'color', cb_PCx, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .25,...
    'box_col_match', 1);
h2 = raincloud_plot(d{2},'band_width',0.1, 'box_on', 1, 'color', cb_PRh, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .45, 'dot_dodge_amount', .55, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';
%% SIGNAL & Noise Correlations


load('cfg_PCx_THESIS.mat');
load('database_PCx_THESIS.mat');
nResponsiveEventsFS_PCx = [];
nResponsiveEventsRS_PCx = [];
normTuningCurveFS_PCx = [];
normTuningCurveRS_PCx = [];
zs_TuningCurve_PCx = [];
tuningCurve = [];
for idxExp = 1:size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        eventResponse = 0;
        tuningCurve = nan(1, size(exp.expID(idxExp).unit(idxUnit).event, 2)- 1);
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2) - 1
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
            tuningCurve(idxEvent) = mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse);    
        end
        if eventResponse > 0
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'FS')
                nResponsiveEventsFS_PCx = [nResponsiveEventsFS_PCx; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveFS_PCx = [normTuningCurveFS_PCx; app];
            else
                nResponsiveEventsRS_PCx = [nResponsiveEventsRS_PCx; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveRS_PCx = [normTuningCurveRS_PCx; app];
                
            end
            zs_TuningCurve_PCx = [zs_TuningCurve_PCx; zscore(tuningCurve)];
        end
    end
end

load('database_PRh_THESIS.mat');
load('cfg_PRh_THESIS.mat');

nResponsiveEventsFS_PRh = [];
nResponsiveEventsRS_PRh = [];
normTuningCurveFS_PRh = [];
normTuningCurveRS_PRh = [];
zs_TuningCurve_PRh = [];
tuningCurve = [];
for idxExp = 1:size(exp.expID,2) %size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        eventResponse = 0;
        tuningCurve = nan(1, size(exp.expID(idxExp).unit(idxUnit).event, 2) - 1);
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2) - 1
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
            tuningCurve(idxEvent) = mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse);    
        end
        if eventResponse > 0
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'FS')
                nResponsiveEventsFS_PRh = [nResponsiveEventsFS_PRh; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveFS_PRh = [normTuningCurveFS_PRh; app];
            else
                nResponsiveEventsRS_PRh = [nResponsiveEventsRS_PRh; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveRS_PRh = [normTuningCurveRS_PRh; app];
                
            end
            zs_TuningCurve_PRh = [zs_TuningCurve_PRh; zscore(tuningCurve)];
        end
    end
end

U_PRh = corr(zs_TuningCurve_PRh');
U_PCx = corr(zs_TuningCurve_PCx');
At_PRh = U_PRh.';
m  = (1:size(At_PRh,1)).' >= (1:size(At_PRh,2));
v_PRh  = At_PRh(m);
v_PRh = v_PRh(v_PRh~=1);

    
At_PCx = U_PCx.';
m  = (1:size(At_PCx,1)).' >= (1:size(At_PCx,2));
v_PCx  = At_PCx(m);
v_PCx = v_PCx(v_PCx~=1);

% 
% for ii = 1:189
%     U_PRh_randomiz(ii,:) = zs_TuningCurve_PRh(ii, randperm(15));
% end

% U_PRh_randomiz_corr = corr(zs_TuningCurve_PRh');
% At = U_PRh_randomiz_corr.'
% m  = (1:size(At,1)).' >= (1:size(At,2));
% v_PRh_shuff  = At(m);
% v_PRh_shuff = v_PRh_shuff(v_PRh_shuff~=1)
mean(v_PCx);mean(v_PRh)
%% Noise Correlations

load('database_PCx_THESIS.mat');
load('cfg_PCx_THESIS.mat');
tempMean = [];
all_NoiseVectors_PCx = [];
eventResponse = 0;
l = 0;
for idxExp = 1:size(exp.expID,2) %size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2)-1
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
        end
        if eventResponse > 0
            l = l + 1
            %if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
                for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2) -1 
                    tempMean = [tempMean exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse - mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse)];
                end
                all_NoiseVectors_PCx(l,:) = zscore(tempMean);
                tempMean = [];
            %end
            eventResponse = 0;
        end
    end
end

load('database_PRh_Th1.mat');
load('cfg_PRh_Th.mat');
tempMean = [];
all_NoiseVectors_PRh = [];
eventResponse = 0;
l=0;
for idxExp = 1:size(exp.expID,2) %size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2) - 1
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
        end
        if eventResponse > 0
            l = l+1
           % if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
                for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2) - 1
                    tempMean = [tempMean exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse - mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse)];
                end
                all_NoiseVectors_PRh(l,:) = zscore(tempMean);
                tempMean = [];
            %end
            eventResponse = 0;
        end
    end
end

U_PRh = corr(all_NoiseVectors_PRh');
U_PCx = corr(all_NoiseVectors_PCx');
At = U_PRh.';
m  = (1:size(At,1)).' >= (1:size(At,2));
vN_PRh  = At(m);
vN_PRh = vN_PRh(vN_PRh~=1);
At = U_PCx.';
m  = (1:size(At,1)).' >= (1:size(At,2));
vN_PCx  = At(m);
vN_PCx = vN_PCx(vN_PCx~=1);

%%
% example 1
f7 = figure( 'Position', [10 10 1400 800])
tiledlayout(2,2,'TileSpacing','compact')
nexttile([2 1])
 [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
d{1} = v_PCx; d{2} = v_PRh;
h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cb_PCx, 'alpha', 0.7,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
     'box_col_match', 1);
h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cb_PRh, 'alpha', 0.7,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';
h1{1, 2}.SizeData = 1; h2{1, 2}.SizeData = 1; h1{1, 2}.MarkerFaceAlpha = .7;; h2{1, 2}.MarkerFaceAlpha = .7;
title(['Signal Correlations']);
set(gca,'XLim', [-1 1]); set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
ax = gca;   %or as appropriate
yticks = get(ax, 'YTick'); yticks=yticks(yticks>=0); set(ax, 'YTick', yticks);
hold on 
ylab = ylabel('PDF')
ylab.Position(2) = 0.75
%ylim([-1.3 1.5])
%h3 = line([0 0], [0 1.5], 'LineStyle', '--');
%h3.Color = [0 0 0 .5];

nexttile([2 1])
d{3} = vN_PCx; d{4} = vN_PRh;
h1 = raincloud_plot(d{3}, 'box_on', 1, 'color', cb_PCx, 'alpha', 0.7,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
     'box_col_match', 1);
h2 = raincloud_plot(d{4}, 'box_on', 1, 'color', cb_PRh, 'alpha', 0.7,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';
h1{1, 2}.SizeData = 1; h2{1, 2}.SizeData = 1; h1{1, 2}.MarkerFaceAlpha = .7;; h2{1, 2}.MarkerFaceAlpha = .7;
title(['Noise Correlations']);
set(gca,'XLim', [-1 1]);set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
ax = gca;   %or as appropriate
yticks = get(ax, 'YTick'); yticks=yticks(yticks>=0); set(ax, 'YTick', yticks);
%yticklabels = get(ax, 'YTickLabel');
%yticklabels{1} = ''; yticklabels{2} = ''; yticklabels{3} = ''; set(ax, 'YTickLabel', yticklabels);
legend([h1{1} h2{1}], {'PCx', 'PRh'});
%h3 = line([0 0], [0 4], 'LineStyle', '--');
%h3.Color = [0 0 0 .5];
legend boxoff
ylab = ylabel('PDF')
ylab.Position(2) = 1.85
%ylim([-3.25 3.5])
[h1,p1]= ranksum(d{1},d{2});
mean(d{1}),mean(d{2})

[h2,p2]= ttest2(d{3},d{4});
mean(d{3}),mean(d{4})

%% Dendrograms
colormap(brewermap([],'*PuBuGn'));
nexttile
imagesc(At_PCx,[-1 1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');set(gca, 'YTick', []);
ax1 = gca;            
ax1.YAxis.Visible = 'off';
ax1.XAxis.Visible = 'off';
set(gca, 'XTick', []);
axis square
title('PCx')
nexttile
imagesc(At_PRh,[-1 1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'YTick', []);
ax1 = gca;            
ax1.YAxis.Visible = 'off';
set(gca, 'XTick', []);
set(gca, 'box', 'off', 'tickDir', 'out')%set(gca, 'XAxis', []);
axis square
xlabel('Neuron ID')
title('PRh')

%%
pd_PRh = linkage(squareform(pdist(U_PRh,'correlation')));
pd_PCx = linkage(squareform(pdist(U_PCx,'correlation')));

dendro_PRh = dendrogram(pd_PRh,0)
dendro_PCx = dendrogram(pd_PCx,0)


%% PSTH analysis - Latency, Peak Intensity and 
cb = brewermap(9,'*Set3');
cb_1 = cb(9,:); cb_2= cb(7,:); cb_3= cb(4,:);
load('cfg_PCx_THESIS'); load('database_PCx_THESIS');
p = 0; indices = []; eventResponseAll = []; i = 0; KK_PRh = []; RPRh=[]; z=0; auROC_PCx = []; auROC_PCx_overall = [];
zz = 0; KK_anti = []; RPRh_anti= []; stdBSLPCx = [];

for idxExp = 1:size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        for idxEventType = 10%:length(exp.expID(idxExp).unit(idxUnit).event)
            if exp.expID(idxExp).unit(idxUnit).event(idxEventType).excitatoryResponse == 1
                z = z + 1;
                [R,t] = psth(exp.expID(idxExp).unit(idxUnit).event(idxEventType).spikesPerTrial, ...
                    0.025,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                meanBSLPCx(z) = mean(R(21:1620)); stdBSLPCx(z) = std(R(21:1620));
                RPRh = [RPRh; R];
                K = (R - meanBSLPCx(z))./stdBSLPCx(z);
                KK_PRh = [KK_PRh ; K]; [KK_peaks(z) KK_peaks_indx] = max(K(1620:end));
                KK_ts(z) = t(KK_peaks_indx + 1619);
                %KK_tso_indx = find((K(1620:end)>3),1);
                %                     if isempty(KK_tso_indx)
                %                        KK_tso(z) = nan
                %                     else
                %                     KK_tso(z) = t(KK_tso_indx + 1619);
                %             end
                %                     auROC_PCx = [auROC_PCx exp.expID(temp(1)).unit(indices).event(idxEventType).auROC];
                %                     sR_PCx(z,:) = exp.expID(temp(1)).unit(indices).event(idxEventType).spikeResponse;
            elseif exp.expID(idxExp).unit(idxUnit).event(idxEventType).inhibitoryResponse == 1
                zz = zz + 1;
                [R,t] = psth(exp.expID(idxExp).unit(idxUnit).event(idxEventType).spikesPerTrial, ...
                    0.025,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                meanBSLPCx_anti(zz) = mean(R(21:1620)); stdBSLPCx_anti(zz) = std(R(21:1620));
                RPRh_anti = [RPRh_anti; R];
                K_anti = (R - meanBSLPCx_anti(zz))./stdBSLPCx_anti(zz);
                KK_anti = [KK_anti ; K_anti];
            end
        end
    end
end
% end
%             KK_peaks_bh(x) = mean(KK_peaks);
%             KK_ts_bh(x) = mean(KK_ts);
%             KK_peaks = []; KK_ts = [];
% end
peakInf_PCx = [KK_peaks' , KK_ts'];
%figure
meanKK = mean(KK_PRh);
stdKK = std(KK_PRh);
errKK = stdKK./sqrt(size(KK_PRh,1));
l1 = shadedErrorBar(t, meanKK, errKK,'transparent',1);
set(l1.edge,'LineWidth',.5,'LineStyle','none');
l1.mainLine.LineWidth = 1; l1.mainLine.Color = cb_3;%k1.mainLine.LineStyle = '--'
    l1.patch.FaceAlpha = .3; l1.patch.FaceColor = cb_3;
hold on
% meanKK_anti = mean(KK_anti);
% stdKK_anti = std(KK_anti);
% errKK_anti = stdKK_anti./sqrt(size(KK_anti,1));
% k2 = shadedErrorBar(t, meanKK_anti, errKK_anti,'transparent',1);
% set(k2.edge,'LineWidth',.5,'LineStyle','none');
% k2.mainLine.LineWidth = 2;k2.patch.FaceColor = cb_1; k2.mainLine.LineStyle = '--'
% k2.patch.FaceAlpha = .3;k2.mainLine.Color = cb_1;

hold on
M = -1*ones(1,3);
ax = plot(0:2, M, 'k', 'LineWidth', 1.5);
ax.Clipping = 'off';
hold on
%legend([k1.mainLine, k2.mainLine, z.mainLine, z2.mainLine],'PCx sign', 'PCx ns', 'PRh sign', 'PRh ns','NumColumns',2,'Location', 'northwest');
%legend('boxoff')
title('PCx');
xlabel('Time (s - Onset Centered)');
ylabel('SD over baseline');
xlim([-1 4]); ylim([-1.5 12]);
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out');
%% all
% PCx

%clear all
load('cfg_PCx_THESIS'); load('database_PCx_THESIS');
p = 0; indices = []; eventResponseAll = []; i = 0; KK_PRh = []; RPRh=[]; z=0; auROC_PCx = []; auROC_PCx_overall = [];
zz = 0; KK_anti = []; RPRh_anti= []; stdBSLPCx = []; KK_tso = [];

for idxExp = 1:size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        for idxEventType = 1:length(exp.expID(idxExp).unit(idxUnit).event)-1
            if exp.expID(idxExp).unit(idxUnit).event(idxEventType).excitatoryResponse == 1
                z = z + 1;
                [R,t] = psth(exp.expID(idxExp).unit(idxUnit).event(idxEventType).spikesPerTrial, ...
                    0.025,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                meanBSLPCx(z) = mean(R(21:1620)); stdBSLPCx(z) = std(R(21:1620));
                RPRh = [RPRh; R];
                K = (R - meanBSLPCx(z))./stdBSLPCx(z);
                KK_PRh = [KK_PRh ; K]; [KK_peaks(z) KK_peaks_indx(z)] = max(K(1621:end));
                %KK_ts(z) = t(KK_peaks_indx + 1619);
                KK_tso_indx = find((K(1621:end)> 2),1);
                if isempty(KK_tso_indx)
                    KK_tso(z) = nan
                else
                    KK_tso(z) = t(KK_tso_indx + 1620);
                end
                %                     auROC_PCx = [auROC_PCx exp.expID(temp(1)).unit(indices).event(idxEventType).auROC];
                %                     sR_PCx(z,:) = exp.expID(temp(1)).unit(indices).event(idxEventType).spikeResponse;
            elseif exp.expID(idxExp).unit(idxUnit).event(idxEventType).inhibitoryResponse == 1
                zz = zz + 1;
                [R,t] = psth(exp.expID(idxExp).unit(idxUnit).event(idxEventType).spikesPerTrial, ...
                    0.025,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                meanBSLPCx_anti(zz) = mean(R(21:1620)); stdBSLPCx_anti(zz) = std(R(21:1620));
                RPRh_anti = [RPRh_anti; R];
                K_anti = (R - meanBSLPCx_anti(zz))./stdBSLPCx_anti(zz);
                KK_anti = [KK_anti ; K_anti];
            end
        end
    end
end
% end
%             KK_peaks_bh(x) = mean(KK_peaks);
%             KK_ts_bh(x) = mean(KK_ts);
%             KK_peaks = []; KK_ts = [];
% end
%peakInf_PCx = [KK_peaks' , KK_ts'];
%figure( 'Position', [10 10 1200 600])
figure
set(gcf, 'PaperPosition', [0 0 15 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 5]); %Set the paper to have width 5 and height 5.
 %Save figure
subplot(1,3,[1 2])
meanKK = mean(KK_PRh);
stdKK = std(KK_PRh);
errKK = stdKK./sqrt(size(KK_PRh,1));
k1 = shadedErrorBar(t, meanKK, errKK);
set(k1.edge,'LineWidth',.5,'LineStyle','none');
k1.mainLine.LineWidth = 1.5; k1.mainLine.Color = cb_PCx;
    k1.patch.FaceAlpha = .4; k1.patch.FaceColor = cb_PCx;
hold on
meanKK_anti = mean(KK_anti);
stdKK_anti = std(KK_anti);
errKK_anti = stdKK_anti./sqrt(size(KK_anti,1));
k2 = shadedErrorBar(t, meanKK_anti, errKK_anti);
set(k2.edge,'LineWidth',.5,'LineStyle','none');
k2.mainLine.LineWidth = 1.5;k2.mainLine.Color = cb_PCx;k2.mainLine.LineStyle = '--';
k2.patch.FaceAlpha = .4;k2.patch.FaceColor = cb_PCx;

% PRh

load('cfg_PRh_THESIS'); load('database_PRh_THESIS');
p = 0; indices = []; eventResponseAll = []; i = 0; ZZ = []; RPRh=[]; v_PRh=0; auROC_PRh = []; auROC_PRh_overall = [];
vv = 0; RPRh_anti = []; ZZ_anti = []; ZZ_peaks = []; ZZ_peaks_indx = []; ZZ_tso =[]; 

for idxExp = 1:size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        for idxEventType = 1:length(exp.expID(idxExp).unit(idxUnit).event) -1
            if exp.expID(idxExp).unit(idxUnit).event(idxEventType).excitatoryResponse == 1
                v_PRh = v_PRh + 1;
                [R,t] = psth(exp.expID(idxExp).unit(idxUnit).event(idxEventType).spikesPerTrial, ...
                    0.025,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                meanBSLPRh(v_PRh) = mean(R(21:1620)); stdBSLPRh(v_PRh) = std(R(21:1620));
                RPRh = [RPRh; R];
                Z = (R - meanBSLPRh(v_PRh))./stdBSLPRh(v_PRh);
                ZZ = [ZZ ; Z]; [ZZ_peaks(v_PRh) ZZ_peaks_indx] = max(Z(1621:end));
                %ZZ_ts(v_PRh) = t(ZZ_peaks_indx + 1619);
                ZZ_tso_indx = find((Z(1621:end)> 2),1);
                if isempty(ZZ_tso_indx)
                    ZZ_tso(v_PRh) = nan
                else
                    ZZ_tso(v_PRh) = t(ZZ_tso_indx + 1620);
                end
                %                     auROC_PRh = [auROC_PRh exp.expID(temp(1)).unit(indices).event(idxEventType).auROC];
                %                     sR_PRh(v_PRh,:) = exp.expID(temp(1)).unit(indices).event(idxEventType).spikeResponse;
            elseif exp.expID(idxExp).unit(idxUnit).event(idxEventType).inhibitoryResponse == 1
                [R,t] = psth(exp.expID(idxExp).unit(idxUnit).event(idxEventType).spikesPerTrial, ...
                    0.025,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                vv = vv + 1;
                meanBSLPRh_anti(vv) = mean(R(21:1620)); stdBSLPRh_anti(vv) = std(R(21:1620));
                RPRh_anti = [RPRh_anti; R];
                Z_anti = (R - meanBSLPRh_anti(vv))./stdBSLPRh_anti(vv);
                ZZ_anti = [ZZ_anti ; Z_anti];
            end
        end
    end
end

%             ZZ_peaks_bh(x) = mean(ZZ_peaks);
%             ZZ_ts_bh(x) = mean(ZZ_ts);
%             ZZ_peaks = []; ZZ_ts = [];
% end
% peakInf_PRh = [ZZ_peaks' , ZZ_ts'];
% set(gcf,'color','white', 'PaperPositionMode', 'auto');
% set(gca, 'box', 'off', 'tickDir', 'out');
meanZZ = mean(ZZ);
stdZZ = std(ZZ);
errZZ = stdZZ./sqrt(size(ZZ,1));
z = shadedErrorBar(t, meanZZ, errZZ);
set(z.edge,'LineWidth',.5,'LineStyle','none');
z.mainLine.LineWidth = 1.5; z.mainLine.Color = cb_PRh;
    z.patch.FaceAlpha = .4; z.patch.FaceColor = cb_PRh;
meanZZ_anti = nanmean(ZZ_anti);
stdZZ_anti = nanstd(ZZ_anti);
errZZ_anti = stdZZ_anti./sqrt(size(ZZ_anti,1));
z2 = shadedErrorBar(t, meanZZ_anti, errZZ_anti);
set(z2.edge,'LineWidth',.5,'LineStyle','none');
z2.mainLine.LineWidth = 1.5; z2.mainLine.Color = cb_PRh; z2.mainLine.LineStyle = '--';
z2.patch.FaceAlpha = .4; z2.patch.FaceColor = cb_PRh;
hold on
title('Grand-Average PSTH');
xlabel('Time (s - Onset Aligned)');
ylabel('SD over baseline');
xlim([-1 4]); ylim([-2.5 8.5]);
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out');
M = -2*ones(1,3);
ax = plot(0:2, M, 'k', 'LineWidth', 1.5);
ax.Clipping = 'off';
hold on
legend([k1.mainLine, k2.mainLine, z.mainLine, z2.mainLine],'PCx Excitatory', 'PCx Inhibitory', 'PRh Excitatory', 'PRh Inhibitory','NumColumns',2,'Location', 'north');
legend('boxoff')

axes('Position',[.45 .55 .15 .20])
plot(t,meanKK./max(meanKK), 'Color',cb_PCx, 'LineWidth', 1.5);
hold on
plot(t,meanZZ./max(meanZZ), 'Color',cb_PRh,'LineWidth', 1.5);
M = -0.1*ones(1,3);
ax = plot(0:2, M, 'k', 'LineWidth', 1.5);
xlim([-1 4]); ylim([-0.2 1.2]);
title('Peak-Normalized Excitatory PSTH');
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out');

%%
figure( 'Position', [10 10 1000 800])
meanKK = mean(KK_PRh)./max(mean(KK_PRh));
stdKK = std(KK_PRh);
errKK = stdKK./sqrt(size(KK_PRh,1));
k1 = shadedErrorBar(t, meanKK, zeros(length(t),1),'transparent',1);
set(k1.edge,'LineWidth',.5,'LineStyle','none');
k1.mainLine.LineWidth = 1.5; k1.mainLine.Color = cb_PCx;
    k1.patch.FaceAlpha = .3; k1.patch.FaceColor = cb_PCx;
hold on

meanZZ = mean(ZZ)./max(mean(ZZ));
stdZZ = std(ZZ);
errZZ = stdZZ./sqrt(size(ZZ,1));
z = shadedErrorBar(t, meanZZ, zeros(length(t),1),'transparent',1);
set(z.edge,'LineWidth',.5,'LineStyle','none');
z.mainLine.LineWidth = 1.5; z.mainLine.Color = cb_PRh;
    z.patch.FaceAlpha = .3; z.patch.FaceColor = cb_PRh;

hold on
M = -.05*ones(1,3);
ax = plot(0:2, M, 'k', 'LineWidth', 1);
ax.Clipping = 'off';
hold on
legend([k1.mainLine, z.mainLine],'PCx Excitatory','PRh Excitatory','Location', 'northwest');
legend('boxoff')
title('Grand-Averaged PSTH Peak-Normalized');
xlabel('Time (s - Onset Aligned)');
ylabel('Scaled Response');
xlim([-1 4]); ylim([-.1 1.4]); yticks([0 1]);
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out');

%% Peak Latency

d{1} =  KK_tso; d{2} = ZZ_tso
subplot(1,3,3)
%figure( 'Position', [10 10 600 800])
h1 = raincloud_plot(d{1},'box_on', 1, 'color', cb_PCx, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
    'box_col_match', 1);
h2 = raincloud_plot(d{2},'box_on', 1, 'color', cb_PRh, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .60, 'dot_dodge_amount', .80, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';
h1{2}.SizeData = 5; h2{2}.SizeData = 5;
ylim([-6 6.5]); xlim([0 1]);title('Onset Latencies Distributions (sd $>$ 2)');
ylab= ylabel('PDF'); xlabel('Time after Stimulus Onset (s)');
legend([h1{1},h2{1}],'PCx','PRh'); yticks([0:2:6]);
legend boxoff 
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
text(0.13, 5.7, '***','FontSize',14);
ylab.Position(2) = 3;
[ppl,hpl]=ranksum(KK_tso,ZZ_tso)
[hpl,ppl]=ttest2(KK_tso,ZZ_tso)
[hpl,ppl]=kstest2(KK_tso,ZZ_tso)

%% Odor Set differences

nNeuroDec = [0:10:200];
nNeuroDec(1) = 1;
% PCx
for i = 1:21
    stringa = ['database_PRh_RSEXCINH_SetA__' num2str(nNeuroDec(i)) '_neurons.mat']
    load(stringa)
    performancesPCx_mu(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(2)];
    performancesPCx_stdev(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(2)];
    performancesPCx_mu_succ(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(3)];
    performancesPCx_stdev_succ(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(3)];
end

f7 = figure( 'Position', [10 10 800 400])
%tiledlayout(2,7)
%nexttile([2 3])
numOfStimuli = 7;
decSetPCx = shadedErrorBar(nNeuroDec, performancesPCx_mu, performancesPCx_stdev, 'k',1)%,{'color', 'r', 'linewidth', 2});
set(decSetPCx .edge,'LineWidth',.5,'LineStyle','none');
decSetPCx.mainLine.LineWidth = 2; decSetPCx.mainLine.Color= cb_setA;
decSetPCx.patch.FaceAlpha = .3; decSetPCx.patch.FaceColor=cb_setA;


line([0 nNeuroDec(end)+5], [1/numOfStimuli 1/numOfStimuli], 'color', [.3 .3 .3], 'linestyle', '--')
text(172,1/numOfStimuli+0.03,'Chance Level','color', [.3 .3 .3]);
set(gcf,'color','white', 'PaperPositionMode', 'auto');
ylim([0 1])
xlim([0 nNeuroDec(end)+1]); xticks(nNeuroDec);
ylabel('Accuracy \%')
xlabel('Population Size')
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
title('Mean Decoding Accuracy of Odor Identity')
hold on

% PRh

for i = 1:21
    stringa = ['database_PRh_RSEXCINH_SetB__' num2str(nNeuroDec(i)) '_neurons.mat']
    load(stringa)
    performancesPRh_mu(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(2)];
    performancesPRh_stdev(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(2)];

end

decSetPRh = shadedErrorBar(nNeuroDec, performancesPRh_mu, performancesPRh_stdev, 'k',1)%,{'color', 'r', 'linewidth', 2});
set(decSetPRh .edge,'LineWidth',.5,'LineStyle','none');
decSetPRh.mainLine.LineWidth = 2; decSetPRh.mainLine.Color= cb_setB;
decSetPRh.patch.FaceAlpha = .3; decSetPRh.patch.FaceColor= cb_setB;

legend([decSetPCx.mainLine,decSetPRh.mainLine] ,{'Set A','Set B'},'Location','north','NumColumns',2)
legend boxoff

%% CM sets
load('database_PRh_RSEXCINH_SetA__200_neurons')
f7 = figure( 'Position', [10 10 1000 800])
set(gcf,'color','white', 'PaperPositionMode', 'auto');

periods{1} = 'Interval: [-1 0]';
periods{2} = 'Interval: [0 1]';
periods{3} = 'Interval: [1 2]';
periods{4} = 'Interval: [2 3]';
colormap(brewermap([],'*YlGnBu'));

for idxPeriod = 1:4
    CM = squeeze(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix(:,:,idxPeriod));
%     nexttile(brah(idxPeriod-1))
     subplot(3,4,idxPeriod)
    clims = [0 1];
    imagesc(CM ./ sum(CM(:,1)), clims)
    axis square
    xlabel('Predicted Stimulus')
    ylabel('Actual Stimulus')
    title(periods{idxPeriod})
end
%sgtitle('PCx')

load('database_PRh_RSEXCINH_SetB__200_neurons')
brah = [11 12];
for idxPeriod = 1:4
    CM = squeeze(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix(:,:,idxPeriod));
%     nexttile(brah(idxPeriod-1))
     subplot(3,4,idxPeriod+4)
    clims = [0 1];
    imagesc(CM ./ sum(CM(:,1)), clims)
    axis square
    xlabel('Predicted Stimulus')
    ylabel('Actual Stimulus')
    title(periods{idxPeriod})
end
h=colorbar;
set(h, 'Location','south','Position', [.2 .29 .6 0.05])
%% Decoder
nNeuroDec = [0:10:110];
nNeuroDec(1) = 1;
% PCx
for i = 1:12
    stringa = ['database_PCx_RSEXC_' num2str(nNeuroDec(i)) '_neurons.mat']
    load(stringa)
    performancesPCx_mu(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(2)];
    performancesPCx_stdev(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(2)];
    performancesPCx_mu_succ(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(3)];
    performancesPCx_stdev_succ(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(3)];
end

f7 = figure( 'Position', [10 10 800 400])
%tiledlayout(2,7)
%nexttile([2 3])
numOfStimuli = 14;
decSetPCx = shadedErrorBar(nNeuroDec, performancesPCx_mu, performancesPCx_stdev, 'k',1)%,{'color', 'r', 'linewidth', 2});
set(decSetPCx .edge,'LineWidth',.5,'LineStyle','none');
decSetPCx.mainLine.LineWidth = 2; decSetPCx.mainLine.Color= cb_PCx;
decSetPCx.patch.FaceAlpha = .3; decSetPCx.patch.FaceColor= cb_PCx;

hold on
decSetPCx_succ = shadedErrorBar(nNeuroDec, performancesPCx_mu_succ, performancesPCx_stdev_succ, 'k',1)%,{'color', 'r', 'linewidth', 2});
set(decSetPCx_succ.edge,'LineWidth',.5,'LineStyle','none');
decSetPCx_succ.mainLine.LineWidth = 2; decSetPCx_succ.mainLine.Color= cb_PCx; decSetPCx_succ.mainLine.LineStyle= '--';
decSetPCx_succ.patch.FaceAlpha = .3; decSetPCx_succ.patch.FaceColor= cb_PCx;

line([0 nNeuroDec(end)+5], [1/numOfStimuli 1/numOfStimuli], 'color', [.3 .3 .3], 'linestyle', '--')
text(96,1/numOfStimuli+0.03,'Chance Level','color', [.3 .3 .3]);
set(gcf,'color','white', 'PaperPositionMode', 'auto');
ylim([0 1])
xlim([0 nNeuroDec(end)+1]); xticks(nNeuroDec);
ylabel('Accuracy \%')
xlabel('Population Size')
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
title('Mean Decoding Accuracy of Odor Identity')
hold on

% PRh

for i = 1:12
    stringa = ['database_PRh_RSEXC_' num2str(nNeuroDec(i)) '_neurons.mat']
    load(stringa)
    performancesPRh_mu(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(2)];
    performancesPRh_stdev(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(2)];
    performancesPRh_mu_succ(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(3)];
    performancesPRh_stdev_succ(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(3)];
end

decSetPRh = shadedErrorBar(nNeuroDec, performancesPRh_mu, performancesPRh_stdev, 'k',1)%,{'color', 'r', 'linewidth', 2});
set(decSetPRh .edge,'LineWidth',.5,'LineStyle','none');
decSetPRh.mainLine.LineWidth = 2; decSetPRh.mainLine.Color= cb_PRh;
decSetPRh.patch.FaceAlpha = .3; decSetPRh.patch.FaceColor= cb_PRh;

decSetPRh_succ = shadedErrorBar(nNeuroDec, performancesPRh_mu_succ, performancesPRh_stdev_succ, 'k',1)%,{'color', 'r', 'linewidth', 2});
set(decSetPRh_succ.edge,'LineWidth',.5,'LineStyle','none');
decSetPRh_succ.mainLine.LineWidth = 2; decSetPRh_succ.mainLine.Color= cb_PRh; decSetPRh_succ.mainLine.LineStyle = '--';
decSetPRh_succ.patch.FaceAlpha = .3; decSetPRh_succ.patch.FaceColor= cb_PRh;

legend([decSetPCx.mainLine,decSetPRh.mainLine,decSetPCx_succ.mainLine,decSetPRh_succ.mainLine] ,{'0-1s','0-1s','1-2s - PCx','1-2s - PRh'},'Location','northwest','NumColumns',2)
legend boxoff

%% CM

load('database_PCx_RSEXC_110_neurons')
f7 = figure( 'Position', [10 10 1000 800])
set(gcf,'color','white', 'PaperPositionMode', 'auto');

periods{1} = 'Interval: [-1 0]';
periods{2} = 'Interval: [0 1]';
periods{3} = 'Interval: [1 2]';
periods{4} = 'Interval: [2 3]';
colormap(brewermap([],'*YlGnBu'));

for idxPeriod = 1:4
    CM = squeeze(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix(:,:,idxPeriod));
%     nexttile(brah(idxPeriod-1))
     subplot(3,4,idxPeriod)
    clims = [0 1];
    imagesc(CM ./ sum(CM(:,1)), clims)
    axis square
    xlabel('Predicted Stimulus')
    ylabel('Actual Stimulus')
    title(periods{idxPeriod})
end
%sgtitle('PCx')

load('database_PRh_RSEXC_110_neurons')
brah = [11 12];
for idxPeriod = 1:4
    CM = squeeze(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix(:,:,idxPeriod));
%     nexttile(brah(idxPeriod-1))
     subplot(3,4,idxPeriod+4)
    clims = [0 1];
    imagesc(CM ./ sum(CM(:,1)), clims)
    axis square
    xlabel('Predicted Stimulus')
    ylabel('Actual Stimulus')
    title(periods{idxPeriod})
end
h=colorbar;
set(h, 'Location','south','Position', [.2 .29 .6 0.05])

%% Information Neurons
figure( 'Position', [10 10 900 300])
%nexttile(brah(idxPeriod-1))
neuron_information_PCx=load('SNI_PCx.mat')
neuron_information_PCx =neuron_information_PCx.neuron_information;
neuron_information_PRh=load('SNI_PRh.mat')
neuron_information_PRh =neuron_information_PRh.neuron_information;
neuron_information_All = [neuron_information_PCx; neuron_information_PRh]
Area = cellstr(["PCx","PRh"])
K    = cell(1, 118);
K(:) = {'PCx'};
C    = cell(1, 169);
C(:) = {'PRh'};
P = [K C]
% hold on
% scatter(2+randn(1,length(neuron_information_PCx))/10 ,neuron_information_PCx, 12,'filled','MarkerEdgeColor', cb_PCx,'MarkerFaceColor',cb_PCx)
% scatter(2,mean(neuron_information_PCx), 's','k','filled')
vs = violinplot(neuron_information_All,P)
vs(1, 1).ViolinColor = cb_PCx;vs(1, 2).EdgeColor = 'none';vs(1, 1).ViolinAlpha = .2;vs(1, 1).ScatterPlot.MarkerFaceAlpha =.5
vs(1, 2).ViolinColor = cb_PRh;vs(1, 1).EdgeColor = 'none';vs(1, 2).ViolinAlpha = .2;vs(1, 2).ScatterPlot.MarkerFaceAlpha =.5
xlim([.5 2.5])
ylim([0 2])
line([1 2], [1.6 1.6], 'color', 'k', 'linestyle', '-')
% scatter(5+randn(1,length(neuron_information_PRh))/10 ,neuron_information_PRh, 12, 'filled','MarkerEdgeColor', cb_PRh,'MarkerFaceColor',cb_PRh);
% scatter(5,mean(neuron_information_PRh), 30,'s','k','filled')

ylabel('bits')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
 text(1.45,1.7, '***','FontSize',14);
 title('Single Unit Mutual Information')
[pni,hni]=ranksum(neuron_information_PCx,neuron_information_PRh)


%% Decoding TimeCourse
load('PCx_TimeCourse_200by100ms')
figure( 'Position', [10 10 800 400])
numOfStimuli = 14;
tcperf_mu_PCx = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_PCx = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
k = shadedErrorBar(1:49, tcperf_mu_PCx, tcperf_stdev_PCx,'k',1);
hold on
set(k.edge,'LineWidth',.5,'LineStyle','none');
k.mainLine.LineWidth = 2; k.mainLine.Color= cb_PCx;
k.patch.FaceAlpha = .3; k.patch.FaceColor= cb_PCx;
load('PRh_TimeCourse_200by100ms')
tcperf_mu_PRh = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_PRh = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
z = shadedErrorBar(1:49, tcperf_mu_PRh, tcperf_stdev_PRh,'k',1);
hold on
set(z.edge,'LineWidth',.5,'LineStyle','none');
z.mainLine.LineWidth = 2; z.mainLine.Color= cb_PRh;
z.patch.FaceAlpha = .3; z.patch.FaceColor= cb_PRh;

line([0 49], [1/numOfStimuli 1/numOfStimuli], 'color', [.3 .3 .3], 'linestyle', '--');
text(2,1/numOfStimuli+0.06,'Chance Level','color', [.3 .3 .3]);
M = -.0*ones(1,21);
ax = plot(9:29, M, 'k', 'LineWidth', 2);
ax.Clipping = 'off';
hold on
ylabel('Accuracy \%')
title('Decoding over Timecourse')
ylim([0 1])
xlim([1 49]); xticks(9:10:49); xticklabels(0:4)
xlabel('Time (s - Onset Aligned)')

legend([k.mainLine,z.mainLine] ,{'PCx','PRh'},'Location','northwest')
legend boxoff
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')

%%
load('PCx_-1_4_TimeCourse_200by100ms')
f7 = figure( 'Position', [10 10 1400 1200])
numOfStimuli = 14;
tcperf_mu_PCx = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_PCx = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
k = shadedErrorBar(1:49, tcperf_mu_PCx./max(tcperf_mu_PCx), zeros(1,length(tcperf_stdev_PCx)),'k',1);
hold on
set(k.edge,'LineWidth',.5,'LineStyle','none');
k.mainLine.LineWidth = 2; k.mainLine.Color= cb_PCx;
k.patch.FaceAlpha = .3; k.patch.FaceColor= cb_PCx;
load('PRh_-1_4_TimeCourse_200by100ms')
tcperf_mu_PRh = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_PRh = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
z = shadedErrorBar(1:49, tcperf_mu_PRh./max(tcperf_mu_PRh), zeros(1,length(tcperf_stdev_PRh)),'k',1);
hold on
set(z.edge,'LineWidth',.5,'LineStyle','none');
z.mainLine.LineWidth = 2; z.mainLine.Color= cb_PRh;
z.patch.FaceAlpha = .3; z.patch.FaceColor= cb_PRh;
M = .015*ones(1,20);
ax = plot(10:29, M, 'k', 'LineWidth', 1.5);
ax.Clipping = 'off';
hold on
ylabel('Scaled Accuracy')
title('Decoding Over Timecourse - Peak Normalized - 200ms Bin $/$ 100ms Step')
ylim([0 1.2]); yticks([0:0.5:1]);
xlim([1 49]); xticks([]); 
ax1 = gca
ax1.XAxis.Visible = 'off';

legend([k.mainLine,z.mainLine] ,{'PCx','PRh'},'Location','northwest')
legend boxoff
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')

%% Rappresentazione dei TOP25 tra i topi.

% PCx

load('sortedBest20_PCx_svmtot_0-2');
subplot(2,2,1)
b1 = bar(histcounts(sortedBest20(1:15,1),[0.5:1:8.5]), 'b'); ylim([0 10]); title('Best 15 Neurons in PCx'); xlabel('Mouse ID'); ylabel('# of neurons');
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out');
b1.FaceAlpha = 0.6;

% PRh

subplot(2,2,2)
load('sortedBest20_PRh_svmtot_0-2');
b2 = bar(histcounts(sortedBest20(1:15,1)',[0.5:1:10.5]), 'r'); ylim([0 10]); title('Best 15 Neurons in PRh'); xlabel('Mouse ID'); ylabel('# of neurons');
b2.FaceAlpha = 0.6;
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out');


load('sortedBest20_PCx_svm_0-2');
subplot(2,2,3)
b1 = bar(histcounts(sortedBest20(1:25,1),[0.5:1:8.5]), 'b'); ylim([0 10]); title('Best 25 Neurons in PCx'); xlabel('Mouse ID'); ylabel('# of neurons');
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out');
b1.FaceAlpha = 0.6;

% PRh

subplot(2,2,4)
load('sortedBest20_PRh_svm_0-2');
b2 = bar(histcounts(sortedBest20(1:25,1),[0.5:1:10.5]), 'r'); ylim([0 10]); title('Best 25 Neurons in PRh'); xlabel('Mouse ID'); ylabel('# of neurons');
b2.FaceAlpha = 0.6;
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out');

%% Lifetime Sparsenees

edges = [0:0.2:1];

% PCx

subplot(1,2,1);
load('sortedBest20_PCx_svm_0-2');
histogram(sortedBest20(1:15,3) ,edges, 'FaceColor', 'b', 'EdgeColor','k', 'Normalization', 'probability')
hold on

% PRh

load('sortedBest20_PRh_svm_0-2');
histogram(sortedBest20(1:15,3), edges, 'FaceColor', 'r', 'EdgeColor','k', 'Normalization', 'probability');
xlabel('Lifetime Sparseness')
ylabel('Fraction of Neurons')
ylim([0 1]);
legend('PCx','PRh'); legend('boxoff');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
title('Best 15 Units');

subplot(1,2,2)
load('sortedBest20_PCx_svm_0-2');
histogram(sortedBest20(1:25,3) ,edges, 'FaceColor', 'b', 'EdgeColor','k', 'Normalization', 'probability')
hold on

% PRh

load('sortedBest20_PRh_svm_0-2');
histogram(sortedBest20(1:25,3), edges, 'FaceColor', 'r', 'EdgeColor','k', 'Normalization', 'probability');
xlabel('Lifetime Sparseness')
ylabel('Fraction of Neurons')
ylim([0 1]);
legend('PCx','PRh'); legend('boxoff');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
title('Best 25 Units');


%% LS for best neurons

for i=0:19
    load('sortedBest20_PCx_svmtot_0-2');
    A_PCx(i+1) = mean(sortedBest20(1:5*(i+1),3));
    Astd_PCx(i+1) = std(sortedBest20(1:5*(i+1),3))
    Aerr_PCx(i+1) = Astd_PCx(i+1)./sqrt(size(A_PCx,2)*5);
    load('sortedBest20_PRh_svmtot_0-2');
    A_PRh(i+1) = mean(sortedBest20(1:5*(i+1),3));
    Astd_PRh(i+1) = std(sortedBest20(1:5*(i+1),3))
    Aerr_PRh(i+1) = Astd_PRh(i+1)./sqrt(size(A_PRh,2)*5);
end
yy1 = shadedErrorBar(1:20,A_PCx,Aerr_PCx,'b');
hold on
yy2 = shadedErrorBar(1:20,A_PRh,Aerr_PRh,'r');
set(yy1.edge,'LineWidth',.5,'LineStyle','none');
yy1.mainLine.LineWidth = 2;
yy1.patch.FaceAlpha = .95;
set(yy2.edge,'LineWidth',.5,'LineStyle','none');
yy2.mainLine.LineWidth = 2;
yy2.patch.FaceAlpha = .95;
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out'); labelsss = 1:5:96;
xlabel('Top N neurons included'); xticks(1:19); xticklabels({'5','10','15','20','25','30','35','40','45','50','55','60','65','70','75','80','85','90','95','100'});
ylabel('Lifetime Sparseness'); xlim([1,19]); ylim([0.3 1]); xline(3, 'k--');
title('Lifetime Sparseness including Top N Neurons');
legend([yy1.mainLine, yy2.mainLine], 'PCx','PRh');
legend('boxoff');

%% Rappresentazione del numero di odori responsivi

load('sortedBest20_PCx_svmtot_0-2');
sortPCx = sortedBest20(1:15,4)
numOdor(:,1) = histcounts(sortedBest20(1:15,4),[0.5:1:14.5]);
load('sortedBest20_PRh_svmtot_0-2');
sortPRh = sortedBest20(1:15,4);
numOdor(:,2) = histcounts(sortedBest20(1:15,4),[0.5:1:14.5]);
sortAll = [sortPCx sortPRh];
% Histogram

figure
subplot(2,1,1)
c = bar(numOdor,'FaceColor','flat','BarWidth', 1);
c(1).FaceColor = [0 0 1]; c(2).FaceColor = [1 0 0]; c(1).FaceAlpha = .6; c(2).FaceAlpha = .6; ylim([0 6]);
legend('PCx','PRh'); legend('boxoff');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
xlabel('Number of Odors'); ylabel('# of neurons'); title('Number of Excitatory Responses for the Odor Set');
mean(numOdor);

% Boxplot

subplot(2,1,2)
boxplot(sortAll,'Labels',{'PCx','PRh'},'Orientation','horizontal');

%% Rappresentazione dei singoli odori

clear all
load('cfg_PCx_dm'); load('database_PCx_dm'); load('sortedBest20_PCx_svmtot_0-2');
p = 0; indices = []; v_PCx = zeros(25,14);


for x = 1:15
    temp = sortedBest20(x,:);
    for idxUnit = 1:size(exp.expID(temp(1)).unit, 2)
        if temp(2) == exp.expID(temp(1)).unit(idxUnit).ID
            indices = idxUnit;
            for idxEventType = 1:length(exp.expID(temp(1)).unit(indices).event)
                if exp.expID(temp(1)).unit(idxUnit).event(idxEventType).excitatoryResponse == 1
                    v_PCx(x,idxEventType) = 1;
                end
            end
        end
    end
end

% PRh

load('cfg_PRh_dm'); load('database_PRh_dm');  load('sortedBest20_PRh_svmtot_0-2');
p = 0; indices = []; v_PRh = zeros(25,14);

for x = 1:15
    temp = sortedBest20(x,:);
    for idxUnit = 1:size(exp.expID(temp(1)).unit, 2)
        if temp(2) == exp.expID(temp(1)).unit(idxUnit).ID
            indices = idxUnit;
            for idxEventType = 1:length(exp.expID(temp(1)).unit(indices).event)
                if exp.expID(temp(1)).unit(idxUnit).event(idxEventType).excitatoryResponse == 1
                    v_PRh(x,idxEventType) = 1;
                end
            end
        end
    end
end

figure
represOdor(:,1) = sum(v_PCx)/15;
represOdor(:,2) = sum(v_PRh)/15;
d = bar(represOdor,'FaceColor','flat','BarWidth', 1);
d(1).FaceColor = [0 0 1]; d(2).FaceColor = [1 0 0]; d(1).FaceAlpha = .6; d(2).FaceAlpha = .6;
legend('PCx','PRh'); legend('boxoff');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
xlabel('Odor ID'); ylabel('% of Neurons'); title('Percentage of significant responses for each Odor');

%% svm odors


clear all
load('cfg_PCx_dm'); load('database_PCx_dm'); load('sortedBest20_PCx_svmtot_0-2');
p = 0; indices = []; v_PCx = zeros(25,14);
% A = squeeze(mean(mean(squeeze(mean(DECODING_RESULTS.FP_INFO{1, 2}.the_p_values_org_order(:,:,2:3,:),3)),1),2));
% [output,indx] = sort(A);
% sortedBestAll = [binned_site_info.sessionID(indx) binned_site_info.clusterID(indx)];

for x = 1:15
    temp = sortedBest20(x,:);
    for idxUnit = 1:size(exp.expID(temp(1)).unit, 2)
        if temp(2) == exp.expID(temp(1)).unit(idxUnit).ID
            indices = idxUnit;
            for idxEventType = 1:length(exp.expID(temp(1)).unit(indices).event)
                if exp.expID(temp(1)).unit(idxUnit).event(idxEventType).excitatoryResponse == 1
                    v_PCx(x,idxEventType) = 1;
                    auroc_PC(x,idxEventType) = exp.expID(temp(1)).unit(idxUnit).event(idxEventType).auROC
                end
            end
        end
    end
end
for x = 1:100
    temp = sortedBest20(x,:);
    for idxUnit = 1:size(exp.expID(temp(1)).unit, 2)
        if temp(2) == exp.expID(temp(1)).unit(idxUnit).ID
            indices = idxUnit;
            for idxEventType = 1:length(exp.expID(temp(1)).unit(indices).event)
                if exp.expID(temp(1)).unit(idxUnit).event(idxEventType).excitatoryResponse == 1
                    v_PCx2(x,idxEventType) = 1;
                end
            end
        end
    end
end
% PRh

load('cfg_PRh_dm'); load('database_PRh_dm'); load('sortedBest20_PRh_svmtot_0-2');
p = 0; indices = []; v_PRh = zeros(25,14);
% A = squeeze(mean(mean(squeeze(mean(DECODING_RESULTS.FP_INFO{1, 2}.the_p_values_org_order(:,:,2:3,:),3)),1),2));
% [output,indx] = sort(A);
% sortedBestAll = [binned_site_info.sessionID(indx) binned_site_info.clusterID(indx)];

for x = 1:15
    temp = sortedBest20(x,:);
    for idxUnit = 1:size(exp.expID(temp(1)).unit, 2)
        if temp(2) == exp.expID(temp(1)).unit(idxUnit).ID
            indices = idxUnit;
            for idxEventType = 1:length(exp.expID(temp(1)).unit(indices).event)
                if exp.expID(temp(1)).unit(idxUnit).event(idxEventType).excitatoryResponse == 1
                    v_PRh(x,idxEventType) = 1;
                    auroc_PRh(x,idxEventType) = exp.expID(temp(1)).unit(idxUnit).event(idxEventType).auROC
                end
            end
        end
    end
end
for x = 1:100
    temp = sortedBest20(x,:);
    for idxUnit = 1:size(exp.expID(temp(1)).unit, 2)
        if temp(2) == exp.expID(temp(1)).unit(idxUnit).ID
            indices = idxUnit;
            for idxEventType = 1:length(exp.expID(temp(1)).unit(indices).event)
                if exp.expID(temp(1)).unit(idxUnit).event(idxEventType).excitatoryResponse == 1
                    v_PRh2(x,idxEventType) = 1;
                end
            end
        end
    end
end

figure
represOdor(:,1) = sum(v_PCx)/100;
represOdor(:,2) = sum(v_PRh)/100;
d = bar(represOdor,'FaceColor','flat','BarWidth', 1);
d(1).FaceColor = [0 0 1]; d(2).FaceColor = [1 0 0]; d(1).FaceAlpha = .6; d(2).FaceAlpha = .6; d(1).EdgeColor = 'none'; d(2).EdgeColor = 'none';
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
xlabel('Odor ID'); ylabel('% of Neurons'); title('Percentage of Significant Responses for each Odor');
hold on
represOdor2(:,1) = sum(v_PCx2)/100;
represOdor2(:,2) = sum(v_PRh2)/100;
d2 = bar(represOdor2,'FaceColor','flat','BarWidth', 1); ylim([0 1])
d2(1).FaceColor = 'none'; d2(2).FaceColor = 'none'; d2(1).EdgeColor = [0 0 1]; d2(2).EdgeColor = [1 0 0];
[ii,~,v] = find(auroc_PC');
out_PCx = accumarray(ii,v,[],@mean);
[ii,~,v] = find(auroc_PRh');
out_PRh = accumarray(ii,v,[],@mean);
plot(1:14,out_PCx, 'b')
plot(1:14,out_PRh, 'r')
hold off
legend('PCx 15N','PRh 15N','PCx 100N','PRh 100N', 'Location','east'); legend('boxoff');
%% RasterPlots

% PCx 

figure( 'Renderer', 'painters','Position', [10 10 1200 400])
tiledlayout(2,4,'TileSpacing','compact')
nexttile
colormap(brewermap([],'*YlGnBu'));
load('SetBasic_Binned_data_1000ms_PCx_exc_10ms_bins_10ms_sampled')
esempio = (binned_data{1, 60})
crop = esempio(11:20,1:400);
imagesc(crop, [0 0.1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
h = gca; h.XAxis.Visible = 'off';
ylabel('Unit \# 60'); title('Odor \#2')
hold on
M = 10.8*ones(1,201);
ax = plot(100:300, M, 'k', 'LineWidth', 3);
ax.Clipping = 'off';
hold off

nexttile
crop = esempio(131:140,1:400);
imagesc(crop, [0 0.1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
h = gca; h.XAxis.Visible = 'off'; title('Odor \#14')
hold on
M = 10.8*ones(1,201);
ax = plot(100:300, M, 'k', 'LineWidth', 3);
ax.Clipping = 'off';
hold off


nexttile([5])
colormap(brewermap([],'*YlGnBu'));
load('SetBasic_Binned_data_1000ms_PCx_inh_10ms_bins_10ms_sampled')
esempio = (binned_data{1, 2})
crop = esempio(11:20,1:400);
imagesc(crop, [0 0.1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
h = gca; h.XAxis.Visible = 'off';
ylabel('Unit \# 2');
hold on
M = 10.8*ones(1,201);
ax = plot(100:300, M, 'k', 'LineWidth', 3);
ax.Clipping = 'off';
hold off

nexttile(6)
crop = esempio(131:140,1:400);
imagesc(crop, [0 0.1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
h = gca; h.XAxis.Visible = 'off';
hold on
M = 10.8*ones(1,201);
ax = plot(100:300, M, 'k', 'LineWidth', 3);
ax.Clipping = 'off';
hold off


nexttile
colormap(brewermap([],'*YlGnBu'));
load('SetBasic_Binned_data_1000ms_PRh_exc_10ms_bins_10ms_sampled')
esempio = (binned_data{1, 93})
crop = esempio(11:20,1:400);
imagesc(crop, [0 0.1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
h = gca; h.XAxis.Visible = 'off';
ylabel('Unit \# 93'); title('Odor \#2')
hold on
M =10.8*ones(1,201);
ax = plot(100:300, M, 'k', 'LineWidth', 3);
ax.Clipping = 'off';
hold off

nexttile
crop = esempio(131:140,1:400);
imagesc(crop, [0 0.1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
h = gca; h.XAxis.Visible = 'off'; title('Odor \#14')
hold on
M = 10.8*ones(1,201);
ax = plot(100:300, M, 'k', 'LineWidth', 3);
ax.Clipping = 'off';
hold off


nexttile([7])
colormap(brewermap([],'*YlGnBu'));
load('SetBasic_Binned_data_1000ms_PRh_inh_10ms_bins_10ms_sampled')
esempio = (binned_data{1, 130})
crop = esempio(11:20,1:400);
imagesc(crop, [0 0.1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
h = gca; h.XAxis.Visible = 'off';
ylabel('Unit \# 130');
hold on
M = 10.8*ones(1,201);
ax = plot(100:300, M, 'k', 'LineWidth', 3);
ax.Clipping = 'off';
hold off

nexttile(8)
crop = esempio(131:140,1:400);
imagesc(crop, [0 0.1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
h = gca; h.XAxis.Visible = 'off';
hold on
M = 10.8*ones(1,201);
ax = plot(100:300, M, 'k', 'LineWidth', 3);
ax.Clipping = 'off';
hold off
%% Automated Raster
clear all

% PCx

figure
load('sortedBest20_PCx_svmtot_0-2');
load('Binned_data_200by100ms_BEST25_maxCorr_PCx_200ms_bins_100ms_sampled')
Labels = {'Od.1', 'Od.2','Od.3','Od.4','Od.5', 'Od.6', 'Od.7', 'Od.8', 'Od.9', 'Od.10', 'Od.11', 'Od.12', 'Od.13', 'Od.14'};
for ii = 1:15
    subplot(3,5,ii)
    ax = imagesc(squeeze(mean(reshape(binned_data{1,sortedBest20(ii,5)},[10,14,49]))), [0 0.04]); set(gca, 'XTick', linspace(0, 49, 11), 'XTickLabel', -1:0.5:4);
    title(['M# ' num2str(sortedBest20(ii,1)) ' ID# ' num2str(sortedBest20(ii,2))  ', LS ' num2str(round(sortedBest20(ii,3),1))]);
    set(gca, 'YTick', 1:14, 'YTickLabel', Labels);
    set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out'); set(gca,'xtick',[]);
    hold on
    M = 147*ones(1,20);
    ax = plot(10:29, M, 'k', 'LineWidth', 3);
    ax.Clipping = 'off';
    hold off
end
sgtitle('PCx')


% PRh

figure
load('sortedBest20_PRh_svmtot_0-2');
load('Binned_data_200by100ms_BEST25_maxCorr_PRh_200ms_bins_100ms_sampled')
Labels = {'Od.1', 'Od.2','Od.3','Od.4','Od.5', 'Od.6', 'Od.7', 'Od.8', 'Od.9', 'Od.10', 'Od.11', 'Od.12', 'Od.13', 'Od.14'};
for ii = 1:15
    subplot(3,5,ii)
    ax = imagesc(squeeze(mean(reshape(binned_data{1,sortedBest20(ii,5)},[10,14,49]))), [0 0.04]); set(gca, 'XTick', linspace(0, 49, 11), 'XTickLabel', -1:0.5:4);
    title(['M# ' num2str(sortedBest20(ii,1)) ' ID# ' num2str(sortedBest20(ii,2)) ', LS ' num2str(round(sortedBest20(ii,3),1))]);
    set(gca, 'YTick', 1:14, 'YTickLabel', Labels);
    set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out'); set(gca,'xtick',[]);
    hold on
    M = 147*ones(1,20);
    ax = plot(10:29, M, 'k', 'LineWidth', 3);
    ax.Clipping = 'off';
    hold off
end
sgtitle('PRh')

%% Average PSTH

% PCx

clear all
load('cfg_PCx_dm'); load('database_PCx_dm'); load('sortedBest20_PCx_svmtot_0-2.mat');
p = 0; indices = []; eventResponseAll = []; i = 0; KK_PRh = []; RPRh=[]; z=0; auROC_PCx = []; auROC_PCx_overall = [];
zz = 0; KK_anti = []; RPRh_anti= []; stdBSLPCx = [];

for x = 1:15
    temp = sortedBest20(x,:);
    for idxUnit = 1:size(exp.expID(temp(1)).unit, 2)
        if temp(2) == exp.expID(temp(1)).unit(idxUnit).ID
            indices = idxUnit;
            for idxEventType = 1:length(exp.expID(temp(1)).unit(indices).event)
                auROC_PCx_overall = [auROC_PCx_overall exp.expID(temp(1)).unit(indices).event(idxEventType).auROC];
                if exp.expID(temp(1)).unit(idxUnit).event(idxEventType).excitatoryResponse == 1
                    z = z + 1;
                    [R,t] = psth(exp.expID(temp(1)).unit(indices).event(idxEventType).spikesPerTrial, ...
                        0.025,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                    meanBSLPCx(z) = mean(R(21:1620)); stdBSLPCx(z) = std(R(21:1620));
                    RPRh = [RPRh; R];
                    K = (R - meanBSLPCx(z))./stdBSLPCx(z);
                    KK_PRh = [KK_PRh ; K]; [KK_peaks(z) KK_peaks_indx] = max(K(1620:end));
                    KK_ts(z) = t(KK_peaks_indx + 1619);
                    KK_tso_indx = find((K(1620:end)>3),1);
                    if isempty(KK_tso_indx)
                       KK_tso(z) = nan
                    else
                    KK_tso(z) = t(KK_tso_indx + 1619);
                    end
                    auROC_PCx = [auROC_PCx exp.expID(temp(1)).unit(indices).event(idxEventType).auROC];
                    sR_PCx(z,:) = exp.expID(temp(1)).unit(indices).event(idxEventType).spikeResponse;
                else
                    zz = zz + 1;
                    [R,t] = psth(exp.expID(temp(1)).unit(indices).event(idxEventType).spikesPerTrial, ...
                        0.025,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                    meanBSLPCx_anti(zz) = mean(R(21:1620)); stdBSLPCx_anti(zz) = std(R(21:1620));
                    RPRh_anti = [RPRh_anti; R];
                    K_anti = (R - meanBSLPCx_anti(zz))./stdBSLPCx_anti(zz);
                    KK_anti = [KK_anti ; K_anti];
                end 
            end
        end
    end
%             KK_peaks_bh(x) = mean(KK_peaks);
%             KK_ts_bh(x) = mean(KK_ts);
%             KK_peaks = []; KK_ts = [];   
end
peakInf_PCx = [KK_peaks' , KK_ts'];
figure
meanKK = mean(KK_PRh);
stdKK = std(KK_PRh);
errKK = stdKK./sqrt(size(KK_PRh,1));
k1 = shadedErrorBar(t, meanKK, errKK, '-b',1);
set(k1.edge,'LineWidth',.5,'LineStyle','none');
k1.mainLine.LineWidth = 2;
k1.patch.FaceAlpha = .3;
hold on
meanKK_anti = mean(KK_anti);
stdKK_anti = std(KK_anti);
errKK_anti = stdKK_anti./sqrt(size(KK_anti,1));
k2 = shadedErrorBar(t, meanKK_anti, errKK_anti, '--b',1);
set(k2.edge,'LineWidth',.5,'LineStyle','none');
k2.mainLine.LineWidth = 2;
k2.patch.FaceAlpha = .3;

% PRh

load('cfg_PRh_dm'); load('database_PRh_dm'); load('sortedBest20_PRh_svmtot_0-2.mat');
p = 0; indices = []; eventResponseAll = []; i = 0; ZZ = []; RPRh=[]; v_PRh=0; auROC_PRh = []; auROC_PRh_overall = [];
vv = 0; RPRh_anti = []; ZZ_anti = []; ZZ_peaks = []; ZZ_peaks_indx = [];

for x = 1:15
    temp = sortedBest20(x,:);
    for idxUnit = 1:size(exp.expID(temp(1)).unit, 2)
        if temp(2) == exp.expID(temp(1)).unit(idxUnit).ID
            indices = idxUnit;
            for idxEventType = 1:length(exp.expID(temp(1)).unit(indices).event)
                auROC_PRh_overall = [auROC_PRh_overall exp.expID(temp(1)).unit(indices).event(idxEventType).auROC];
                if exp.expID(temp(1)).unit(idxUnit).event(idxEventType).excitatoryResponse == 1
                    v_PRh = v_PRh + 1;
                    [R,t] = psth(exp.expID(temp(1)).unit(indices).event(idxEventType).spikesPerTrial, ...
                        0.025,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                    meanBSLPCx(v_PRh) = mean(R(21:1620)); stdBSLPRh(v_PRh) = std(R(21:1620));
                    RPRh = [RPRh; R];
                    Z = (R - meanBSLPCx(v_PRh))./stdBSLPRh(v_PRh);
                    ZZ = [ZZ ; Z]; [ZZ_peaks(v_PRh) ZZ_peaks_indx] = max(Z(1620:end));
                    ZZ_ts(v_PRh) = t(ZZ_peaks_indx + 1619);
                    ZZ_tso_indx = find((Z(1620:end)>3),1);
                    if isempty(ZZ_tso_indx)
                       ZZ_tso(v_PRh) = nan
                    else
                    ZZ_tso(v_PRh) = t(ZZ_tso_indx + 1619);
                    end
                    auROC_PRh = [auROC_PRh exp.expID(temp(1)).unit(indices).event(idxEventType).auROC];
                    sR_PRh(v_PRh,:) = exp.expID(temp(1)).unit(indices).event(idxEventType).spikeResponse;
                else
                    vv = vv + 1;
                    [R,t] = psth(exp.expID(temp(1)).unit(indices).event(idxEventType).spikesPerTrial, ...
                        0.025,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                    meanBSLPRh_anti(vv) = mean(R(21:1620)); stdBSLPRh_anti(vv) = std(R(21:1620));
                    RPRh_anti = [RPRh_anti; R];
                    Z_anti = (R - meanBSLPRh_anti(vv))./stdBSLPRh_anti(vv);
                    ZZ_anti = [ZZ_anti ; Z_anti];
                end
            end
        end
    end
%             ZZ_peaks_bh(x) = mean(ZZ_peaks);
%             ZZ_ts_bh(x) = mean(ZZ_ts);
%             ZZ_peaks = []; ZZ_ts = [];
end
peakInf_PRh = [ZZ_peaks' , ZZ_ts'];
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out');
meanZZ = mean(ZZ);
stdZZ = std(ZZ);
errZZ = stdZZ./sqrt(size(ZZ,1));
z = shadedErrorBar(t, meanZZ, errZZ, '-r',1);
set(z.edge,'LineWidth',.5,'LineStyle','none');
z.mainLine.LineWidth = 2;
z.patch.FaceAlpha = .3;
meanZZ_anti = mean(ZZ_anti);
stdZZ_anti = std(ZZ_anti);
errZZ_anti = stdZZ_anti./sqrt(size(ZZ_anti,1));
z2 = shadedErrorBar(t, meanZZ_anti, errZZ_anti, '--r',1);
set(z2.edge,'LineWidth',.5,'LineStyle','none');
z2.mainLine.LineWidth = 2;
z2.patch.FaceAlpha = .3;

M = -1*ones(1,3);
ax = plot(0:2, M, 'k', 'LineWidth', 3);
ax.Clipping = 'off';
hold on
legend([k1.mainLine, k2.mainLine, z.mainLine, z2.mainLine],'PCx sign', 'PCx ns', 'PRh sign', 'PRh ns','NumColumns',2,'Location', 'northwest');
legend('boxoff')
title('Average PSTH Best 25 Neurons');
xlabel('Time (s)');
ylabel('SD over baseline');
xlim([-1 4]); ylim([-2 12]);
axes('Position',[.7 .55 .20 .30])
plot(t,stdKK,'b', 'LineWidth', 1.5);
hold on
plot(t,stdZZ, 'r','LineWidth', 1.5);
M = zeros(1,3);
ax = plot(0:2, M, 'k', 'LineWidth', 1.5);
ax.Clipping = 'off'; ylim([-0.5 13]); xlim([-1 4]); 
hold off
title("Responses' SD Over Time"); ylabel('SDs'); set(gca, 'box', 'off'); 

figure
plot(t,stdKK./meanKK)
ylim([0 5]); xlim([0 4]);
%% Peaks, Latency and Noise Distributions

% Peak Latency

figure
subplot(2,1,1)
edges = [0:.2:4];
f1 = histogram(ZZ_ts, edges, 'FaceColor', 'r', 'EdgeColor','none', 'Normalization', 'probability');
hold on
f2 = histogram(KK_ts, edges, 'FaceColor', 'b', 'EdgeColor','none', 'Normalization', 'probability');
xlabel('Peak Onset (s)')
ylabel('Fraction of Cell-Odor Pairs')
ylim([0 0.4]); xlim([-0.05 4]);
legend('PRh','PCx'); legend('boxoff');
M = .4*ones(1,3);
ax = plot(0:2, M, 'k', 'LineWidth', 2);
ax.Clipping = 'off';
legend('PRh','PCx'); legend('boxoff');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
hold off
[p1 ,h1] = ranksum(ZZ_ts, KK_ts);


% Peak Latency

subplot(2,1,2)
edges = [0:.05:2];
j1 = histogram(ZZ_tso, edges, 'FaceColor', 'r', 'EdgeColor','none', 'Normalization', 'probability');
hold on
j2 = histogram(KK_tso, edges, 'FaceColor', 'b', 'EdgeColor','none', 'Normalization', 'probability');
xlabel('3SD over Baseline (s)')
ylabel('Fraction of Cell-Odor Pairs')
ylim([0 0.25]); xlim([-.01 1]);
legend('PRh','PCx'); legend('boxoff');
legend('PRh','PCx'); legend('boxoff');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
hold off
[p2 ,h2] = ranksum(ZZ_tso, KK_tso);
hold off
sgtitle('Timing Analysis PSTH');


figure
edges = [0:2.5:70];
k1 = histogram(ZZ_peaks, edges, 'FaceColor', 'r', 'EdgeColor','none', 'Normalization', 'probability');
hold on
k2 = histogram(KK_peaks, edges, 'FaceColor', 'b', 'EdgeColor','none', 'Normalization', 'probability');
xlabel('Peak Amplitude (SD over Baseline)');
ylabel('Fraction of Cell-Odor Pairs');
ylim([0 0.4]); xlim([0,70]);
legend('PRh','PCx'); legend('boxoff');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
hold off
title('Amplitude Analysis PSTH');
[p3 ,h3] = ranksum(ZZ_peaks, KK_peaks);


%% Reliability Over Trial

sR_PCx_std = std(sR_PCx,[],2);
sR_PCx_mean = mean(sR_PCx,2);
rel_PCx = sR_PCx_std./sR_PCx_mean; 
sR_PRh_std = std(sR_PRh,[],2);
sR_PRh_mean = mean(sR_PRh,2);
rel_PRh = sR_PRh_std./sR_PRh_mean;
histogram(rel_PCx, [0:0.2:3], 'FaceColor', 'b', 'EdgeColor','none','Normalization', 'probability')
hold on
histogram(rel_PRh, [0:0.2:3], 'FaceColor', 'r', 'EdgeColor','none','Normalization', 'probability')
xlabel('s.d./mean'); ylabel('Fraction of Cell-Odor Pairs');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
hold off
legend('PCx','PRh');
[p4 ,h4] = ranksum(rel_PCx, rel_PRh);
%% Test on Reliability

zsR_PCx = zscore(sR_PCx,[],2); 
for i = 1:size(zsR_PCx,1)
    plot(1:2, [mean(zsR_PCx(i,[1:5])) mean(zsR_PCx(i,[6:10]))]);
    hold on
end
PCx_first = mean(mean(zsR_PCx(:,[1:5]))); PCx_second = mean(mean(zsR_PCx(:,[6:10])));

zsR_PRh = zscore(sR_PRh,[],2); 
for i = 1:size(zsR_PRh,1)
    plot(1:2, [mean(zsR_PRh(i,[1:5])) mean(zsR_PRh(i,[6:10]))]);
    hold on
end
PRh_first = mean(mean(zsR_PRh(:,[1:5]))); PRh_second = mean(mean(zsR_PRh(:,[6:10])));

%% auROC Excitatory Responses

figure
subplot(1,2,1)
h1 = histogram(auROC_PCx_overall, [0:0.1:1], 'Normalization', 'probability')
h1.FaceColor = 'none'; h1.EdgeColor = 'b'; h1.EdgeAlpha = .8;
ylim([0 .3])
hold on
sign1 = histogram([auROC_PCx nan(1,length(auROC_PCx_overall) - length(auROC_PCx))], [0:0.1:1], 'Normalization', 'probability')
sign1.FaceColor = 'b'; sign1.EdgeColor = 'none';
title('PCx');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
subplot(1,2,2)
h2 =  histogram(auROC_PRh_overall, [0:0.1:1], 'Normalization', 'probability')
h2.FaceColor = 'none'; h2.EdgeColor = 'r'; h2.EdgeAlpha = .8;
ylim([0 .3])
hold on
sign2 = histogram([auROC_PRh nan(1,length(auROC_PRh_overall) - length(auROC_PRh))], [0:0.1:1], 'Normalization', 'probability')
sign2.FaceColor = 'r'; sign2.EdgeColor = 'none'; 
title('PRh');
sgtitle('Best 15 Neurons auROCs');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
hold off

%% Tuning Curves

% PCx

clear all
load('cfg_PCx_dm'); load('database_PCx_dm'); load('sortedBest20_PCx_svmtot_0-2');

for x = 1:100
    temp = sortedBest20(x,:);
    for idxUnit = 1:size(exp.expID(temp(1)).unit, 2)
        if temp(2) == exp.expID(temp(1)).unit(idxUnit).ID  
            indices = idxUnit;
            eventResponse = 0;
            for idxEventType = 1:length(exp.expID(temp(1)).unit(indices).event)
                tuningCurve(idxEventType) = mean(exp.expID(temp(1)).unit(indices).event(idxEventType).spikeResponse);
            end
            app_PCx = tuningCurve ./ max(tuningCurve);
            appl_notsorted_PCx(x,:) = app_PCx;
            app_PCx = sort(app_PCx,'descend');
            appl_PCx(x,:) = app_PCx;
            
            Zappl_notsorted_PCx(x,:) = zscore(tuningCurve);
            Zapp_PCx = sort(tuningCurve,'descend');
            Zappl_PCx(x,:) = zscore(tuningCurve);
            
            end
        end
end

% PRh
load('cfg_PRh_dm'); load('database_PRh_dm'); load('sortedBest20_PRh_svmtot_0-2');

for x = 1:100
    temp = sortedBest20(x,:);
    for idxUnit = 1:size(exp.expID(temp(1)).unit, 2)
        if temp(2) == exp.expID(temp(1)).unit(idxUnit).ID  
            indices = idxUnit;
            for idxEventType = 1:length(exp.expID(temp(1)).unit(indices).event)
                tuningCurve(idxEventType) = mean(exp.expID(temp(1)).unit(indices).event(idxEventType).spikeResponse);
            end
            app_PRh = tuningCurve ./ max(tuningCurve);
            appl_notsorted_PRh(x,:) = app_PRh;
            app_PRh = sort(app_PRh,'descend');
            appl_PRh(x,:) = app_PRh;
            
            Zappl_notsorted_PRh(x,:) = zscore(tuningCurve);
            Zapp_PRh = sort(tuningCurve,'descend');
            Zappl_PRh(x,:) = zscore(tuningCurve);
            
            end
        end
end

subplot(1,2,1)
n1 = plot(1:14, appl_PCx, 'Color',[0.1 0.1 0.9 0.3], 'LineWidth', 0.5);
hold on
plot(1:14, mean(appl_PCx), 'Color',[0.1 0.1 0.9 0.7], 'LineWidth', 2); xticks(1:14);
ylim([-1 1.1]); xlim([0 15]); title('PCx'); 
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
hold off
xlabel('Odors (Preference Ordered)'); ylabel('Scaled Response')

subplot(1,2,2)
n2 = plot(1:14, appl_PRh, 'Color',[0.9 0.1 0.1 0.3], 'LineWidth', 0.5);
hold on
plot(1:14, mean(appl_PRh), 'Color',[0.9 0.1 0.1 0.7], 'LineWidth', 2); xticks(1:14);
ylim([-1 1.1]); xlim([0 15]); title('PRh');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
hold off
xlabel('Odors (Preference Ordered)'); ylabel('Scaled Response')
sgtitle('Tuning Curves - Best 15 Neurons');

figure
n1 = shadedErrorBar(1:14, mean(appl_PCx), std(appl_PCx)./sqrt(size(appl_PCx,1)), '-b');
set(n1.edge,'LineWidth',.5,'LineStyle','none');
n1.mainLine.LineWidth = 2;
n1.patch.FaceAlpha = .6;
xticks(1:14); ylim([-1 1.1]); xlim([0 15]); title('Average Tuning Curve');  xlabel('Odors (Preference Ordered)'); ylabel('Scaled Response')
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
hold on
n2 = shadedErrorBar(1:14, mean(appl_PRh), std(appl_PRh)./sqrt(size(appl_PRh,1)), '-r');
set(n2.edge,'LineWidth',.5,'LineStyle','none');
n2.mainLine.LineWidth = 2;
n2.patch.FaceAlpha = .6;
legend([n1.mainLine, n2.mainLine], 'PCx','PRh');
hold off
%% Correlations between neurons
subplot(1,2,1)
imagesc(abs(corr(Zappl_notsorted_PCx')), [0 1]); axis square; xticks(1:15),yticks(1:15);set(gca,'xaxisLocation','top'); title('PCx');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
subplot(1,2,2)
imagesc(abs(corr(Zappl_notsorted_PRh')), [0 1]); axis square; xticks(1:15),yticks(1:15);set(gca,'xaxisLocation','top'); title('PRh');
sgtitle('Pairwise Correlations Tuning Curve');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
%%
s = mean(corr(Zappl_notsorted_PCx'))
s2 = mean(corr(Zappl_notsorted_PRh'))
ranksum(s,s2)

%%

load('cfg_PCx_dm.mat');
load('database_PCx_dm.mat');
nResponsiveEventsFS_PCx = [];
nResponsiveEventsRS_PCx = [];
normTuningCurveFS_PCx = [];
normTuningCurveRS_PCx = [];
zs_TuningCurve_PCx = [];

for idxExp = 1:size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        eventResponse = 0;
        tuningCurve = nan(1, size(exp.expID(idxExp).unit(idxUnit).event, 2));
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2)
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
            tuningCurve(idxEvent) = mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse);    
        end
        if eventResponse > 0
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'FS')
                nResponsiveEventsFS_PCx = [nResponsiveEventsFS_PCx; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveFS_PCx = [normTuningCurveFS_PCx; app];
            else
                nResponsiveEventsRS_PCx = [nResponsiveEventsRS_PCx; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveRS_PCx = [normTuningCurveRS_PCx; app];
                zs_TuningCurve_PCx = [zs_TuningCurve_PCx; zscore(tuningCurve)];
            end
        end
    end
end

load('database_PRh_dm.mat');
load('cfg_PRh_dm.mat');

nResponsiveEventsFS_PRh = [];
nResponsiveEventsRS_PRh = [];
normTuningCurveFS_PRh = [];
normTuningCurveRS_PRh = [];
zs_TuningCurve_PRh = [];

for idxExp = 1:size(exp.expID,2) %size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        eventResponse = 0;
        tuningCurve = nan(1, size(exp.expID(idxExp).unit(idxUnit).event, 2));
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2)
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
            tuningCurve(idxEvent) = mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse);    
        end
        if eventResponse > 0
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'FS')
                nResponsiveEventsFS_PRh = [nResponsiveEventsFS_PRh; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveFS_PRh = [normTuningCurveFS_PRh; app];
            else
                nResponsiveEventsRS_PRh = [nResponsiveEventsRS_PRh; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveRS_PRh = [normTuningCurveRS_PRh; app];
                zs_TuningCurve_PRh = [zs_TuningCurve_PRh; zscore(tuningCurve)];
            end
        end
    end
end

U_PRh = corr(zs_TuningCurve_PRh');
U_PCx = corr(zs_TuningCurve_PCx');
At = U_PRh.';
m  = (1:size(At,1)).' >= (1:size(At,2));
v_PRh  = At(m);
v_PRh = v_PRh(v_PRh~=1);
At = U_PCx.';
m  = (1:size(At,1)).' >= (1:size(At,2));
v_PCx  = At(m);
v_PCx = v_PCx(v_PCx~=1);
 
%% Noise Correlations

load('database_PCx_dm.mat');
load('cfg_PCx_dm.mat');
tempMean = [];
all_NoiseVectors_PCx = [];
eventResponse = 0;
for idxExp = 1:size(exp.expID,2) %size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2)
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
        end
        if eventResponse > 0
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
                for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2)
                    tempMean = zscore([tempMean exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse - mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse)]);
                end
                all_NoiseVectors_PCx = [all_NoiseVectors_PCx ; tempMean];
                tempMean = [];
            end
            eventResponse = 0;
        end
    end
end

load('database_PRh_dm.mat');
load('cfg_PRh_dm.mat');
tempMean = [];
all_NoiseVectors_PRh = [];
eventResponse = 0;
for idxExp = 1:size(exp.expID,2) %size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2)
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
        end
        if eventResponse > 0
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
                for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2)
                    tempMean = zscore([tempMean exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse - mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse)]);
                end
                all_NoiseVectors_PRh = [all_NoiseVectors_PRh ; tempMean];
                tempMean = [];
            end
            eventResponse = 0;
        end
    end
end

U_PRh = corr(all_NoiseVectors_PRh');
U_PCx = corr(all_NoiseVectors_PCx');
At = U_PRh.';
m  = (1:size(At,1)).' >= (1:size(At,2));
vN_PRh  = At(m);
vN_PRh = vN_PRh(vN_PRh~=1);
At = U_PCx.';
m  = (1:size(At,1)).' >= (1:size(At,2));
vN_PCx  = At(m);
vN_PCx = vN_PCx(vN_PCx~=1);

%%
% example 1
f7 = figure;
subplot(1,2,1)
 [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
d{1} = v_PCx; d{2} = v_PRh;
h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cb(5,:), 'alpha', 0.7,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
     'box_col_match', 1);
h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.7,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';
h1{1, 2}.SizeData = 1; h2{1, 2}.SizeData = 1; h1{1, 2}.MarkerFaceAlpha = .7;; h2{1, 2}.MarkerFaceAlpha = .7;
title(['Signal Correlations']);
set(gca,'XLim', [-1 1]); set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
ax = gca;   %or as appropriate
yticks = get(ax, 'YTick'); yticks=yticks(yticks>=0); set(ax, 'YTick', yticks);
hold on 
%h3 = line([0 0], [0 1.5], 'LineStyle', '--');
%h3.Color = [0 0 0 .5];

subplot(1,2,2)
d{3} = vN_PCx; d{4} = vN_PRh;
h1 = raincloud_plot(d{3}, 'box_on', 1, 'color', cb(5,:), 'alpha', 0.7,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
     'box_col_match', 1);
h2 = raincloud_plot(d{4}, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.7,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';
h1{1, 2}.SizeData = 1; h2{1, 2}.SizeData = 1; h1{1, 2}.MarkerFaceAlpha = .7;; h2{1, 2}.MarkerFaceAlpha = .7;
title(['Noise Correlations']);
set(gca,'XLim', [-1 1]);set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
ax = gca;   %or as appropriate
yticks = get(ax, 'YTick'); yticks=yticks(yticks>=0); set(ax, 'YTick', yticks);
%yticklabels = get(ax, 'YTickLabel');
%yticklabels{1} = ''; yticklabels{2} = ''; yticklabels{3} = ''; set(ax, 'YTickLabel', yticklabels);
legend([h1{1} h2{1}], {'PCx', 'PRh'});
%h3 = line([0 0], [0 4], 'LineStyle', '--');
%h3.Color = [0 0 0 .5];
box off
[h,p]= ttest2(d{1},d{2});
mean(d{1}),mean(d{2});

%%
% Create sample data:
correlations = corr(Zappl_notsorted_PCx');
minValue = -1;
maxValue = 1;
% Scale the data to between -1 and +1.
correlations = (correlations-minValue) * 2 / (maxValue - minValue) - 1;
% Display - will use some weird color map to start with.
imagesc(correlations, [-1 1]);
colorbar
% Create colormap that is green for negative, red for positive,
% and a chunk inthe middle that is black.
greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
% Apply the colormap.
colormap(colorMap);
%% Tuning not Sorted

figure
subplot(1,2,1)
n1 = plot(1:14, appl_notsorted_PCx, 'Color',[0.1 0.1 0.9 0.3], 'LineWidth', 0.5);
hold on
plot(1:14, mean(appl_notsorted_PCx), 'Color',[0.1 0.1 0.9 0.8], 'LineWidth', 2); xticks(1:14);
ylim([-1 1.1]); xlim([0 15]); title('PCx'); 
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
hold off
subplot(1,2,2)
n2 = plot(1:14, appl_notsorted_PRh, 'Color',[0.9 0.1 0.1 0.3], 'LineWidth', 0.5);
hold on
plot(1:14, mean(appl_notsorted_PRh), 'Color',[0.9 0.1 0.1 0.8], 'LineWidth', 2); xticks(1:14);
ylim([-1 1.1]); xlim([0 15]); title('PRh');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
hold off
sgtitle('Best 15 Neurons Tuning Curves');


%% Confusion Matrices

figure
subplot(1,2,1)
load('PCx_-1_4_results_200by100ms_BEST25_maxCorr.mat');
e1 = imagesc(mean(DECODING_RESULTS.NORMALIZED_RANK_RESULTS.confusion_matrix_results.rank_confusion_matrix(:,:,10:19),3));
axis square;
subplot(1,2,2)
load('PRh_-1_4_results_200by100ms_BEST25_maxCorr.mat');
e2 = imagesc(mean(DECODING_RESULTS.NORMALIZED_RANK_RESULTS.confusion_matrix_results.rank_confusion_matrix(:,:,10:19),3));
axis square;

%%
figure; load('database_PCx_Principled_1Step_TOPn_def_15_neurons.mat');
set(gcf,'color','white', 'PaperPositionMode', 'auto');
periods{1} = '[0 1]';
periods{2} = '[1 2]';
for idxPeriod = 1:2
    CM = squeeze(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix(:,:,idxPeriod+1));
    subplot(2,2,idxPeriod)
    clims = [0 1];
    imagesc(CM ./ sum(CM(:,1)), clims)
    axis square
    xlabel('predicted stimulus')
    ylabel('actual stimulus');xticks(1:14); yticks(1:14);
    colorbar
    title(['PCx:' periods{idxPeriod}])
    set(gca, 'box', 'off', 'tickDir', 'out')
end
load('database_PRh_Principled_1Step_TOPn_def_15_neurons.mat');
set(gcf,'color','white', 'PaperPositionMode', 'auto');
for idxPeriod = 1:2
    CM = squeeze(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix(:,:,idxPeriod+1));
    subplot(2,2,idxPeriod+2)
    clims = [0 1];
    imagesc(CM ./ sum(CM(:,1)), clims)
    axis square
    xlabel('predicted stimulus')
    ylabel('actual stimulus'); xticks(1:14); yticks(1:14);
    colorbar
    title(['PRh:' periods{idxPeriod}])
    set(gca, 'box', 'off', 'tickDir', 'out')
end
%%
load('database_PCx_dm.mat');
load('cfg_PCx_dm.mat');
eventResponse = 0;
trl = 0;
sResp_PCx2 = [];
y= 0;
for idxExp = 1:size(exp.expID,2) %size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2)
                if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                    trl = trl + 1;
                    sResp_PCx(trl,:)= zscore(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse);
                    %temp2 = mean(reshape(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse, 2,5,1));
                    temp2 = exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse;
                    if temp2(1)>0
                       temp2 = (temp2./temp2(1));
                       sResp_PCx2 = [sResp_PCx2; temp2];
                    else
                       y = y+1;
                    end
                end
            end
        end  
    end
end

load('database_PRh_dm.mat');
load('cfg_PRh_dm.mat');
eventResponse = 0;
trl = 0;
sResp_PRh2 = [];
x= 0;
for idxExp = 1:size(exp.expID,2) %size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2)
                if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                    trl = trl + 1;
                    sResp_PRh(trl,:)= zscore(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse);
                    %temp2 = mean(reshape(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse, 2,5,1));
                    temp2 = exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse;
                    if temp2(1)>0
                       temp2 = (temp2./temp2(1));
                       sResp_PRh2 = [sResp_PRh2; temp2];
                    else
                        x = x +1;
                    end
                end
            end
        end  
    end
end

for i=1:size(sResp_PCx2,1)
    [g] = polyfit(1:10,sResp_PCx2(i,:),1);
    gcoeff_PCx(i) = g(1);
end
for i=1:size(sResp_PRh2,1)
    [g] = polyfit(1:10,sResp_PRh2(i,:),1);
    gcoeff_PRh(i) = g(1);
end

[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
d{1} =  gcoeff_PCx; d{2} = gcoeff_PRh;
h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cb(5,:), 'alpha', 0.7,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
     'box_col_match', 1);
h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.7,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';
h1{1, 2}.SizeData = 3; h2{1, 2}.SizeData = 3; h1{1, 2}.MarkerFaceAlpha = .7;; h2{1, 2}.MarkerFaceAlpha = .7;
title(['Adaptation Gradients']);
set(gca,'XLim', [-1 1],'YLim', [-9 9]); set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
ax = gca;   %or as appropriate
yticks = get(ax, 'YTick'); yticks=yticks(yticks>=0); set(ax, 'YTick', yticks);
hold off
%%
sResp_PCx2boot = [];sResp_PRh2boot = []; gboots_PRh = []; gboots_PCx = [];

for l = 1:200
    for i=1:size(sResp_PRh2,1)
        sResp_PRh2boot(i,:) = sResp_PRh2(i,randperm(10));
        [grand] = polyfit(1:10,sResp_PRh2boot(i,:),1);
        gboots_PRh = [gboots_PRh grand(1)];
    end
    for i=1:size(sResp_PCx2,1)
        sResp_PCx2boot(i,:) = sResp_PCx2(i,randperm(10));
        [grand] = polyfit(1:10,sResp_PCx2boot(i,:),1);
        gboots_PCx = [gboots_PCx grand(1)];
    end
end

figure
d{1} =  gboots_PCx; d{2} = gboots_PRh;
h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cb(5,:), 'alpha', 0.7,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
     'box_col_match', 1);
h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.7,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';
h1{1, 2}.SizeData = 1; h2{1, 2}.SizeData = 1; h1{1, 2}.MarkerFaceAlpha = .7;; h2{1, 2}.MarkerFaceAlpha = .7;
title(['Signal Correlations']);
set(gca,'XLim', [-1 1],'YLim', [-10 10]); set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
ax = gca;   %or as appropriate
yticks = get(ax, 'YTick'); yticks=yticks(yticks>=0); set(ax, 'YTick', yticks);