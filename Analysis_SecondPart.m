set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize',12)
set(0,'defaultLegendFontSize',20)
[cb_pca]=brewermap(8,'Set2');
cbFan_PCx = brewermap([],'*Blues');
cbFan_PRh = brewermap([],'*Reds');
cb = brewermap(9,'*Set3');
cb_d1 = cb(5,:);
cb_d5 = cb(4,:);
cb_d6 = cb(9,:);
colors_days = [cb_d1;cb_d5;cb_d6];
    
%% Behavior

load('cfg_A1_THESIS.mat');
load('database_A1_THESIS');

all_trials_A1 = [];
for idxExp = 1:size(exp.expID,2)
    for idxEventType = 1 : size(exp.expID(1).breathing,2)
        all_trials_A1 = [all_trials_A1 exp.expID(idxExp).breathing(idxEventType).reponse_peak ./ exp.expID(idxExp).breathing(idxEventType).baseline_mean];
    end 
end

load('cfg_A5_THESIS.mat');
load('database_A5_THESIS');

all_trials_A5 = [];
for idxExp = 1:size(exp.expID,2)
    for idxEventType = 1 : size(exp.expID(1).breathing,2)
        all_trials_A5 = [all_trials_A5 exp.expID(idxExp).breathing(idxEventType).reponse_peak ./ exp.expID(idxExp).breathing(idxEventType).baseline_mean];
    end
end

load('cfg_B1_THESIS.mat');
load('database_B1_THESIS');
all_trials = [];
all_trials_B1 = [];
for idxExp = 1:size(exp.expID,2)
    for idxEventType = 1 : size(exp.expID(1).breathing,2)
        all_trials_B1 = [all_trials_B1 exp.expID(idxExp).breathing(idxEventType).reponse_peak ./ exp.expID(idxExp).breathing(idxEventType).baseline_mean];
    end
end

figure
figure( 'Renderer', 'painters','Position', [10 10 900 500])
sniffs_A1 = mean(all_trials_A1,2);
sniffs_A1_sd = std(all_trials_A1,[],2)./sqrt(size(all_trials_A1,2));
a1sn = shadedErrorBar(1:15, sniffs_A1, sniffs_A1_sd ,'-k',1);
set(a1sn.edge,'LineWidth',.5,'LineStyle','none');
a1sn .mainLine.LineWidth = 1.5; a1sn .mainLine.Color= cb_d1;
a1sn .patch.FaceAlpha = .1; a1sn .patch.FaceColor= cb_d1;

hold on; 
sniffs_A5 = mean(all_trials_A5,2);
sniffs_A5_sd = std(all_trials_A5,[],2)./sqrt(size(all_trials_A5,2));
a5sn = shadedErrorBar(1:15, sniffs_A5, sniffs_A5_sd ,'-k',1);
set(a5sn.edge,'LineWidth',.5,'LineStyle','none');
a5sn .mainLine.LineWidth = 1.5; a5sn .mainLine.Color= cb_d5;
a5sn .patch.FaceAlpha = .1; a5sn .patch.FaceColor= cb_d5;
hold on
sniffs_B1 = mean(all_trials_B1,2);
sniffs_B1_sd = std(all_trials_B1,[],2)./sqrt(size(all_trials_B1,2));
b1sn = shadedErrorBar(1:15, sniffs_B1, sniffs_B1_sd ,'-k',1);
set(b1sn.edge,'LineWidth',.5,'LineStyle','none');
b1sn .mainLine.LineWidth = 1.5; b1sn .mainLine.Color= cb_d6;
b1sn .patch.FaceAlpha = .1; b1sn .patch.FaceColor= cb_d6;

title('Behavioural Read-Out of Recognition Memory');
ylim([1 4]); xlim([0 16]); xticks(1:2:15);
xlabel('Trial Number');
ylabel('Fold-Increase relative to Baseline');
legend([a1sn.mainLine, a5sn.mainLine, b1sn.mainLine],{'First Exposure - Odor Set A (A1)','Fifth Exposure - Odor Set A (A5)','First Exposure - Odor Set B (B1)'},'Location','northeast')
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
legend boxoff

for ii=1:42
gradA1(ii,:)=polyfit(1:15,all_trials_A1(:,ii)',1)
gradA5(ii,:)=polyfit(1:15,all_trials_A5(:,ii)',1)
gradB1(ii,:)=polyfit(1:15,all_trials_B1(:,ii)',1)
end

for ii=1:15
[p,h]=ranksum(all_trials_A5(ii,:),all_trials_B1(ii,:))
s(ii) = h 
end

mean(gradA1(:,2))
mean(gradA5(:,2))
mean(gradB1(:,2))
[h,p]=ranksum(gradA1(:,2),gradB1(:,2))
[h,p]=ranksum(gradA1(:,2),gradA5(:,2))
[h,p]=ranksum(gradA5(:,2),gradB1(:,2))
[h,p]=ttest(gradA1(:,2),gradB1(:,2))
[h,p]=ranksum(gradA1(:,2),gradA5(:,2))
[h,p]=ranksum(gradA5(:,2),gradB1(:,2))
%%

figure('Renderer', 'painters', 'Position', [10 10 900 300])
tiledlayout(1,2,'TileSpacing','compact')
nexttile
load('cfg_A1_THESIS.mat');
load('database_A1_THESIS');

listNeuron_D1= [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        excCounter = 0;
        inhCounter = 0;
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                inhCounter = inhCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse;
                excCounter = excCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse;
        end
        listNeuron_D1 = [listNeuron_D1 ; excCounter inhCounter];
    end
end


xRange = 0:1:7;
[N, edges] = histcounts(listNeuron_D1(:,1), -0.5:1:7.5);
N = N ./ numel(listNeuron_D1(:,1));
t1 = plot(xRange, N, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', cb_d1, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 7)
hold on
errExc_D1 = sqrt((N.*(1-N))/numel(listNeuron_D1(:,1)));
errorbar(xRange,  N, errExc_D1, 'color', cb_d1,'LineStyle', 'none')
xlim([-1 8]);
xticks([0:1:7]);
title('Excitatory');
xlabel('Number of Odors');
ylabel('Fraction of Neurons');
ylim([-0.05 1]);
hold on

load('cfg_A5_THESIS.mat');
load('database_A5_THESIS');
listNeuron_D5= [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        excCounter = 0;
        inhCounter = 0;
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                inhCounter = inhCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse;
                excCounter = excCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse;
        end
        listNeuron_D5 = [listNeuron_D5 ; excCounter inhCounter];
    end
end

[N, edges] = histcounts(listNeuron_D5(:,1), -0.5:1:7.5);
N = N ./ numel(listNeuron_D5(:,1));
t2 = plot(xRange, N, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', cb_d5, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 7)
hold on
errExc_D5 = sqrt((N.*(1-N))/numel(listNeuron_D5(:,1)));
errorbar(xRange,  N, errExc_D5, 'color', cb_d5,'LineStyle', 'none')
xlim([-1 8]);
xticks([0:1:7]);
title('Excitatory');
xlabel('Number of Odors');
ylabel('Fraction of Neurons');
ylim([-0.05 1]);

load('cfg_B1_THESIS.mat');
load('database_B1_THESIS');
listNeuron_D6= [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        excCounter = 0;
        inhCounter = 0;
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                inhCounter = inhCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse;
                excCounter = excCounter + exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse;
        end
        listNeuron_D6 = [listNeuron_D6 ; excCounter inhCounter];
    end
end


xRange = 0:1:7;
[N, edges] = histcounts(listNeuron_D6(:,1), -0.5:1:7.5);
N = N ./ numel(listNeuron_D6(:,1));
t3 = plot(xRange, N, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', cb_d6, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 7)
hold on
errExc_D6 = sqrt((N.*(1-N))/numel(listNeuron_D6(:,1)));
errorbar(xRange,  N, errExc_D6, 'color', cb_d6,'LineStyle', 'none')
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');



nexttile

xRange = 0:1:7;
[N, edges] = histcounts(listNeuron_D1(:,2), -0.5:1:7.5);
N = N ./ numel(listNeuron_D1(:,2));
t1 = plot(xRange, N, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', cb_d1, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 7)
hold on
errExc_D1 = sqrt((N.*(1-N))/numel(listNeuron_D1(:,2)));
errorbar(xRange,  N, errExc_D1, 'color', cb_d1,'LineStyle', 'none')
hold on

[N, edges] = histcounts(listNeuron_D5(:,2), -0.5:1:7.5);
N = N ./ numel(listNeuron_D5(:,2));
t2 = plot(xRange, N, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', cb_d5, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 7)
hold on
errExc_D5 = sqrt((N.*(1-N))/numel(listNeuron_D5(:,2)));
errorbar(xRange,  N, errExc_D5, 'color', cb_d5,'LineStyle', 'none')
xlim([-1 8]);
xticks([0:1:7]);
title('Inhibitory');
xlabel('Number of Odors');
ylabel('Fraction of Neurons');
ylim([-0.05 1]);

[N, edges] = histcounts(listNeuron_D6(:,2), -0.5:1:7.5);
N = N ./ numel(listNeuron_D6(:,2));
t3 = plot(xRange, N, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', cb_d6, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 7)
hold on
errExc_D6 = sqrt((N.*(1-N))/numel(listNeuron_D6(:,2)));
errorbar(xRange,  N, errExc_D6, 'color', cb_d6,'LineStyle', 'none')

set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');set(gca, 'YTick', []);
ax1 = gca;            
ax1.YAxis.Visible = 'off';
legend([t1, t2,t3],{'Day 1 - Odor Set A', 'Day 5 - Odor Set A','Day 1 - Odor Set B'});
legend boxoff

[h,p]=kstest2(listNeuron_D1(:,1),listNeuron_D5(:,1))
[h,p]=kstest2(listNeuron_D6(:,1),listNeuron_D5(:,1))
[h,p]=kstest2(listNeuron_D1(:,1),listNeuron_D6(:,1))
%%

%% auROCs distributions



load('cfg_A1_THESIS.mat');
load('database_A1_THESIS');
auROCs_A1 = [];
inH = [];
exC = [];

for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
       % if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                auROCs_A1 = [auROCs_A1 exp.expID(idxExp).unit(idxUnit).event(idxEvent).auROC];
                inH = [ inH exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse];
                exC = [ exC exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse];
            end
       % end
    end
end
figure('Renderer', 'painters', 'Position', [10 10 900 500])
%subplot(1,2,1)
tiledlayout(5,3, 'TileSpacing','compact')
nexttile([3 1])
pl = 0:0.05:1;
h1 = histogram(auROCs_A1, pl, 'Normalization', 'probability');
h1.FaceColor = 'none';
h1.EdgeColor = cb_d1;
h1.EdgeAlpha = 1;
h1.LineWidth = 1;

excNaN = length(auROCs_A1) - length(auROCs_A1(exC==1));
hold on
exc1 = histogram([auROCs_A1(exC==1) nan(1,excNaN)], pl, 'Normalization', 'probability');
exc1.FaceColor = cb_d1;
exc1.FaceAlpha = .85;
exc1.EdgeColor = [0 0 0];
exc1.EdgeAlpha = 0;

inhNaN = length(auROCs_A1) - length(auROCs_A1(inH==1));
hold on
inh = histogram([auROCs_A1(inH==1) nan(1,inhNaN)], pl, 'Normalization', 'probability');
inh.FaceColor = cb_d1;
inh.FaceAlpha = .85;
inh.EdgeColor = [0 0 0];
inh.EdgeAlpha = 0;
title('Day 1 - Odor Set A');
xlabel('auROC'); ylabel('Fraction of Cell-Odor Pairs');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
% legend(exc1, 'Significant Responses')
% legend boxoff
ylim([0 0.25]);
exca1 = length(find(exC==1))
inha1 = length(find(inH==1))

load('cfg_A5_THESIS.mat');
load('database_A5_THESIS');
auROCs_D5 = [];
inH = [];
exC = [];

for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
      %  if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                auROCs_D5 = [auROCs_D5 exp.expID(idxExp).unit(idxUnit).event(idxEvent).auROC];
                inH = [ inH exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse];
                exC = [ exC exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse];
            end
      %  end
    end
end
%subplot(1,2,2)
nexttile([3 1])
pl = 0:0.05:1;
h2 = histogram(auROCs_D5, pl, 'Normalization', 'probability');
h2.FaceColor = 'none';
h2.EdgeColor = cb_d5;
h2.EdgeAlpha = 1;
h2.LineWidth = 1;

excNaN = length(auROCs_D5) - length(auROCs_D5(exC==1));
hold on
exc2 = histogram([auROCs_D5(exC==1) nan(1,excNaN)], pl, 'Normalization', 'probability');
exc2.FaceColor = cb_d5;
exc2.FaceAlpha = .85;
exc2.EdgeColor = [0 0 0];
exc2.EdgeAlpha = 0;

inhNaN = length(auROCs_D5) - length(auROCs_D5(inH==1));
hold on
inh = histogram([auROCs_D5(inH==1) nan(1,inhNaN)], pl, 'Normalization', 'probability');
inh.FaceColor = cb_d5;
inh.FaceAlpha = .85;
inh.EdgeColor = [0 0 0];
inh.EdgeAlpha = 0;
title('Day 5 - Odor Set A');
xlabel('auROC');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');set(gca, 'YTick', []);
ax1 = gca;            
ax1.YAxis.Visible = 'off';
% legend(exc2, 'Significant Responses')
% legn = legend([h1, h2],'PCx','PRh')
% legn.Position = [0.8 0.8 0.1 0.1]
ylim([0 0.25]);
exca5 = length(find(exC==1))
inha5 = length(find(inH==1))

load('cfg_B1_THESIS.mat');
load('database_B1_THESIS');
auROCs_D6 = [];
inH = [];
exC = [];

for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
       % if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                auROCs_D6 = [auROCs_D6 exp.expID(idxExp).unit(idxUnit).event(idxEvent).auROC];
                inH = [ inH exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse];
                exC = [ exC exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse];
            end
      %  end
    end
end
%subplot(1,2,2)
nexttile([3 1])
pl = 0:0.05:1;
h2 = histogram(auROCs_D6, pl, 'Normalization', 'probability');
h2.FaceColor = 'none';
h2.EdgeColor = cb_d6;
h2.EdgeAlpha = 1;
h2.LineWidth = 1;

excNaN = length(auROCs_D6) - length(auROCs_D6(exC==1));
hold on
exc2 = histogram([auROCs_D6(exC==1) nan(1,excNaN)], pl, 'Normalization', 'probability');
exc2.FaceColor = cb_d6;
exc2.FaceAlpha = .85;
exc2.EdgeColor = [0 0 0];
exc2.EdgeAlpha = 0;

inhNaN = length(auROCs_D6) - length(auROCs_D6(inH==1));
hold on
inh = histogram([auROCs_D6(inH==1) nan(1,inhNaN)], pl, 'Normalization', 'probability');
inh.FaceColor = cb_d6;
inh.FaceAlpha = .85;
inh.EdgeColor = [0 0 0];
inh.EdgeAlpha = 0;
title('Day 1 - Odor Set B');
xlabel('auROC');
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');set(gca, 'YTick', []);
ax1 = gca;            
ax1.YAxis.Visible = 'off';
% legend(exc2, 'Significant Responses')
% legn = legend([h1, h2],'PCx','PRh')
% legn.Position = [0.8 0.8 0.1 0.1]
ylim([0 0.25]);
excb1 = length(find(exC==1))
inhb1 = length(find(inH==1))

 n1 = exca1; N1 = length(auROCs_A1);
 n2 = excb1; N2 = length(auROCs_D6);
 x1 = [repmat('a',N1,1); repmat('b',N2,1)];
 x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
 [tbl,chi2stat,pval] = crosstab(x1,x2)



nexttile([2 3])
sbar=bar([exca1/length(auROCs_A1),exca5/length(auROCs_D5),excb1/length(auROCs_D6); ...
    inha1/length(auROCs_A1),inha5/length(auROCs_D5),inhb1/length(auROCs_D6)])

error_a1 = sqrt((exca1/length(auROCs_A1))*(1-(exca1/length(auROCs_A1)))/length(auROCs_A1))
error_a5 = sqrt((exca5/length(auROCs_D5))*(1-(exca5/length(auROCs_D5)))/length(auROCs_D5))
error_b1 = sqrt((excb1/length(auROCs_D6))*(1-(excb1/length(auROCs_D6)))/length(auROCs_D6))
error_a1I = sqrt((inha1/length(auROCs_A1))*(1-(inha1/length(auROCs_A1)))/length(auROCs_A1))
error_a5I = sqrt((inha5/length(auROCs_D5))*(1-(inha5/length(auROCs_D5)))/length(auROCs_D5))
error_b1I = sqrt((inhb1/length(auROCs_D6))*(1-(inhb1/length(auROCs_D6)))/length(auROCs_D6))
ylim([0 0.20])
%Get the x coordinate of the bars
hold on
x = [];
for i = 1:3
    x = [x ; sbar(i).XEndPoints];
end
%Plot the errorbars
errb = errorbar(x',[exca1/length(auROCs_A1) exca5/length(auROCs_D5) excb1/length(auROCs_D6); inha1/length(auROCs_A1) inha5/length(auROCs_D5) inhb1/length(auROCs_D6)],[error_a1 error_a5 error_b1 ; error_a1I error_a5I error_b1I],'k','linestyle','none')'


sbar(1, 1).FaceColor =  cb_d1;sbar(1,2).FaceColor =  cb_d5;sbar(1,3).FaceColor =  cb_d6;
sbar(1, 1).LineStyle =  'none'; sbar(1, 2).LineStyle =  'none';sbar(1, 3).LineStyle =  'none';
errb(3,1).Color = cb_d6
errb(2,1).Color = cb_d5
errb(1,1).Color = cb_d1
% hold off
lg  = legend('A1', 'A5','B1', 'Location', 'northwest'); 
legend boxoff
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
xticklabels({'Excitatory','Inhibitory'})
title('Fraction of Significant Cell-Odor Pairs')
ylabel('Fraction of Cell-Odor Pairs')
hold on

[h,p]=ranksum(auROCs_D5,auROCs_A1)
[h,p]=ranksum(auROCs_D5,auROCs_D6)
[h,p]=ranksum(auROCs_A1,auROCs_D6)
kD1 = kurtosis(auROCs_A1)
kD5 = kurtosis(auROCs_D5)
kD6 = kurtosis(auROCs_D6)
mean(auROCs_D5),mean(auROCs_A1),mean(auROCs_D6)
median(auROCs_D5),median(auROCs_A1),median(auROCs_D6)
std(auROCs_D5),std(auROCs_A1),std(auROCs_D6)
% [p,diff] = permutationTest(auROCs_A1,auROCs_D6,50000,'plotresult',1)
% [p,diff] = permutationTest(auROCs_A1,auROCs_D5,50000,'plotresult',1)
% [p,diff] = permutationTest(auROCs_D6,auROCs_D5,50000,'plotresult',1)

%%
% Excitator
n1 = exca1; N1 = length(auROCs_A1);
n2 = exca5; N2 = length(auROCs_D5);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)

n1 = exca1; N1 = length(auROCs_A1);
n2 = excb1; N2 = length(auROCs_D6);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)

n1 = exca5; N1 = length(auROCs_D5);
n2 = excb1; N2 = length(auROCs_D6);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)


% Inhibitor
n1 = inha1; N1 = length(auROCs_A1);
n2 = inha5; N2 = length(auROCs_D5);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)

n1 = inha1; N1 = length(auROCs_A1);
n2 = inhb1; N2 = length(auROCs_D6);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)

n1 = inha5; N1 = length(auROCs_D5);
n2 = inhb1; N2 = length(auROCs_D6);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)

%%
SizeDots = 10
figure('Render','painters', 'Position',[10 10 900 300])
hold on
load('database_A1_THESIS');
load('cfg_A1_THESIS.mat');
subplot(1,3,1)
LS_A1 = []; LS_A1_FS = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        i = 0;
        if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse==1
                    i = i + 1;
                end
            end
            if i > 0
                LS_A1 = [LS_A1 exp.expID(idxExp).unit(idxUnit).lifetimeSparseness];
            end
        else
           for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)
                if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse==1
                    i = i + 1;
                end
            end
            if i > 0
                LS_A1_FS = [LS_A1_FS exp.expID(idxExp).unit(idxUnit).lifetimeSparseness];
            end 
       end
    end
end


load('database_A5_THESIS');
load('cfg_A5_THESIS.mat');
subplot(1,3,2)
LS_A5 = []; LS_A5_FS =[];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        i = 0;
       if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1 
                if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse==1
                    i = i + 1;
                end
            end
            if i > 0
                LS_A5 = [LS_A5 exp.expID(idxExp).unit(idxUnit).lifetimeSparseness];
            end
        else
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)
                if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse==1
                    i = i + 1;
                end
            end
            if i > 0
                LS_A5_FS = [LS_A5_FS exp.expID(idxExp).unit(idxUnit).lifetimeSparseness];
            end
        end
   end
end


load('database_B1_THESIS');
load('cfg_B1_THESIS.mat');
LS_B1 = [];LS_B1_FS=[];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        i = 0;
        if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
                if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse==1 %|| exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse==1
                    i = i + 1;
                end
            end
            if i > 0
                LS_B1 = [LS_B1 exp.expID(idxExp).unit(idxUnit).lifetimeSparseness];
            end
        else
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)
                if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse==1 %|| exp.expID(idxExp).unit(idxUnit).event(idxEvent).inhibitoryResponse==1
                    i = i + 1;
                end
            end
           if i > 0
                LS_B1_FS = [LS_B1_FS exp.expID(idxExp).unit(idxUnit).lifetimeSparseness];
            end
        end
   end
end

SizeDots = 10
figure('Render','painters', 'Position',[10 10 900 300])
hold on
subplot(1,3,1)
edges = [0:0.2:1];
histogram(LS_A1, edges, 'FaceColor', cb_d1, 'EdgeColor','none', 'Normalization', 'probability');
hold on
scatter(LS_A1, 0.42 +0.02*rand(1,length(LS_A1)), 'filled', 'MarkerFaceColor', cb_d1,'MarkerEdgeColor','none',...
    'SizeData',SizeDots)
scatter(LS_A1_FS, 0.42 +0.02*rand(1,length(LS_A1_FS)),'MarkerEdgeColor',cb_d1,...
    'SizeData',SizeDots)
ylabel('Fraction of Neurons')
ylim([0 0.47]); yticks(0:0.1:0.4)
title('A1')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')

subplot(1,3,2)
histogram(LS_A5, edges, 'FaceColor', cb_d5, 'EdgeColor','none', 'Normalization', 'probability');
hold on
scatter(LS_A5, 0.42 +0.02*rand(1,length(LS_A5)), 'filled', 'MarkerFaceColor', cb_d5,'MarkerEdgeColor','none',...
    'SizeData',SizeDots)
scatter(LS_A5_FS, 0.42 +0.02*rand(1,length(LS_A5_FS)),'MarkerEdgeColor', cb_d5,...
    'SizeData',SizeDots)
ylim([0 0.47]);yticks(0:0.1:0.3);
title('A5')
xlabel('Lifetime Sparseness')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
ax1 = gca;       
ax1.YAxis.Visible = 'off';

subplot(1,3,3)
histogram(LS_B1, edges, 'FaceColor', cb_d6, 'EdgeColor','none', 'Normalization', 'probability');
hold on
scatter(LS_B1, 0.42 +0.02*rand(1,length(LS_B1)), 'filled', 'MarkerFaceColor', cb_d6,'MarkerEdgeColor','none',...
    'SizeData',SizeDots)
scatter(LS_B1_FS, 0.42 +0.02*rand(1,length(LS_A5_FS)),'MarkerEdgeColor',cb_d6,...
    'SizeData',SizeDots)
ylim([0 0.47]);yticks(0:0.1:0.3);
title('B1')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
ax1 = gca;       
ax1.YAxis.Visible = 'off';
[h,p]=kstest2(LS_A1,LS_A5)
[h,p]=kstest2(LS_A1,LS_B1)
[h,p]=kstest2(LS_B1,LS_A5)

[h,p]=ranksum(LS_A1,LS_A5)
[h,p]= ranksum(LS_A1,LS_B1)
[h,p]=ranksum(LS_B1,LS_A5)

median(LS_A1)
median(LS_A5)
median(LS_B1)
%% Tuning Curves

load('cfg_A1_THESIS.mat');
load('database_A1_THESIS');
nResponsiveEventsFS_A1 = [];
nResponsiveEventsRS_A1 = [];
normTuningCurveFS_A1 = [];
normTuningCurveRS_A1 = [];
zs_TuningCurve_A1 = [];
grad_A1 = [];

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
                nResponsiveEventsFS_A1 = [nResponsiveEventsFS_A1; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveFS_A1 = [normTuningCurveFS_A1; app];
            else
                nResponsiveEventsRS_A1 = [nResponsiveEventsRS_A1; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveRS_A1 = [normTuningCurveRS_A1; app];
                app_grad_A1 = polyfit(1:7,fliplr(app),1);
                grad_A1 = [grad_A1 app_grad_A1(1)];
                
            end
            zs_TuningCurve_A1 = [zs_TuningCurve_A1; zscore(tuningCurve)];
        end
    end
end

load('cfg_A5_THESIS.mat');
load('database_A5_THESIS');
nResponsiveEventsFS_A5 = [];
nResponsiveEventsRS_A5 = [];
normTuningCurveFS_A5 = [];
normTuningCurveRS_A5 = [];
zs_TuningCurve_A5 = [];
grad_A5 = [];

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
                nResponsiveEventsFS_A5 = [nResponsiveEventsFS_A5; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveFS_A5 = [normTuningCurveFS_A5; app];
            else
                nResponsiveEventsRS_A5 = [nResponsiveEventsRS_A5; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveRS_A5 = [normTuningCurveRS_A5; app];
                app_grad_A5 = polyfit(1:7,fliplr(app(1:end)),1);
                grad_A5 = [grad_A5 app_grad_A5(1)];
                
            end
            zs_TuningCurve_A5 = [zs_TuningCurve_A5; zscore(tuningCurve)];
        end
    end
end


load('cfg_B1_THESIS.mat');
load('database_B1_THESIS');
nResponsiveEventsFS_B1 = [];
nResponsiveEventsRS_B1 = [];
normTuningCurveFS_B1 = [];
normTuningCurveRS_B1 = [];
zs_TuningCurve_B1 = [];
grad_B1 = [];

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
                nResponsiveEventsFS_B1 = [nResponsiveEventsFS_B1; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveFS_B1 = [normTuningCurveFS_B1; app];
            else
                nResponsiveEventsRS_B1 = [nResponsiveEventsRS_B1; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveRS_B1 = [normTuningCurveRS_B1; app];
                app_grad_B1 = polyfit(1:7,fliplr(app),1);
                grad_B1 = [grad_B1 app_grad_B1(1)];
                
            end
            zs_TuningCurve_B1 = [zs_TuningCurve_B1; zscore(tuningCurve)];
        end
    end
end
figure('Render','painters', 'Position', [10 10 900 500])
 %tiledlayout(1,3)
 %nexttile([1 2])
KK_A1 = [normTuningCurveRS_A1];
meanKK_A1 = fliplr(mean(KK_A1));
stdKK_A1 = fliplr(std(KK_A1));
% [maxKK, imaxKK] = max(meanKK);
errKK_A1 = stdKK_A1./sqrt(size(KK_A1,1));
k2 = shadedErrorBar(1:7, meanKK_A1, errKK_A1)%,'transparent',1);
set(k2.edge,'LineWidth',.5,'LineStyle','none');
k2.mainLine.LineWidth = 2; k2.mainLine.Color= cb_d1;
k2.patch.FaceAlpha = .3; k2.patch.FaceColor= cb_d1;

hold on

KK_A5 = normTuningCurveRS_A5;
meanKK_A5 = fliplr(mean(KK_A5));
stdKK_A5 = fliplr(std(KK_A5));
% [maxKK, imaxKK] = max(meanKK);
errKK_A5 = stdKK_A5./sqrt(size(KK_A5,1));
k1 = shadedErrorBar(1:7, meanKK_A5, errKK_A5)%,'transparent',1);
set(k1.edge,'LineWidth',.5,'LineStyle','none');
k1.mainLine.LineWidth = 2; k1.mainLine.Color= cb_d5;
k1.patch.FaceAlpha = .3;  k1.patch.FaceColor= cb_d5;

KK_B1 = normTuningCurveRS_B1;
meanKK_B1 = fliplr(mean(KK_B1));
stdKK_B1 = fliplr(std(KK_B1));
% [maxKK, imaxKK] = max(meanKK);
errKK_B1 = stdKK_B1./sqrt(size(KK_B1,1));
k3 = shadedErrorBar(1:7, meanKK_B1, errKK_B1)%,'transparent',1);
set(k3.edge,'LineWidth',.5,'LineStyle','none');
k3.mainLine.LineWidth = 2; k3.mainLine.Color= cb_d6;
k3.patch.FaceAlpha = .3;  k3.patch.FaceColor= cb_d6;

xlabel('Number of Stimuli')
ylabel('Scaled response')
legend([k2.mainLine,k1.mainLine,k3.mainLine],{'A1','A5','B1'})
legend('boxoff')
title('Scaled average tuning curve')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
xticks(1:7)
xlim([0 8]); ylim([-.5 1.1]);
set(gca, 'box', 'off', 'tickDir', 'out')
hold on

clear pp
clear hh
for ii = 1:7
[h,p]=ranksum(fliplr(normTuningCurveRS_A1(:,ii)),fliplr(normTuningCurveRS_A5(:,ii)))
hh(ii)=h
pp(ii)=p
end

clear pp
clear hh
for ii = 1:7
[h,p]=ranksum(fliplr(normTuningCurveRS_A1(:,ii)),fliplr(normTuningCurveRS_B1(:,ii)))
hh(ii)=h
pp(ii)=p
end

clear pp
clear hh
for ii = 1:7
[h,p]=ranksum(fliplr(normTuningCurveRS_A5(:,ii)),fliplr(normTuningCurveRS_B1(:,ii)))
hh(ii)=h
pp(ii)=p
end
%%

nexttile
plot(2, mean(nResponsiveEventsRS_A1), 'o', 'markersize', 8, 'markeredgecolor', cb_d1, 'markerfacecolor', cb_d1)
hold on
%plot(3, mean(nResponsiveEventsFS_A1), 'o', 'markersize', 8, 'markeredgecolor', cb_d1)
hold on
plot(4, mean(nResponsiveEventsRS_A5), 'o', 'markersize', 8, 'markeredgecolor', cb_d5, 'markerfacecolor', cb_d5)
hold on
%plot(5, mean(nResponsiveEventsFS_A5), 'o', 'markersize', 8, 'markeredgecolor', cb_d5)
hold on
plot(6, mean(nResponsiveEventsRS_B1), 'o', 'markersize', 8, 'markeredgecolor', cb_d6, 'markerfacecolor', cb_d6)
hold on
%plot(7, mean(nResponsiveEventsFS_B1), 'o', 'markersize', 8, 'markeredgecolor', cb_d6)
hold on
%errbar(3, mean(nResponsiveEventsFS_A1), std(nResponsiveEventsFS_A1)./sqrt(length(nResponsiveEventsFS_A1)), 'color', cb_d1, 'linewidth', 1); %
hold on
errbar(2, mean(nResponsiveEventsRS_A1), std(nResponsiveEventsRS_A1)./sqrt(length(nResponsiveEventsRS_A1)), 'color', cb_d1, 'linewidth', 1); %
hold on
%errbar(5, mean(nResponsiveEventsFS_A5), std(nResponsiveEventsFS_A5)./sqrt(length(nResponsiveEventsFS_A5)), 'color', cb_d5, 'linewidth', 1); %
hold on
errbar(4, mean(nResponsiveEventsRS_A5), std(nResponsiveEventsRS_A5)./sqrt(length(nResponsiveEventsRS_A5)), 'color', cb_d5, 'linewidth', 1); %
hold on
%errbar(7, mean(nResponsiveEventsFS_B1), std(nResponsiveEventsFS_B1)./sqrt(length(nResponsiveEventsFS_B1)), 'color', cb_d6, 'linewidth', 1); %
hold on
errbar(6, mean(nResponsiveEventsRS_B1), std(nResponsiveEventsRS_B1)./sqrt(length(nResponsiveEventsRS_B1)), 'color', cb_d6, 'linewidth', 1); %
xlim([1 7])
ylim([0 7])
%legend('A1 RS','A1 FS', 'A5 RS', 'A5 FS','B1 RS','B1 FS', 'NumColumns',3)
legend('A1 RS', 'A5 RS','B1 RS')
legend('boxoff')
set(gca, 'XColor', 'w', 'box','off')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
title('Number of responses per neuron')
ylabel('Number of stimuli')

%% Noise

load('cfg_A1_THESIS.mat');
load('database_A1_THESIS');
tempMean = [];
all_NoiseVectors_A1 = [];
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
                all_NoiseVectors_A1(l,:) = zscore(tempMean);
                tempMean = [];
            %end
            eventResponse = 0;
        end
    end
end

load('cfg_A5_THESIS.mat');
load('database_A5_THESIS');
tempMean = [];
all_NoiseVectors_A5 = [];
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
                all_NoiseVectors_A5(l,:) = zscore(tempMean);
                tempMean = [];
            %end
            eventResponse = 0;
        end
    end
end

load('cfg_B1_THESIS.mat');
load('database_B1_THESIS');
tempMean = [];
all_NoiseVectors_B1 = [];
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
                all_NoiseVectors_B1(l,:) = zscore(tempMean);
                tempMean = [];
            %end
            eventResponse = 0;
        end
    end
end

U_A1 = corr(all_NoiseVectors_A1');
U_A5 = corr(all_NoiseVectors_A5');
U_B1 = corr(all_NoiseVectors_B1');

At = U_A1.';
m  = (1:size(At,1)).' >= (1:size(At,2));
vN_A1  = At(m);
vN_A1 = vN_A1(vN_A1~=1);

At = U_A5.';
m  = (1:size(At,1)).' >= (1:size(At,2));
vN_A5  = At(m);
vN_A5 = vN_A5(vN_A5~=1);

At = U_B1.';
m  = (1:size(At,1)).' >= (1:size(At,2));
vN_B1  = At(m);
vN_B1 = vN_B1(vN_B1~=1);


%% Signal Noise
U_A1 = corr(zs_TuningCurve_A1');
U_A5 = corr(zs_TuningCurve_A5');
U_B1 = corr(zs_TuningCurve_B1');

At_A1 = U_A1.';
m  = (1:size(At_A1,1)).' >= (1:size(At_A1,2));
v_A1  = At_A1(m);
v_A1 = v_A1(v_A1~=1);

At_A5 = U_A5.';
m  = (1:size(At_A5,1)).' >= (1:size(At_A5,2));
v_A5  = At_A5(m);
v_A5 = v_A5(v_A5~=1);

At_B1 = U_B1.';
m  = (1:size(At_B1,1)).' >= (1:size(At_B1,2));
v_B1  = At_B1(m);
v_B1 = v_B1(v_B1~=1);

f7 = figure( 'Position', [10 10 400 400])
% tiledlayout(3,12,'TileSpacing','compact')
% nexttile([3 1])
d{1} = v_A1; d{2} = v_A5; d{3} = v_B1
h1 = raincloud_plot(d{1}, 'box_on', 1, 'color', cb_d1, 'alpha', 0.6,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
     'box_col_match', 1);
h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cb_d5, 'alpha', 0.6,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75, 'box_col_match', 1);
 h3 = raincloud_plot(d{3}, 'box_on', 1, 'color', cb_d6, 'alpha', 0.6,...
     'box_dodge', 1, 'box_dodge_amount', .95, 'dot_dodge_amount', 1.15, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';  h3{1, 1}.LineStyle = 'none';
h1{1, 2}.SizeData = 2; h2{1, 2}.SizeData = 2; h3{1, 2}.SizeData = 2;
h1{1, 2}.MarkerFaceAlpha = .7; h2{1, 2}.MarkerFaceAlpha = .7; h3{1, 2}.MarkerFaceAlpha = .7;
title(['Signal Correlations']); ylim([-1.3 1.5])
set(gca,'XLim', [-1 1]); set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
ax = gca;   %or as appropriate
yticks = get(ax, 'YTick'); yticks=yticks(yticks>=0); set(ax, 'YTick', yticks);
hold on 
ylab = ylabel('PDF')
ylab.Position(2) = 0.75
ylab.Position(1) = -1.2

hold off

% nexttile([2 1])
% d{4} = vN_A1; d{5} = vN_A5; d{6} = vN_B1;
% h4 = raincloud_plot(d{4}, 'box_on', 1, 'color', cb_d1, 'alpha', 0.6,...
%      'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
%      'box_col_match', 1);
% h5 = raincloud_plot(d{5}, 'box_on', 1, 'color', cb_d5, 'alpha', 0.6,...
%      'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75, 'box_col_match', 1);
%  h6 = raincloud_plot(d{6}, 'box_on', 1, 'color', cb_d6, 'alpha', 0.6,...
%      'box_dodge', 1, 'box_dodge_amount', .95, 'dot_dodge_amount', 1.15, 'box_col_match', 1);
% h4{1, 1}.LineStyle = 'none'; h5{1, 1}.LineStyle = 'none';  h6{1, 1}.LineStyle = 'none';
% h4{1, 2}.SizeData = 8; h5{1, 2}.SizeData = 8; h6{1, 2}.SizeData = 8;
% h4{1, 2}.MarkerFaceAlpha = .7; h5{1, 2}.MarkerFaceAlpha = .7; h6{1, 2}.MarkerFaceAlpha = .7;
% title(['Noise Correlations']);
% ylim([-5.1 6])
% set(gca,'XLim', [-1 1]);set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
% ax = gca;   %or as appropriate
% yticks = get(ax, 'YTick'); yticks=yticks(yticks>=0); set(ax, 'YTick', yticks);
% yticklabels = get(ax, 'YTickLabel');
% yticklabels{1} = ''; yticklabels{2} = ''; yticklabels{3} = ''; set(ax, 'YTickLabel', yticklabels);
% legend([h4{1} h5{1} h6{1}], {'Odor Set A - Day 1', 'Odor Set A - Day 5', 'Odor Set B - Day 1'});
% h3 = line([0 0], [0 4], 'LineStyle', '--');
% h3.Color = [0 0 0 .5];
% legend boxoff
% ylab = ylabel('PDF')
% ylab.Position(2) = 2
% ylim([-3.25 3.5])

[h1,p1]= ranksum(d{1},d{2});
mean(d{1})
mean(d{2})
median(d{3})

[h2,p2]= ttest2(d{6},d{5});
mean(d{3}),mean(d{4})


mean(d{4})
mean(d{5})
mean(d{6})

[H,P] = ranksum(d{1},d{2})
[H,P] = ranksum(d{1},d{3})
[H,P] = ranksum(d{2},d{3})
%%
figure( 'Position', [10 10 900 300])
%nexttile(brah(idxPeriod-1))
neuron_information_A1=load('SNI_A1.mat')
neuron_information_A1 =neuron_information_A1.neuron_information;
neuron_information_A5=load('SNI_A5.mat')
neuron_information_A5 =neuron_information_A5.neuron_information;
neuron_information_B1=load('SNI_B1.mat')
neuron_information_B1 =neuron_information_B1.neuron_information;
neuron_information_All = [neuron_information_A1;neuron_information_A5; neuron_information_B1]
Area = cellstr(["A1","A5"])
K    = cell(1, 107);
K(:) = {'A1'};
C    = cell(1, 102);
C(:) = {'A5'};
S   = cell(1, 113);
S(:) = {'B1'};
P = [K C S]
% hold on
% scatter(2+randn(1,length(neuron_information_PCx))/10 ,neuron_information_PCx, 12,'filled','MarkerEdgeColor', cb_PCx,'MarkerFaceColor',cb_PCx)
% scatter(2,mean(neuron_information_PCx), 's','k','filled')
vs = violinplot(neuron_information_All,P)
vs(1, 1).ViolinColor = cb_d1;vs(1, 2).EdgeColor = 'none';vs(1, 1).ViolinAlpha = .35;vs(1, 1).ScatterPlot.MarkerFaceAlpha =.5
vs(1, 2).ViolinColor = cb_d5;vs(1, 1).EdgeColor = 'none';vs(1, 2).ViolinAlpha = .35;vs(1, 2).ScatterPlot.MarkerFaceAlpha =.5
vs(1, 3).ViolinColor = cb_d6;vs(1, 3).EdgeColor = 'none';vs(1, 3).ViolinAlpha = .35;vs(1, 3).ScatterPlot.MarkerFaceAlpha =.5
xlim([.5 3.5])
ylim([0 .85])
line([1 3], [0.75 0.75], 'color', 'k', 'linestyle', '-')
line([1 2], [0.68 0.68], 'color', 'k', 'linestyle', '-')
% scatter(5+randn(1,length(neuron_information_PRh))/10 ,neuron_information_PRh, 12, 'filled','MarkerEdgeColor', cb_PRh,'MarkerFaceColor',cb_PRh);
% scatter(5,mean(neuron_information_PRh), 30,'s','k','filled')

ylabel('bits')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
 text(1.45,0.69, '*','FontSize',14);
  text(1.985,0.78, '*','FontSize',14);
 title('Single Unit Mutual Information')
[pni,hni]=ranksum(neuron_information_A1,neuron_information_B1)
[pni,hni]=ranksum(neuron_information_A1,neuron_information_A5)
[pni,hni]=ranksum(neuron_information_A5,neuron_information_B1)
%%

f7 = figure( 'Position', [10 10 900 400])
colormap(brewermap([],'*PuBuGn'));
subplot(1,3,1)
imagesc(At_A1,[-1 1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');set(gca, 'YTick', []);
ax1 = gca;            
ax1.XAxis.Visible = 'off';
set(gca, 'XTick', []);
ylabel('Neuron ID')
axis square
title('A1')
subplot(1,3,2)
imagesc(At_A5,[-1 1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'YTick', []);
ax1 = gca;            
ax1.YAxis.Visible = 'off';
ax1.XAxis.Visible = 'off';
set(gca, 'XTick', []);
set(gca, 'box', 'off', 'tickDir', 'out')%set(gca, 'XAxis', []);
axis square
title('A5')
subplot(1,3,3)
imagesc(At_B1,[-1 1])
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'YTick', []);
ax1 = gca;            
ax1.YAxis.Visible = 'off';
ax1.XAxis.Visible = 'off';
set(gca, 'XTick', []);
set(gca, 'box', 'off', 'tickDir', 'out')%set(gca, 'XAxis', []);
axis square
title('B1')
h=colorbar;
hold off
set(h, 'Location','south','Position', [.20 .12 .6 0.05])
%% DECO SIZE

nNeuroDec = [0:10:100];
nNeuroDec(1) = 1;

% PCx
for i = 1:11
    stringa = ['database_A1_RSEXC_' num2str(nNeuroDec(i)) '_neurons.mat']
    load(stringa)
    performancesA1_mu(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(2)];
    performancesA1_stdev(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(2)];
    performancesA1_mu_succ(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(3)];
    performancesA1_stdev_succ(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(3)];
end



% A5

for i = 1:11
    stringa = ['database_A5_RSEXC_' num2str(nNeuroDec(i)) '_neurons.mat']
    load(stringa)
    performancesA5_mu(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(2)];
    performancesA5_stdev(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(2)];
    performancesA5_mu_succ(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(3)];
    performancesA5_stdev_succ(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(3)];
end



% b1

for i = 1:11
    stringa = ['database_B1_RSEXC_' num2str(nNeuroDec(i)) '_neurons.mat']
    load(stringa)
    performancesB1_mu(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(2)];
    performancesB1_stdev(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(2)];
    performancesB1_mu_succ(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results(3)];
    performancesB1_stdev_succ(i) = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples(3)];
end

f7 = figure( 'Position', [10 10 800 400]);
tiledlayout(1,2)
nexttile
numOfStimuli = 7;
decSetA1 = shadedErrorBar(nNeuroDec, performancesA1_mu, performancesA1_stdev, 'k',1)%,{'color', 'r', 'linewidth', 2});
set(decSetA1 .edge,'LineWidth',.5,'LineStyle','none');
decSetA1.mainLine.LineWidth = 2; decSetA1.mainLine.Color= cb_d1;
decSetA1.patch.FaceAlpha = .3; decSetA1.patch.FaceColor= cb_d1;

hold on

decSetA5 = shadedErrorBar(nNeuroDec, performancesA5_mu, performancesA5_stdev, 'k',1)%,{'color', 'r', 'linewidth', 2});
set(decSetA5 .edge,'LineWidth',.5,'LineStyle','none');
decSetA5.mainLine.LineWidth = 2; decSetA5.mainLine.Color= cb_d5;
decSetA5.patch.FaceAlpha = .3; decSetA5.patch.FaceColor= cb_d5;

decSetB1 = shadedErrorBar(nNeuroDec, performancesB1_mu, performancesB1_stdev, 'k',1)%,{'color', 'r', 'linewidth', 2});
set(decSetB1 .edge,'LineWidth',.5,'LineStyle','none');
decSetB1.mainLine.LineWidth = 2; decSetB1.mainLine.Color= cb_d6;
decSetB1.patch.FaceAlpha = .3; decSetB1.patch.FaceColor= cb_d6;


line([0 nNeuroDec(end)+5], [1/numOfStimuli 1/numOfStimuli], 'color', [.3 .3 .3], 'linestyle', '--')
text(66,1/numOfStimuli-0.04,'Chance Level','color', [.3 .3 .3]);
set(gcf,'color','white', 'PaperPositionMode', 'auto');
ylim([0 1])
xlim([0 nNeuroDec(end)+1]); xticks(nNeuroDec);
ylabel('Accuracy \%')
xlabel('Population Size')
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
title({'Mean Decoding Accuracy of Odor Identity','0-1s Window'})
hold on
legend([decSetA1.mainLine,decSetA5.mainLine,decSetB1.mainLine] ,{'A1','A5','B1'},'Location','northwest')
legend boxoff

nexttile
decSetA1_succ = shadedErrorBar(nNeuroDec, performancesA1_mu_succ, performancesA1_stdev_succ, 'k',1)%,{'color', 'r', 'linewidth', 2});
set(decSetA1_succ.edge,'LineWidth',.5,'LineStyle','none');
decSetA1_succ.mainLine.LineWidth = 2; decSetA1_succ.mainLine.Color= cb_d1; decSetA1_succ.mainLine.LineStyle = '--';
decSetA1_succ.patch.FaceAlpha = .3; decSetA1_succ.patch.FaceColor= cb_d1;

hold on

decSetA5_succ = shadedErrorBar(nNeuroDec, performancesA5_mu_succ, performancesA5_stdev_succ, 'k',1)%,{'color', 'r', 'linewidth', 2});
set(decSetA5_succ.edge,'LineWidth',.5,'LineStyle','none');
decSetA5_succ.mainLine.LineWidth = 2; decSetA5_succ.mainLine.Color= cb_d5; decSetA5_succ.mainLine.LineStyle = '--';
decSetA5_succ.patch.FaceAlpha = .3; decSetA5_succ.patch.FaceColor= cb_d5;

decSetB1_succ = shadedErrorBar(nNeuroDec, performancesB1_mu_succ, performancesB1_stdev_succ, 'k',1)%,{'color', 'r', 'linewidth', 2});
set(decSetB1_succ.edge,'LineWidth',.5,'LineStyle','none');
decSetB1_succ.mainLine.LineWidth = 2; decSetB1_succ.mainLine.Color= cb_d6; decSetB1_succ.mainLine.LineStyle = '--';
decSetB1_succ.patch.FaceAlpha = .3; decSetB1_succ.patch.FaceColor= cb_d6;

line([0 nNeuroDec(end)+5], [1/numOfStimuli 1/numOfStimuli], 'color', [.3 .3 .3], 'linestyle', '--')
text(66,1/numOfStimuli-0.04,'Chance Level','color', [.3 .3 .3]);
set(gcf,'color','white', 'PaperPositionMode', 'auto');
ylim([0 1])
xlim([0 nNeuroDec(end)+1]); xticks(nNeuroDec);
ylabel('Accuracy \%')
xlabel('Population Size')
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
title({'Mean Decoding Accuracy of Odor Identity','1-2s Window'})


%% DECO TIMECOURSE
load('A1_TimeCourse_200by100ms')
figure( 'Position', [10 10 800 400])
%subplot(2,1,1);
numOfStimuli = 7;
tcperf_mu_A1 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_A1 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
k = shadedErrorBar(1:49, tcperf_mu_A1, tcperf_stdev_A1,'k',1);
hold on
set(k.edge,'LineWidth',.5,'LineStyle','none');
k.mainLine.LineWidth = 2; k.mainLine.Color= cb_d1;
k.patch.FaceAlpha = .3; k.patch.FaceColor= cb_d1;
load('A5_TimeCourse_200by100ms')
tcperf_mu_A5 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_A5 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
z = shadedErrorBar(1:49, tcperf_mu_A5, tcperf_stdev_A5,'k',1);
hold on
set(z.edge,'LineWidth',.5,'LineStyle','none');
z.mainLine.LineWidth = 2; z.mainLine.Color= cb_d5;
z.patch.FaceAlpha = .3; z.patch.FaceColor= cb_d5;
load('B1_TimeCourse_200by100ms')
tcperf_mu_B1 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_B1 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
t = shadedErrorBar(1:49, tcperf_mu_B1, tcperf_stdev_B1,'k',1);
hold on
set(t.edge,'LineWidth',.5,'LineStyle','none');
t.mainLine.LineWidth = 2; t.mainLine.Color= cb_d6;
t.patch.FaceAlpha = .3; t.patch.FaceColor= cb_d6;

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

legend([k.mainLine,z.mainLine,t.mainLine] ,{'A1','A5','B1'},'Location','northwest')
legend boxoff
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
%%
subplot(2,1,2);
numOfStimuli = 7;
load('A5+M1_TimeCourse_200by100ms')
tcperf_mu_A5 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_A5 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
z = shadedErrorBar(1:49, tcperf_mu_A5, tcperf_stdev_A5,'k',1);
hold on
set(z.edge,'LineWidth',.5,'LineStyle','none');
z.mainLine.LineWidth = 2; z.mainLine.Color= cb_d5;
z.patch.FaceAlpha = .3; z.patch.FaceColor= cb_d5;
subplot(2,1,2);
load('A5-M1_TimeCourse_200by100ms')
tcperf_mu_A5m = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_A5m = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
zm = shadedErrorBar(1:49, tcperf_mu_A5m, tcperf_stdev_A5m,'k',1);
hold on
set(zm.edge,'LineWidth',.5,'LineStyle','none');
zm.mainLine.LineWidth = 2; zm.mainLine.Color= cb_d5; zm.mainLine.LineStyle='--'
zm.patch.FaceAlpha = .3; zm.patch.FaceColor= cb_d5;
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

legend([z.mainLine,zm.mainLine] ,{'A5 with M1','A5 without M1'},'Location','northwest')
legend boxoff
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
%% PCA

load('database_A1_THESIS.mat');
load('cfg_A1_THESIS.mat');
a = 0;
s = 0;
sR_acrossTrials = [];
pseudoTrials = [];
rep = [];
matPCA = [];
Trial = 5;
dotSize = 40;
AlphaValue = .8;
EdgeAlpha = 0.5;

% Average across three trials
for idxExp = 1:6 %:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        s = s + 1;
        count = 0;
        a = 0;
        for idxEvents = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
            for idxTrial = 1:length(exp.expID(idxExp).unit(idxUnit).event(idxEvents).spikeResponse)
                sR = exp.expID(idxExp).unit(idxUnit).event(idxEvents).spikeResponse(idxTrial);
                a = a + 1;
                matPCA(a,s) = exp.expID(idxExp).unit(idxUnit).event(idxEvents).spikeResponse(idxTrial);
            end
        end
    end
end

s= 0;
z= 0;
ind = [];
for idxExp = 1:6 %:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        s = s + 1;
        for idxEvents = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
            if exp.expID(idxExp).unit(idxUnit).event(idxEvents).excitatoryResponse == 1
                z = z + 1;
            end
        end
        if z > 0
            ind = [ind s];
        end
        z = 0;
    end
end

switch Trial
    case 5
        excmatPCA_A1 = matPCA(:,ind);
        B1_A1 = squeeze(mean(reshape(excmatPCA_A1,3,35,[]),1)); %49
    case 15
        excmatPCA_A1 = matPCA(:,ind);
        B1_A1 = excmatPCA_A1;
    case 3
        excmatPCA_A1 = matPCA(:,ind);
        B1_A1 = squeeze(mean(reshape(excmatPCA_A1,5,21,[]),1));
end

Bz_A1 = zscore(B1_A1,[],1);
[coeff,score_A1,latent,tsquared, explained] = pca(Bz_A1);

A1 = figure('Render','painters', 'Position', [10 10 1200 1200])
for i = 1:7
    rep(i,:) = (1+Trial*(i-1)):(Trial + Trial*(i-1));
    scatter3(score_A1(rep(i,:),1), score_A1(rep(i,:),2), score_A1(rep(i,:),3), dotSize, cb_pca(i,:), 'filled', 'MarkerFaceAlpha', AlphaValue);
    hold on
end
hold off
zlim([-6 6]);
ylim([-6 6]);
xlim([-6 6]);
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('Odor Set A - Novel');
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')

%%
load('database_A5_THESIS.mat');
load('cfg_A5_THESIS.mat');
a = 0;
s = 0;
sR_acrossTrials = [];
pseudoTrials = [];
matPCA = [];
B1 = [];
Bz = [];

for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        s = s + 1;
        count = 0;
        a = 0;
        for idxEvents = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
            for idxTrial = 1:length(exp.expID(idxExp).unit(idxUnit).event(idxEvents).spikeResponse)
                sR = exp.expID(idxExp).unit(idxUnit).event(idxEvents).spikeResponse(idxTrial);
                a = a + 1;
                matPCA(a,s) = exp.expID(idxExp).unit(idxUnit).event(idxEvents).spikeResponse(idxTrial);
            end
        end
    end
end

s= 0;
z= 0;
ind = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        s = s + 1;
        for idxEvents = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
            if exp.expID(idxExp).unit(idxUnit).event(idxEvents).excitatoryResponse == 1
                z = z + 1;
            end
        end
        if z > 0
            ind = [ind s];
        end
        z = 0;
    end
end

switch Trial
    case 5
        excmatPCA_A5 = matPCA(:,ind);
        B1_A5 = squeeze(mean(reshape(excmatPCA_A5, 3,35,[]),1)); %38
    case 15
        excmatPCA_A5 = matPCA(:,ind);
        B1_A5 = excmatPCA_A5;
    case 1
        excmatPCA_A5 = matPCA(:,ind);
        B1_A5 = squeeze(mean(reshape(excmatPCA,[5,21,[]]),1));
end

Bz_A5 = zscore(B1_A5,[],1);
[coeff,score_A5,latent,tsquared, explained] = pca(Bz_A5);

A5 = figure('Render','painters', 'Position', [10 10 1200 1200]);
rep = [];

for i = 1:7
    rep(i,:) = (1+Trial*(i-1)):(Trial + Trial*(i-1));
    scatter3(score_A5(rep(i,:),1), score_A5(rep(i,:),2), score_A5(rep(i,:),3), dotSize, cb_pca(i,:), 'filled','MarkerFaceAlpha', AlphaValue );
    hold on
end

zlim([-6 6]);
ylim([-6 6]);
xlim([-6 6]);
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('Odor Set A - Familiar');
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')

%%
load('database_B1_THESIS.mat');
load('cfg_B1_THESIS.mat');
a = 0;
s = 0;
sR_acrossTrials = [];
pseudoTrials = [];
matPCA = [];
B1 = [];
Bz = [];

% Average across three trials
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        s = s + 1;
        count = 0;
        a = 0;
        for idxEvents = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
            for idxTrial = 1:length(exp.expID(idxExp).unit(idxUnit).event(idxEvents).spikeResponse)
                sR = exp.expID(idxExp).unit(idxUnit).event(idxEvents).spikeResponse(idxTrial);
                a = a + 1;
                matPCA(a,s) = exp.expID(idxExp).unit(idxUnit).event(idxEvents).spikeResponse(idxTrial);
            end
        end
    end
end

s= 0;
z= 0;
ind = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        s = s + 1;
        for idxEvents = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
            if exp.expID(idxExp).unit(idxUnit).event(idxEvents).excitatoryResponse == 1
                z = z + 1;
            end
        end
        if z > 0
            ind = [ind s];
        end
        z = 0;
    end
end

switch Trial
    case 5
        excmatPCA_B1 = matPCA(:,ind);
        B1_B1 = squeeze(mean(reshape(excmatPCA_B1,3,35,[]),1)); %49
    case 15
        excmatPCA_B1 = matPCA(:,ind);
        B1_B1 = excmatPCA_B1;
    case 3
        excmatPCA_B1 = matPCA(:,ind);
        B1_B1 = squeeze(mean(reshape(excmatPCA_B1,5,21,[]),1));
end

Bz_B1 = zscore(B1_B1,[],1);
[coeff,score_B1,latent,tsquared, explained] = pca(Bz_B1);

B1 = figure('Render','painters', 'Position', [10 10 1200 1200]);
for i = 1:7
    rep(i,:) = (1+Trial*(i-1)):(Trial + Trial*(i-1));
    scatter3(score_B1(rep(i,:),1), score_B1(rep(i,:),2), score_B1(rep(i,:),3), dotSize, cb_pca(i,:), 'filled', 'MarkerFaceAlpha', AlphaValue);
    hold on
end

hold off
zlim([-6 6]);
ylim([-6 6]);
xlim([-6 6]);
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('Odor Set B - Novel');
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')

%%
pdClusters_A5 = zeros(1,size(score_A5,1));
for i = 1:size(score_A5,1)
    excl = find(~any(rep==i,2));
    pdClusters_A5(i) = mean(pdist2(score_A5(i, 1:3), score_A5(rep(excl,:), 1:3)));
end

pdClusters_mean_A5 = mean(pdClusters_A5);


pdClusters_A1 = zeros(1,size(score_A1,1));
for i = 1:size(score_A1,1)
    excl = find(~any(rep==i,2));
    pdClusters_A1(i) = mean(pdist2(score_A1(i, 1:3), score_A1(rep(excl,:), 1:3)));
end

pdClusters_mean_A1 = mean(pdClusters_A1);


pdClusters_B1 = zeros(1,size(score_B1,1));
for i = 1:size(score_B1,1)
    excl = find(~any(rep==i,2));
    pdClusters_B1(i) = mean(pdist2(score_B1(i, 1:3), score_B1(rep(excl,:), 1:3)));
end

pdClusters_mean_B1 = mean(pdClusters_B1);

%%
Std_meanClusters_A5 = [];
Std_meanClusters_A1 = [];
Std_meanClusters_B1 = [];
Coordinates_meanClusters_A5 = [];
Coordinates_meanClusters_A1 = [];
Coordinates_meanClusters_B1 = [];


for i = 1:7
    Coordinates_meanClusters_A5(i,:) = mean(score_A5(rep(i,:), 1:3));
    Std_meanClusters_A5(i,:) = std(score_A5(rep(i,:), 1:3),[], 1);
    Coordinates_meanClusters_A1(i,:) = mean(score_A1(rep(i,:), 1:3));
    Std_meanClusters_A1(i,:) = std(score_A1(rep(i,:), 1:3),[], 1);
    Coordinates_meanClusters_B1(i,:) = mean(score_B1(rep(i,:), 1:3));
    Std_meanClusters_B1(i,:) = std(score_B1(rep(i,:), 1:3),[], 1);
end

C_A5 = Coordinates_meanClusters_A5;
R_A5 = Std_meanClusters_A5;
figure(A5);

for ii = 1:7
    hold on
    [x,y,z] = ellipsoid(C_A5(ii,1),C_A5(ii,2),C_A5(ii,3),R_A5(ii,1),R_A5(ii,2),R_A5(ii,3), 20);
    surf(x,y,z, 'EdgeColor', [1 1 1], 'FaceAlpha', EdgeAlpha, 'FaceColor', cb_pca(ii,:))
    hold off
end

figure(A1)
C_A1 = Coordinates_meanClusters_A1;
R_A1 = Std_meanClusters_A1;

for ii = 1:7
    hold on
    [x,y,z] = ellipsoid(C_A1(ii,1),C_A1(ii,2),C_A1(ii,3),R_A1(ii,1),R_A1(ii,2),R_A1(ii,3), 20);
    surf(x,y,z, 'EdgeColor', [1 1 1], 'FaceAlpha', EdgeAlpha, 'FaceColor', cb_pca(ii,:))
    hold off
end

figure(B1)
C_B1 = Coordinates_meanClusters_B1;
R_B1 = Std_meanClusters_B1;

for ii = 1:7
    hold on
    [x,y,z] = ellipsoid(C_B1(ii,1),C_B1(ii,2),C_B1(ii,3),R_B1(ii,1),R_B1(ii,2),R_B1(ii,3), 20);
    surf(x,y,z, 'EdgeColor', [1 1 1], 'EdgeAlpha', EdgeAlpha, 'FaceColor', cb_pca(ii,:))
    hold off
end
%%

for i=1:25000
    permutations = randperm(15);
    orderedPerms(i,:) = permutations;
    excmatPCA__A1_shuff = reshape(excmatPCA_A1, 15, 7, []);
    excmatPCA__A1_shuff = excmatPCA__A1_shuff(permutations,:,:);
    Bperm_A1 = squeeze(mean(reshape(excmatPCA__A1_shuff,3,35,[]),1));
    Bz_A1 = zscore(Bperm_A1,[],1);
    [coeff,score_perm_A1,latent,tsquared, explained] = pca(Bz_A1);
    
    excmatPCA__A5_shuff = reshape(excmatPCA_A5, 15, 7, []);
    excmatPCA__A5_shuff = excmatPCA__A5_shuff(orderedPerms(i,:),:,:);
    Bperm_A5 = squeeze(mean(reshape(excmatPCA__A5_shuff,3,35,[]),1));
    Bz_A5 = zscore(Bperm_A5,[],1);
    [coeff,score_perm_A5,latent,tsquared, explained] = pca(Bz_A5);
    
    excmatPCA__B1_shuff = reshape(excmatPCA_B1, 15, 7, []);
    excmatPCA__B1_shuff = excmatPCA__B1_shuff(orderedPerms(i,:),:,:);
    Bperm_B1 = squeeze(mean(reshape(excmatPCA__B1_shuff,3,35,[]),1));
    Bz_B1 = zscore(Bperm_B1,[],1);
    [coeff,score_perm_B1,latent,tsquared, explained] = pca(Bz_B1);
    
    for x = 1:7
        pdOdor_shuff_A1(x,:) = mean(pdist(score_perm_A1(rep(x,:), 1:3)));
        CoG_shuff_A1(x,:,i) = mean(score_perm_A1(rep(x,:), 1:3));
        pdOdor_shuff_A5(x,:) = mean(pdist(score_perm_A5(rep(x,:), 1:3)));
        CoG_shuff_A5(x,:,i) = mean(score_perm_A5(rep(x,:), 1:3));
        pdOdor_shuff_B1(x,:) = mean(pdist(score_perm_B1(rep(x,:), 1:3)));
        CoG_shuff_B1(x,:,i) = mean(score_perm_B1(rep(x,:), 1:3));
    end
    
    pdOdor_meanshuff_A1(i) = mean(pdOdor_shuff_A1);
    pdCoG_meansuff_A1(i) = mean(pdist(CoG_shuff_A1(:,:,i)));
    pdOdor_meanshuff_A5(i) = mean(pdOdor_shuff_A5);
    pdCoG_meansuff_A5(i) = mean(pdist(CoG_shuff_A5(:,:,i)));
    pdOdor_meanshuff_B1(i) = mean(pdOdor_shuff_B1);
    pdCoG_meansuff_B1(i) = mean(pdist(CoG_shuff_B1(:,:,i)));
    
%     pdClusters_shuff_A1 = zeros(1,size(score_perm_A1,1));
%     for z = 1:size(score_perm_A1,1)
%         excl = find(~any(rep==z,2));
%         pdClusters_shuff_A1(z) = mean(pdist2(score_perm_A1(z, 1:3), score_perm_A1(rep(excl,:), 1:3)));
%     end
%     pdClusters_meanshuff_A1(i) = mean(pdClusters_shuff_A1);
end

pdCluster = figure;
histogram(pdOdor_meanshuff_A1, [2:0.1:8], 'Normalization', 'probability');
hold on
histogram(pdOdor_meanshuff_A5, [2:0.1:8], 'Normalization', 'probability');
histogram(pdOdor_meanshuff_B1, [2:0.1:8], 'Normalization', 'probability');
mean(pdOdor_meanshuff_A1)
mean(pdOdor_meanshuff_A5)
title('Average Distance Within the Clusters')
legend('Day 1', 'Day 5','Nov')
legend('boxoff')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')

pdCoG = figure
histogram(pdCoG_meansuff_A1, [2:0.1:6], 'Normalization', 'probability');
hold on
histogram(pdCoG_meansuff_A5, [2:0.1:6], 'Normalization', 'probability');
histogram(pdCoG_meansuff_B1, [2:0.1:6], 'Normalization', 'probability');
mean(pdCoG_meansuff_A1)
mean(pdCoG_meansuff_A5)
title('Average Distance Between the Centers of Mass of the Clusters')
legend('Day 1', 'Day 5', 'Nov')
legend('boxoff')
[h,p]=ttest2(pdCoG_meansuff_A1,pdCoG_meansuff_A5);
[h,p]=ttest2(pdOdor_meanshuff_A1,pdOdor_meanshuff_A5);
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')

%% 

pdCluster = figure('Render','painters', 'Position', [10 10 1800 1200]);;
subplot(1,2,1)
d{1} =  pdOdor_meanshuff_A1; d{2} = pdOdor_meanshuff_A5; d{3} = pdOdor_meanshuff_B1;
h1 = raincloud_plot(d{1},.5,'box_on', 1, 'color', cb_d1, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .25,...
    'box_col_match', 1);
h2 = raincloud_plot(d{2},.5, 'box_on', 1, 'color', cb_d5, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .45, 'dot_dodge_amount', .55, 'box_col_match', 1);
h3 = raincloud_plot(d{3},.5,'box_on', 1, 'color', cb_d6, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .85,...
    'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none'; h3{1, 1}.LineStyle = 'none';
h1{1, 2}.SizeData = 1; h2{1, 2}.SizeData = 1; h3{1, 2}.SizeData = 1;
ylabel('PDF'); xlabel('Euclidean Distance')
ylim([-2 3])
xlim([1.5 6.5]); yticks([0:0.5:3])
% legend([h1{1},h2{1},h3{1}],'Odor Set A - Day 1','Odor Set A - Day 5', 'Odor Set B - Day 1')
% legend boxoff
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');
title('Average Distance Within the Clusters') 

subplot(1,2,2)
d{1} =  pdCoG_meansuff_A1; d{2} = pdCoG_meansuff_A5; d{3} = pdCoG_meansuff_B1;
h1 = raincloud_plot(d{1},'box_on', 1, 'color', cb_d1, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .25,...
    'box_col_match', 1);
h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cb_d5, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .45, 'dot_dodge_amount', .55, 'box_col_match', 1);
h3 = raincloud_plot(d{3},'box_on', 1, 'color', cb_d6, 'alpha', 0.7,...
    'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .85,...
    'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none'; h3{1, 1}.LineStyle = 'none';
h1{1, 2}.SizeData = 1; h2{1, 2}.SizeData = 1; h3{1, 2}.SizeData = 1;
ylim([-2 3]);
xlim([1 6]); yticks([0:0.5:3]); xlabel('Euclidean Distance')
title('Average Distance Between the Centers of Mass of the Clusters')
legend([h1{1},h2{1},h3{1}],'Odor Set A - Day 1','Odor Set A - Day 5', 'Odor Set B - Day 1')
legend('boxoff')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')  


[h1,p1]=ttest2(pdCoG_meansuff_B1,pdCoG_meansuff_A5);
[h2,p2]=ttest2(pdOdor_meanshuff_B1,pdOdor_meanshuff_A5);

%% Tuning Curve

hold on
load('database_NovFam');
load('cfg_NovFam');
nResponsiveEventsFS_A1 = [];
nResponsiveEventsRS_A1 = [];
normTuningCurveFS_A1 = [];
normTuningCurveRS_A1 = [];
zs_TuningCurve_A1 = [];

for idxExp = 1:6
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        eventResponse = 0;
        tuningCurve = nan(1, size(exp.expID(idxExp).unit(idxUnit).event, 2) -1);
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2) -1
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
            tuningCurve(idxEvent) = mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse);    
        end
        if eventResponse > 0
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'FS')
                nResponsiveEventsFS_A1 = [nResponsiveEventsFS_A1; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveFS_A1 = [normTuningCurveFS_A1; app];
            else
                nResponsiveEventsRS_A1 = [nResponsiveEventsRS_A1; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveRS_A1 = [normTuningCurveRS_A1; app];
                zs_TuningCurve_A1 = [zs_TuningCurve_A1; zscore(tuningCurve)];
            end
        end
    end
end
nResponsiveEventsFS_A5 = [];
nResponsiveEventsRS_A5 = [];
normTuningCurveFS_A5 = [];
normTuningCurveRS_A5 = [];
zs_TuningCurve_A5 = [];
for idxExp = 7:12
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        eventResponse = 0;
        tuningCurve = nan(1, size(exp.expID(idxExp).unit(idxUnit).event, 2) -1);
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2) -1
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
            tuningCurve(idxEvent) = mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse);    
        end
        if eventResponse > 0
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'FS')
                nResponsiveEventsFS_A5 = [nResponsiveEventsFS_A5; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveFS_A5 = [normTuningCurveFS_A5; app];
            else
                nResponsiveEventsRS_A5 = [nResponsiveEventsRS_A5; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveRS_A5 = [normTuningCurveRS_A5; app];
                zs_TuningCurve_A5 = [zs_TuningCurve_A5; zscore(tuningCurve)];
                
            end
        end
    end
end

nResponsiveEventsFS_B1 = [];
nResponsiveEventsRS_B1 = [];
normTuningCurveFS_B1 = [];
normTuningCurveRS_B1 = [];
zs_TuningCurve_B1 = [];
for idxExp = 13:18
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        eventResponse = 0;
        tuningCurve = nan(1, size(exp.expID(idxExp).unit(idxUnit).event, 2) -1);
        for idxEvent = 1:size(exp.expID(idxExp).unit(idxUnit).event, 2) -1
            if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                eventResponse = eventResponse + 1;
            end
            tuningCurve(idxEvent) = mean(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse);    
        end
        if eventResponse > 0
            if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'FS')
                nResponsiveEventsFS_B1 = [nResponsiveEventsFS_B1; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveFS_B1 = [normTuningCurveFS_B1; app];
            else
                nResponsiveEventsRS_B1 = [nResponsiveEventsRS_B1; eventResponse];
                app = tuningCurve ./ max(tuningCurve);
                app = sort(app);
                normTuningCurveRS_B1 = [normTuningCurveRS_B1; app];
                zs_TuningCurve_B1 = [zs_TuningCurve_B1; zscore(tuningCurve)];
            end
        end
    end
end
normTuningCurveRS_A1(20,:) = [];
normTuningCurveRS_B1(14,:) = [];
figure
plot(1:7,fliplr(mean(normTuningCurveRS_A1)), 'color', cb_d1, 'LineWidth', 2)
hold on
plot(1:7,fliplr(mean(normTuningCurveRS_A5)), 'color', cb_d5, 'LineWidth', 2)
plot(1:7,fliplr(mean(normTuningCurveRS_B1)), 'color', cb_d6, 'LineWidth', 2)
xlabel('Number of Stimuli')
ylabel('Scaled response')
legend('A1','A5','B1')
legend('boxoff')
title('Scaled average tuning curve')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
xticks(1:7)
xlim([0 8]); ylim([-.5 1.1]);
set(gca, 'box', 'off', 'tickDir', 'out')
hold on

% error bars

%% 
figure;
plot(2, mean(nResponsiveEventsRS_A1), 'o', 'markersize', 8, 'markeredgecolor', cb_d1, 'markerfacecolor', cb_d1)
hold on
plot(3, mean(nResponsiveEventsFS_A1), 'o', 'markersize', 8, 'markeredgecolor', cb_d1)
hold on
plot(4, mean(nResponsiveEventsRS_A5), 'o', 'markersize', 8, 'markeredgecolor', cb_d5, 'markerfacecolor', cb_d5)
hold on
plot(5, mean(nResponsiveEventsFS_A5), 'o', 'markersize', 8, 'markeredgecolor', cb_d5)
hold on
plot(6, mean(nResponsiveEventsRS_B1), 'o', 'markersize', 8, 'markeredgecolor', cb_d6, 'markerfacecolor', cb_d6)
hold on
plot(7, mean(nResponsiveEventsFS_B1), 'o', 'markersize', 8, 'markeredgecolor', cb_d6)
hold on
errbar(3, mean(nResponsiveEventsFS_A1), std(nResponsiveEventsFS_A1)./sqrt(length(nResponsiveEventsFS_A1)), 'color', cb_d1, 'linewidth', 1); %
hold on
errbar(2, mean(nResponsiveEventsRS_A1), std(nResponsiveEventsRS_A1)./sqrt(length(nResponsiveEventsRS_A1)), 'color', cb_d1, 'linewidth', 1); %
hold on
errbar(5, mean(nResponsiveEventsFS_A5), std(nResponsiveEventsFS_A5)./sqrt(length(nResponsiveEventsFS_A5)), 'color', cb_d5, 'linewidth', 1); %
hold on
errbar(4, mean(nResponsiveEventsRS_A5), std(nResponsiveEventsRS_A5)./sqrt(length(nResponsiveEventsRS_A5)), 'color', cb_d5, 'linewidth', 1); %
hold on
errbar(7, mean(nResponsiveEventsFS_B1), std(nResponsiveEventsFS_B1)./sqrt(length(nResponsiveEventsFS_B1)), 'color', cb_d6, 'linewidth', 1); %
hold on
errbar(6, mean(nResponsiveEventsRS_B1), std(nResponsiveEventsRS_B1)./sqrt(length(nResponsiveEventsRS_B1)), 'color', cb_d6, 'linewidth', 1); %
xlim([1 8])
ylim([0 7])
legend('A1 RS','A1 FS', 'A5 RS', 'A5 FS','B1 RS','B1 FS', 'NumColumns',3)
legend('boxoff')
set(gca, 'XColor', 'w', 'box','off')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out', 'fontname', 'helvetica', 'fontsize', 14)
title('Number of responses per neuron')
ylabel('Number of stimuli')

%% PSTH

FaceAlpha = 0.2;
hold on
load('cfg_A1_THESIS.mat');
load('database_A1_THESIS');
A1z = [];
i = 0;
allA1z  = [];
kernel = 0.05;
endlength = 820;
z=0;
AvgN_allA1z =[];  S_A1z=[]; 

for idxExp = 1:size(exp.expID,2)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            z = z+1;q = 0;
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)
                if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                    i = i + 1;
                    q = q + 1;
                    [R,t] = psth(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikesPerTrial, ...
                        kernel,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                    meanBSL_A1(i) = mean(R(21:endlength)); stdBSL_A1(i) = std(R(21:endlength));
                    A1z = (R - meanBSL_A1(i))./stdBSL_A1(i);
                    allA1z(i,:) = A1z;
                    AvgN_allA1z(i,:) = A1z;
                end
            end
            if q > 0
            S_A1z = [S_A1z; mean(AvgN_allA1z,1)];
            AvgN_allA1z=[];
            end
        end
    end
end

mean_A1 = mean(allA1z);
std_A1 = std(allA1z);
%[max_A1, imax_A1] = max(mean_A1);
err_A1 = std_A1./sqrt(size(allA1z,1));
ampl_A1 = max(S_A1z,[],2)

load('cfg_A5_THESIS.mat');
load('database_A5_THESIS');
A5z = [];
i = 0;
allA5z  = [];
z=0;
AvgN_allA5z =[];  S_A5z=[]; 

for idxExp = 1:size(exp.expID,2)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            z = z+1;q = 0;
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)
                if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                    i = i + 1;
                    q = q + 1;
                    [R,t] = psth(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikesPerTrial, ...
                        kernel,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                    meanBSL_A5(i) = mean(R(21:endlength)); stdBSL_A5(i) = std(R(21:endlength));
                    A5z = (R - meanBSL_A5(i))./stdBSL_A5(i);
                    allA5z(i,:) = A5z;
                    AvgN_allA5z(i,:) = allA5z(i,:);
                end
            end
            if q > 0
                S_A5z = [S_A5z; mean(AvgN_allA5z,1)];
                AvgN_allA5z=[];
            end
        end
    end
end

mean_A5 = mean(allA5z);
std_A5 = std(allA5z);
% [max_A1, imax_A1] = max(mean_A1);
err_A5 = std_A5./sqrt(size(allA5z,1));
ampl_A5 = max(S_A5z,[],2)

load('cfg_B1_THESIS.mat');
load('database_B1_THESIS');
B1z = [];
i = 0;
allB1z  = [];
z=0;
AvgN_allB1z =[];  S_B1z=[]; 
for idxExp = 1:size(exp.expID,2)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
            z = z+1;q = 0;
            for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)
                if exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                    i = i + 1;
                    q = q + 1;
                    [R,t] = psth(exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikesPerTrial, ...
                        kernel,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                    meanBSL_B1(i) = mean(R(21:endlength)); stdBSL_B1(i) = std(R(21:endlength));
                    B1z = (R - meanBSL_B1(i))./stdBSL_B1(i);
                    allB1z(i,:) = B1z;
                    
                    AvgN_allB1z(i,:) = B1z;
                end
            end
            if q > 0
                S_B1z = [S_B1z; mean(AvgN_allB1z,1)];
                AvgN_allB1z=[];
            end
        end
    end
end

 mean_B1 = mean(allB1z);
 std_B1 = std(allB1z);
% [max_A1, imax_A1] = max(mean_A1);
 err_B1 = std_B1./sqrt(size(allB1z,1));
 ampl_B1 = max(S_B1z,[],2)
 
load('cfg_PRh_THESIS'); load('database_PRh_THESIS');
p = 0; indices = []; eventResponseAll = []; i = 0; ZZ = []; RPRh=[]; v_PRh=0; auROC_PRh = []; auROC_PRh_overall = [];
vv = 0; RPRh_anti = []; ZZ_anti = []; ZZ_peaks = []; ZZ_peaks_indx = [];

for idxExp = 1:size(exp.expID,2)
    for idxUnit = 1:size(exp.expID(idxExp).unit, 2)
        for idxEventType = 1:length(exp.expID(idxExp).unit(idxUnit).event)
            if exp.expID(idxExp).unit(idxUnit).event(idxEventType).excitatoryResponse == 1
                v_PRh = v_PRh + 1;
                [R,t] = psth(exp.expID(idxExp).unit(idxUnit).event(idxEventType).spikesPerTrial, ...
                    kernel,'n', [-cfg.preEventOnset cfg.postEventOnset], 0);
                meanBSLPRh(v_PRh) = mean(R(21:endlength)); stdBSLPRh(v_PRh) = std(R(21:endlength));
                RPRh = [RPRh; R];
                Z = (R - meanBSLPRh(v_PRh))./stdBSLPRh(v_PRh);
                ZZ = [ZZ ; Z];
            end
        end
    end
end
 mean_PRH = mean(ZZ);
 std_PRH = std(ZZ);
% [max_A1, imax_A1] = max(mean_A1);
 err_PRH = std_PRH./sqrt(size(ZZ,1));
figure('Render','painters', 'Position', [10 10 1400 800])

 a1 = shadedErrorBar(t, mean_A1, err_A1, 'k',1);
set(a1.edge,'LineWidth',.5,'LineStyle','none');
a1.mainLine.LineWidth = 1.5; a1.mainLine.Color= cb_d1;
a1.patch.FaceAlpha = FaceAlpha; a1.patch.FaceColor= cb_d1;
hold on
a5 = shadedErrorBar(t, mean_A5, err_A5, 'k',1);
set(a5.edge,'LineWidth',.5,'LineStyle','none');
a5.mainLine.LineWidth = 1.5;a5.mainLine.Color= cb_d5;
a5.patch.FaceAlpha = FaceAlpha; a5.patch.FaceColor= cb_d5;
     
hold on
b1 = shadedErrorBar(t, mean_B1, err_B1, 'k',1);
set(b1.edge,'LineWidth',.5,'LineStyle','none');
b1.mainLine.LineWidth = 1.5;b1.mainLine.Color= cb_d6;
b1.patch.FaceAlpha = FaceAlpha; b1.patch.FaceColor= cb_d6;
hold on

prh = plot(t, mean_PRH, 'k--');

title('Average PSTH Excitatory Responses of RS Neurons');
xlabel('Time (s - Onset Aligned)');
ylabel('SD over baseline');
xlim([-1 4]); ylim([-1 6]);
hold on
M = -.5*ones(1,3);
ax = plot(0:2, M, 'k', 'LineWidth', 3);
ax.Clipping = 'off';
hold on
legend([a1.mainLine,a5.mainLine,b1.mainLine, prh],'A1','A5','B1', 'Avg PRh');
legend('boxoff')
%xline(t(imaxKK), 'r--');
%xline(t(imaxZZ), 'b--');
set(gcf,'color','white', 'PaperPositionMode', 'auto');

%% Response Intensity

%% Reliability

load('database_A1_Long');
load('cfg_A1_Long.mat');
l = 0;
vect_A1 = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
             if  exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                l = l + 1;
                vect_A1(l,:) = exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse;
             end
        end
        end
    end
end

fanoF_A1 = var(vect_A1,[],2)./mean(vect_A1,2);
fastTr_A1 = mean(vect_A1(:,1:5),2);
lastTr_A1 = mean(vect_A1(:,6:10),2);

load('database_A5_Long');
load('cfg_A5_Long.mat');
l = 0; vect_A5 = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
            if  exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                l = l + 1;
                vect_A5(l,:) = exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse;
            end
        end
        end
    end
end

fanoF_A5 = var(vect_A5,[],2)./mean(vect_A5,2);
fastTr_A5 = mean(vect_A5(:,1:5),2);
lastTr_A5 = mean(vect_A5(:,6:10),2);

load('database_B1_Long');
load('cfg_B1_Long.mat');
l = 0;
vect_B1 = [];
for idxExp = 1:size(cfg.folders.folderList,1)
    for idxUnit = 1: length(exp.expID(idxExp).unit)
        if strcmp(exp.expID(idxExp).unit(idxUnit).cellType, 'RS')
        for idxEvent = 1: length(exp.expID(idxExp).unit(idxUnit).event)-1
             if  exp.expID(idxExp).unit(idxUnit).event(idxEvent).excitatoryResponse == 1
                l = l + 1;
                vect_B1(l,:) = exp.expID(idxExp).unit(idxUnit).event(idxEvent).spikeResponse;
             end
        end
        end
    end
end

fanoF_B1 = var(vect_B1,[],2)./mean(vect_B1,2);
fastTr_B1 = mean(vect_B1(:,1:5),2);
lastTr_B1 = mean(vect_B1(:,6:10),2);

figure( 'Position', [10 10 1400 800])
d{1} = fanoF_A1; d{2} = fanoF_A5; d{3} = fanoF_B1; 
h1 = raincloud_plot(d{1},'box_on', 1, 'color', cb_d1, 'alpha', 0.6,...
    'box_dodge', 1, 'box_dodge_amount', .1, 'dot_dodge_amount', .3,...
    'box_col_match', 1);
h2 = raincloud_plot(d{2}, 'box_on', 1, 'color', cb_d5, 'alpha', 0.6,...
    'box_dodge', 1, 'box_dodge_amount', .5, 'dot_dodge_amount', .7, 'box_col_match', 1);
h3 = raincloud_plot(d{3}, 'box_on', 1, 'color', cb_d6, 'alpha', 0.6,...
    'box_dodge', 1, 'box_dodge_amount', .9, 'dot_dodge_amount', 1.1, 'box_col_match', 1);
h1{1, 1}.LineStyle = 'none'; h2{1, 1}.LineStyle = 'none';h3{1, 1}.LineStyle = 'none';
ylim([-0.18 0.2]); xlim([0 70]);title('Fano Factor Distributions for Excitatory Neurons');
ylabel('PDF'); xlabel('Fano Factor');
legend([h1{1},h2{1},h3{1}],'Day 1 - Odor Set A','Day 5 - Odor Set A','Day 1 - Odor Set B'); yticks([0:0.05:0.2]);
legend boxoff 
set(gcf,'color','white', 'PaperPositionMode', 'auto'); set(gca, 'box', 'off', 'tickDir', 'out');

[H,P]=kstest2( fanoF_A1,fanoF_A5);
[H,P]=kstest2( fanoF_A5,fanoF_B1)
[H,P]=kstest2( fanoF_A1,fanoF_B1)

%% Decoding

load('A1_-1_4_TimeCourse_200by100ms')
f8 = figure( 'Position', [10 10 1400 800])
numOfStimuli = 7;
tcperf_mu_A1 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_A1 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
k = shadedErrorBar(1:49, tcperf_mu_A1, tcperf_stdev_A1,'transparent',1);
hold on
set(k.edge,'LineWidth',.5,'LineStyle','none');
k.mainLine.LineWidth = 2; k.mainLine.Color= cb_d1;
k.patch.FaceAlpha = .3; k.patch.FaceColor= cb_d1;
load('A5_-1_4_TimeCourse_200by100ms')
tcperf_mu_A5 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_A5 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
z = shadedErrorBar(1:49, tcperf_mu_A5, tcperf_stdev_A5,'transparent',1);
hold on
set(z.edge,'LineWidth',.5,'LineStyle','none');
z.mainLine.LineWidth = 2; z.mainLine.Color= cb_d5;
z.patch.FaceAlpha = .3; z.patch.FaceColor= cb_d5;
load('B1_-1_4_TimeCourse_200by100ms')
tcperf_mu_PRh1 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_PRh1 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
z3 = shadedErrorBar(1:49, tcperf_mu_PRh1, tcperf_stdev_PRh1,'transparent',1);
hold on
set(z3.edge,'LineWidth',.5,'LineStyle','none');
z3.mainLine.LineWidth = 2; z3.mainLine.Color= cb_d6;
z3.patch.FaceAlpha = .3; z3.patch.FaceColor= cb_d6;

yline(1/7, '--', {'Chance Level'});
M = .015*ones(1,20);
ax = plot(10:29, M, 'k', 'LineWidth', 1.5);
ax.Clipping = 'off';
hold on
ylabel('Accuracy \%')
title('Decoding Over Timecourse - 200ms Bin $/$ 100ms Step')
ylim([0 1])
xlim([1 49]); xticks([]); 
ax1 = gca
ax1.XAxis.Visible = 'off';

set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
legend([k.mainLine,z.mainLine,z3.mainLine] ,{'Day 1 - Odor Set A','Day 5 - Odor Set A','Day 1 - Odor Set B'},'Location','northwest')
legend boxoff
%% TOP 14

%%

f8 = figure( 'Position', [10 10 1400 800])
numOfStimuli = 7;
tiledlayout(2,5)
nexttile([2,2])

BARTOP14 = bar(data,'stacked')
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')
BARTOP14.FaceColor = cb_d5
BARTOP14.LineStyle = 'none'
title('Most Informative Neurons (p $<$ .01) - Mice Representation')
ylabel('Percentage of Neurons')
xlabel('Mice ID')
nexttile([2,3])
%%
f8 = figure( 'Position', [10 10 1400 800])
load('A5_-1_4_TimeCourse_200by100ms_ONLY33')
tcperf_mu_A5 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_A5 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
z = shadedErrorBar(1:49, tcperf_mu_A5, tcperf_stdev_A5,'transparent',1);
hold on
set(z.edge,'LineWidth',.5,'LineStyle','none');
z.mainLine.LineWidth = 2; z.mainLine.Color= cb_d5;
z.patch.FaceAlpha = .3; z.patch.FaceColor= cb_d5;

load('A5_-1_4_TimeCourse_200by100ms_ONLY33woM1')
tcperf_mu_PRh1 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results];
tcperf_stdev_PRh1 = [DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.stdev.over_resamples];
z1 = shadedErrorBar(1:49, tcperf_mu_PRh1, tcperf_stdev_PRh1,'transparent',1);
hold on
set(z1.edge,'LineWidth',.5,'LineStyle','none');
z1.mainLine.LineWidth = 2; z1.mainLine.Color= cb_d5; z1.mainLine.LineStyle= '--';
z1.patch.FaceAlpha = .3; z1.patch.FaceColor= cb_d5;

yline(1/7, '--', {'Chance Level'});
M = .015*ones(1,20);
ax = plot(10:29, M, 'k', 'LineWidth', 1.5);
ax.Clipping = 'off';
hold on
ylabel('Accuracy \%')
title('Decoding Over Timecourse - 200ms Bin $/$ 100ms Step')
ylim([0 1])
xlim([1 49]); xticks([]); 
ax1 = gca
ax1.XAxis.Visible = 'off';

legend([z.mainLine,z1.mainLine] ,{'A5: with Mice \#1','A5: without Mice \#1'},'Location','northwest')
legend boxoff
set(gcf,'color','white', 'PaperPositionMode', 'auto');
set(gca, 'box', 'off', 'tickDir', 'out')