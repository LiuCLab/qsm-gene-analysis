%% Loading, processing QSM and gene expression

% Directory containing subject QSMs
mydir1 = 'C:\Users\Zoe\Documents\MATLAB\iron_project_mat\AHBA_gene\GeneAnalysis\qsm_subjects\LPI_BB';
% Directory containing segmentation registered to each subject
mydir2 = 'C:\Users\Zoe\Documents\MATLAB\iron_project_mat\AHBA_gene\GeneAnalysis\qsm_subjects\Seg_bysubj';

% Load Avg. QSM across subjects
avg_qsm_all = qsmAvgSubj('ROI_labels.xlsx', mydir1, mydir2); % Average QSM in deep grey nuclei for all Subjects, Median Age 45

% Average QSM across all subjects:
qsm_45 = mean(avg_qsm_all);

% Gene-Probe Reannotation, Probe Filtering, and Probe Selection (via Differential Stability)
% performed using software package from Arnatkevic?i?t? et al., 2019
% .mat file is output from using Arnatkevic?i?t? et al., 2019 software and
% contains the following:
% expressionAll = gene expression across all AHBA subjects
% options = input conditions chosen to produce .mat
% probeInformation = ID information for probes and genes
% sampleInfo = Locations of tissue samples for each AHBA subject
load('MicroarrayDataWITHcustProbesUpdatedXXXDSQC.mat')
all_genes1 = probeInformation.GeneSymbol;

all_genes = strings(size(all_genes1)); % Array of all Gene IDs
for i = 1:1:length(all_genes1)
    all_genes(i) = string(all_genes1{i,1}); % Convert to string
end

% Matches AHBA sampled regions to QSM segmentation, normalizes gene
% expression across ROIs, and averages across AHBA subjects
% requires SampleAnnot.xlsx for AHBA subjects 1-6, and match table (ROI_labels.xlsx) between
% ROIs in QSM atlas segmentation and AHBA sampled regions
[all_rois,avg_all_subs,qsm_45_all,exp_all] = getAllExpressionDG(all_genes,expressionAll,qsm_45);


%% Perform regression between all genes and the average QSM:
% all_rois = all matched ROIs (1:1:34 for all ROIs in deep grey nuclei)
% qsm_45 = average QSM, median age 45
% 1:1:length(all_genes) = indices of genes for regression
% all_genes = names of genes for regression
% avg_all_subs = expression of all genes, averaged across subjects
% qsm_45_all = average QSM across matched regions, by subject
% exp_all = expression of all genes across matched regions, by subject
% "" = indicates no plot

[all_p,slope,~] = regressionGeneDG(all_rois,qsm_45,1:1:length(all_genes),all_genes,avg_all_subs,qsm_45_all,exp_all,"");
all_plog = -log10(all_p);


%% Benjamini-Hochberg
% Perform B-H correction to linear regression p-values:

[B,I] = sort(all_p);
Q = .0001;
bh_cval = zeros(size(all_p));
for i = 1:1:length(all_p)
    bh_cval(i) = I(i)*Q/length(all_p);
end

passed = find(all_p < bh_cval);
[M,~] = max(all_p(passed));

% All p-values from regression with QSM, corrected by B-H: 
bh_sig_p = all_p(find(all_p < M));
disp("Maximum significant p-value, B-H corrected")
max(bh_sig_p)

% All genes with significant relationship with QSM, as corrected by B-H: 
bh_sig_genes = all_genes(find(all_p < M));

disp("False Discovery Rate")
.0001*length(bh_sig_genes)/length(all_p)

% Find which iron genes are significant, following B-H:
iron_genes = ["TF","TFRC","SLC40A1","FTH1","FTL","SLC11A2"];
bh_iron = zeros(size(iron_genes));
for i = 1:1:length(iron_genes)
    idxI = find(bh_sig_genes == iron_genes(i));
    if idxI > 0
        bh_iron(i) = idxI;
    end
end

just_idxI = bh_iron(find(bh_iron > 0));

disp("Iron genes significantly correlated with QSM, following B-H")
bh_sig_genes(just_idxI)

% Find which myelin genes are significant, following B-H:
myelin_genes = ["CNP","ILK","MAG","MAL","MBP","MOBP","MOG","OMG","CLDN11",...
        "PLP1","POU3F1","KLK6","EIF2AK3","GAL3ST1","OLIG2","PLLP","NRG1"];
bh_myelin = zeros(size(myelin_genes));
for i = 1:1:length(myelin_genes)
    idxM = find(bh_sig_genes == myelin_genes(i));
    if idxM > 0
        bh_myelin(i) = idxM;
    end
end

just_idxM = bh_myelin(find(bh_myelin > 0));

disp("Myelin genes significantly correlated with QSM, following B-H")
bh_sig_genes(just_idxM)

%% Plot regression results for specified gene set, AHBA subjects averaged

% Iron gene set
gene_set = ["TF","TFRC","SLC40A1","FTH1","FTL","SLC11A2"];

% Myelin gene set (significant following B-H):
%gene_set = ["CNP","MAG","MAL","MOBP","MOG","CLDN11","PLP1","GAL3ST1","PLLP"];

for i = 1:1:length(gene_set)  
    gene_name = gene_set(i);
    [all_p,slope,mdl] = regressionGeneDG(all_rois,qsm_45,find(all_genes == gene_name),gene_name,avg_all_subs,qsm_45_all,exp_all,"Subject");
    pause
    clf
end

%% Partial Least Squares regression
% Use outputs from getAllExpressionDG: all_rois, avg_all_subs

% options are Iron, Myelin, and IronMyelin
geneSet = 'IronMyelin';

if strcmp(geneSet, 'Iron')
    gene_set = ["TF","FTL","FTH1","SLC40A1"]';
elseif strcmp(geneSet,'Myelin')
    gene_set = ["CNP","MAG","MAL","MOBP","MOG","CLDN11","PLP1","GAL3ST1","PLLP"];
elseif strcmp(geneSet,'IronMyelin')
    gene_set = ["TF","FTL","FTH1","SLC40A1","CNP","MAG","MAL","MOBP","MOG","CLDN11","PLP1","GAL3ST1","PLLP"];
end

gene_idx = zeros(size(gene_set));
for i = 1:1:length(gene_set)
    gene_idx(i) = find(all_genes == gene_set(i));
end

gene_exp_interactions = avg_all_subs(gene_idx,:);

qsm_45_matched = qsm_45(all_rois);
mu = mean(qsm_45_matched);
sdev = std(qsm_45_matched);
qsm_region = (qsm_45_matched - mu*ones(size(qsm_45_matched)))./sdev;    

% PLSR between gene expression matrix and QSM
[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(gene_exp_interactions',qsm_region',size(gene_exp_interactions,1));

% Figures for model selection:

% Mean squared error vs. number of PLS components
figure
subplot(3,1,1)
yyaxis left
plot(0:size(gene_exp_interactions,1),MSE(1,:),'-o')
yyaxis right
plot(0:size(gene_exp_interactions,1),MSE(2,:),'-o')
legend('MSE Predictors','MSE Response')
xlabel('Number of Components')
title('MSE vs. Number of Components')

r_sq = zeros(size(gene_exp_interactions,1),1);
RSS = zeros(size(gene_exp_interactions,1),1);
PRESS = zeros(size(gene_exp_interactions,1),1);
SS = zeros(size(gene_exp_interactions,1),1);

% Cross validation for model selection:

for i = 1:1:size(gene_exp_interactions,1)
    c = cvpartition(34,'KFold',5);
    press = zeros(5,1);
    for j = 1:1:34
        test_set = j; % leave one out
        train_set = setdiff(1:1:34,j); 
        new_X = gene_exp_interactions(:,train_set)';
        new_y = qsm_region(train_set)';
        cv_X = gene_exp_interactions(:,test_set)';
        cv_y = qsm_region(test_set)';
        [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(new_X,new_y,i);
        yfit = [ones(size(cv_X,1),1) cv_X]*beta;
        press(j) = (cv_y - yfit).^2;
    end
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(gene_exp_interactions',qsm_region',i);
    yfit = [ones(size(gene_exp_interactions',1),1) gene_exp_interactions']*beta;
    SS(i) = sqrt(sum((qsm_region'-yfit).^2)/34); % RMSE
    PRESS(i) = sqrt(sum(press)/34); % RMSE cross validation
end

for i = 1:1:size(gene_exp_interactions,1)
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(gene_exp_interactions',qsm_region',i,'cv',5);
    yfit = [ones(size(gene_exp_interactions',1),1) gene_exp_interactions']*beta;
    TSS = sum((qsm_region'-mean(qsm_region)).^2);
    RSS(i) = sum((qsm_region'-yfit).^2);
    Rsquared = 1 - RSS(i)/TSS;
    r_sq(i) = Rsquared;
end

% R^2 vs. number of PLS components
subplot(2,1,2)
plot(1:size(gene_exp_interactions,1),r_sq,'-ko')
xlabel('Number of Components')
ylabel('R^{2}')
title('R^{2} vs. Number of Components')

% Root mean squared error and cross validation vs. number of PLS components
figure
plot(1:1:length(gene_set),PRESS)
hold on
plot(1:1:length(gene_set), SS)
legend('RMSE-CV','RMSE','FontSize',18)
title('Leave-one-out Cross-validation','FontSize',20)
xlabel('Number of PLSR Components','FontSize',18)
ylabel('RMSE (z-score)','FontSize',18)


%% PLSR with 2 components and confidence intervals

% Bootstrapping to calculate standard error on regression coefficients:
betas = zeros(length(qsm_region),size(gene_exp_interactions,1)+1);
for i = 1:1:length(qsm_region)
    test_set = i; % leave one out
    train_set = setdiff(1:1:34,i);
    new_X = gene_exp_interactions(:,train_set)';
    new_y = qsm_region(train_set)';
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(new_X,new_y,2);
    betas(i,:) = beta;
end

beta_avg = mean(betas,1);
% Calculate standard error from bootstrap results
betas_SE = sqrt(sum(((betas - beta_avg).^2),1))*length(qsm_region)/(length(qsm_region)-1);

% Calculate 95 percent confidence interval using standard error and
% Z-distribution
CI_neg = beta_avg - 1.96*betas_SE;
CI_pos = beta_avg + 1.96*betas_SE;

if strcmp(geneSet,'Iron')
    labels = {'y-intercept','TF','FTL','FTH1','SLC40A1'};
elseif strcmp(geneSet,'Myelin')
    labels = {'y-intercept','CNP','MAG','MAL','MOBP','MOG','CLDN11','PLP1','GAL3ST1','PLLP'};
elseif strcmp(geneSet,'IronMyelin')
    labels = {'y-intercept','TF','FTL','FTH1','SLC40A1','CNP','MAG','MAL','MOBP','MOG','CLDN11','PLP1','GAL3ST1','PLLP'};
end

% Bar plot PLSR regression coefficients and confidence intervals
bar(1:1:length(gene_set)+1,beta_avg,'g')
hold on
er = errorbar(1:1:length(gene_set)+1,beta_avg,CI_neg,CI_pos,'-o');
er.LineStyle = 'none';
er.Color = [0 0 0];
xticks(1:1:length(gene_set)+1)
xticklabels(labels)
xtickangle(45)
title('Coefficients on Genes in Relation to QSM','FontSize',24)
ylabel('Coefficient','FontSize',18)

%% Visualize Gene Expression x ROI matrix

imagesc(gene_exp_interactions)
title('Gene Expression x Region of Interest Matrix')
yticklabels({'TF','FTL','FTH1','SLC40A1','CNP','MAG','MAL','MOBP','MOG','CLDN11','PLP1','GAL3ST1','PLLP'});
yticks(1:1:13)
xticklabels(1:1:34)
xticks(1:1:34)
ylabel('Normalized Gene Expression (z-score)')
xlabel('Region of Interest (ROI)')
