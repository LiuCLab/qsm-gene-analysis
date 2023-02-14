% Function for performing regression between specified gene set and average
% QSM across deep grey nuclei regions
% all_rois = all ROIs with matching AHBA tissue samples (all deep grey
% nuclei regions in QSM segmentation have a match)
% qsm_45 = average QSM across deep grey nuclei ROIs
% gene_idx = indices of genes in set
% gene_name = names of genes in set
% avg_all_subs = normalized gene expression in deep grey nuclei ROIs,
% averaged across AHBA subjects
% qsm_45_all = average QSM across regions matched to AHBA samples, by subject
% exp_all = normalized gene expression in deep grey nuclei ROIs, by subject
% plotMode = string to specify plotting of regression results (options:
% "" = no plot, "Avg" = averaged across subjects, "Subject" = by subject
% Returns:
% p-vals, slope, mdl = p-values, slope, and regression model between gene
% expression and QSM
% plot of regression results, if specified

function [p_vals, slope,mdl] = regressionGeneDG(all_rois,qsm_45,gene_idx,gene_name,avg_all_subs,qsm_45_all,exp_all, plotMode)

    avg_all_gene = avg_all_subs(gene_idx,:);
    
    p_vals = zeros(size(gene_idx));
    slope = zeros(size(gene_idx));
    intcep = zeros(size(gene_idx));
    
    qsm_45_matched = qsm_45(all_rois);
    mu = mean(qsm_45_matched);
    sdev = std(qsm_45_matched);
    % Normalize QSM across ROIs
    qsm_norm = (qsm_45_matched - mu*ones(size(qsm_45_matched)))./sdev;
    
    mdl_pltx = zeros(length(gene_idx),length(qsm_norm));
    mdl_plty = zeros(length(gene_idx),length(qsm_norm));
    mdl_pltpred = zeros(length(gene_idx),length(qsm_norm));
    
    % Regression between QSM and gene expression, averaged across subjects
    for i = 1:length(gene_idx)
        % Fit regression model for gene i with normalized QSM
        mdl = fitlm(avg_all_gene(i,:),qsm_norm);
        mdl_pltx(i,:) = mdl.Variables.x1;
        mdl_plty(i,:) = mdl.Variables.y;
        mdl_pltpred(i,:) = mdl.Fitted;
        slope(i) = mdl.Coefficients.Estimate(2);
        p_vals(i) = mdl.Coefficients.pValue(2);
        intcep(i) = round(mdl.Coefficients.Estimate(1),2);
    end
    
    % Plot regression results, by AHBA subject
    if strcmp(plotMode, "Subject")
        p_vals_subj = zeros(length(gene_idx),6);
        slope_subj = zeros(length(gene_idx),6);
        cept_subj = zeros(length(gene_idx),6);

        colors = [1 0 0; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.4940 0.1840 0.5560; 1 0 1];
        for j = 1:1:6
            subj = j;
            expr_sub = exp_all{subj};
            expr_sub_gene = expr_sub(gene_idx,:);
            
            qsmsub = qsm_45_all{subj};
            mu = mean(qsmsub);
            sdev = std(qsmsub);
            qsm_norm = (qsmsub - mu*ones(size(qsmsub)))./sdev;

            for i = 1:1:length(gene_idx)
                mdl_sub = fitlm(expr_sub_gene(i,:),qsm_norm);
                p_vals_subj(i,j) = mdl_sub.Coefficients.pValue(2);
                slope_subj(i,j) = round(mdl_sub.Coefficients.Estimate(2),2);
                cept_subj(i,j) = round(mdl_sub.Coefficients.Estimate(1),2);
                scatter(mdl_sub.Variables.x1,mdl_sub.Variables.y,'x','MarkerFaceColor',colors(j,:),'MarkerEdgeColor',colors(j,:))
                hold on
                plot(mdl_sub.Variables.x1,mdl_sub.Fitted,'Color',colors(j,:))
                hold on
            end  
        end
        
        % Save regression results by subject for each gene of interest
        subj_table = table(p_vals_subj, slope_subj, cept_subj);
        sub_round = [round(subj_table.slope_subj,2)',subj_table.p_vals_subj'];
        save('Subject_regressions ' + gene_name,'subj_table')
        save('Subject_table_' + gene_name,'sub_round')
        
        leg = "";
        for i = 1:1:6
            leg(i) = "Slope = " + slope_subj(1,i) + ", p = " + p_vals_subj(1,i);
        end
        
        legend('Subject 1',leg(1),'Subject 2',leg(2),'Subject 3',leg(3),'Subject 4',leg(4),'Subject 5', leg(5),'Subject 6', leg(6))
        title(gene_name + " Expression vs. Avg. QSM for All Subjects")
        xlabel('Gene expression (z-score)')
        ylabel('Avg. QSM across ROIs (z-score)')
    
    % Plot regression results, gene expression averaged across subjects    
    elseif strcmp(plotMode, "Avg")
        for i = 1:1:length(gene_idx)
            plot(mdl_pltx(i,:),mdl_plty(i,:),'x')
            hold on
            plot(mdl_pltx(i,:),mdl_pltpred(i,:))
          
            xlim([min(mdl_pltx(i,:)),max(mdl_pltx(i,:))])
            legend("Averaged across subjects","Prediction, p = " + mdl.Coefficients.pValue(2),'FontSize',18)
            title({gene_name + " Expression vs. Avg. QSM"; "y = " + round(slope(i),2) + "*x +"+intcep(i)},'FontSize',24 )
            xlabel('Gene expression (z-score)','FontSize',20)
            ylabel("Avg. QSM across ROIs (z-score)",'FontSize',20)
        end

    end
end 