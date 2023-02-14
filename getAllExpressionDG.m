% Match AHBA tissue samples to ROIs in QSM segmentation, normalize gene
% expression across ROIs, and average across AHBA subjects
% all_genes = list of gene names
% expressionAll = gene expression matrix for all subjects, post gene-probe
% reannotation, probe filtering, and probe selection
% qsm_45 = average QSM across ROIs, median age 45
% Returns:
% avg_all_subs = normalized gene expression across ROIs in QSM
% segmentation, averaged across AHBA subjects
% qsm_45_all = average QSM corresponding to matched ROIs, by subject
% exp_all = normalized gene expression across ROIs in QSM segmentation, by
% subject
% all_rois = all matched ROIs, after combining AHBA subjects (1:1:34 
% if all deep grey nuclei ROIs are matched)

function [all_rois,avg_all_subs,qsm_45_all,exp_all] = getAllExpressionDG(all_genes,expressionAll,qsm_45)

    exp_all = cell(6,1);
    match_idx_all = cell(6,1);
    qsm_45_all = cell(6,1);
    
    
    % Step through AHBA subjects:
    for i = 1:1:6
        subj = i; % use Subject i only
        % Match AHBA tissue samples to QSM segmentation
        [me_ROI, match_idx] = matchSamplesDG('ROI_labels.xlsx', "SampleAnnot" + subj + ".xlsx", expressionAll{subj}');
        % Normalize matched gene expression across ROIs
        gene_expression = normalizeExpression(me_ROI, all_genes);
    
        exp_all{i} = gene_expression;
        match_idx_all{i} = match_idx';
        qsm_45_all{i} = qsm_45(match_idx);
        
    end

    % Combine gene expression across subjects
    % Average across subjects, for a given ROI
    a = match_idx_all';
    % Find all matched ROIs, across all subjects
    all_rois = unique([a{:}],'first');
    % Step through all unique regions, average multiple samples, if necessary
    avg_all_subs = zeros(length(all_genes),length(all_rois));
    for i = 1:1:length(all_rois)
        num = 0;
        sum = zeros(length(all_genes),1);
        for j = 1:1:6
            exp_subj = exp_all{j};
            subj_match = match_idx_all{j};
            if find(subj_match == all_rois(i))
                sum = sum + exp_subj(:,find(subj_match == all_rois(i)));
                num = num + 1;
            end
        end
        avg_all_subs(:,i) = sum./num;
    end
    
end