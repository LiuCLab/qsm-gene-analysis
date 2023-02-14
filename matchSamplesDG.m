% Match AHBA tissue samples, listed in SampleAnnot (by subject), to ROI's
% in QSM segmentation
% Match according to matchTable, specified by matchTableName =ROI_labels.xlsx
% expressionSubj = matrix of gene expression across AHBA samples
% A sample located at column j in expressionSubj corresponds to the region
% in row j of SampleAnnot (specified by SamplesName)
% Return:
% matched_exp = gene expression across ROIs specified in QSM segmentation
% match_idx = indices of matched ROIs
function [matched_exp, match_idx] = matchSamplesDG(matchTableName, SamplesName, expressionSubj)

    % Import Match Table from file specified by matchTableName
    opts = detectImportOptions(matchTableName);
    opts.VariableTypes{2} = 'string';
    opts.VariableTypes{3} = 'string';
    opts.VariableTypes{4} = 'string';
    opts.VariableTypes{5} = 'string';
    opts.VariableTypes{6} = 'string';
    opts.VariableTypes{7} = 'string';
    opts.VariableTypes{8} = 'string';
    opts.DataRange = 'A1';
    match = readtable(matchTableName,opts);
    matchTable = match(:,2:8);
    
    % Read SampleAnnot matrix (specified by SamplesName)
    samples_names1 = readmatrix(SamplesName,'OutputType','string');
    samples_names = samples_names1(:,6);

    me_ROI = zeros(size(expressionSubj,1),size(matchTable,1));
    matches = table2array(matchTable);
    for i = 1:1:size(matchTable,1)
        part_means = zeros(size(expressionSubj,1),6);
        for j = 1:1:6
            if not(ismissing(matches(i,j)))
                % Find AHBA samples corresponding to those listed in
                % matchTable
                idx = strfind(samples_names,matches(i,j));
                rowsSA = find(~cellfun(@isempty,idx));
                if ~isempty(rowsSA)
                    % Find gene expression values from matched samples
                    colsME = rowsSA;
                    me_rois = expressionSubj(:,colsME);
                    % Average if multiple AHBA samples from same region
                    part_means(:,j) = mean(me_rois,2);
                end
            end
        end
        % Average AHBA samples corresponding to same ROI
        % Only include in mean calculation if column is nonzero:
        num_cols = nnz(sum(part_means,1)); % find # of nonzero cols
        avg_ROI = sum(part_means,2)/num_cols; % add cols together and divide by number of cols summed
        me_ROI(:,i) = avg_ROI;
    
    end
    
    id = 1:1:size(matchTable,1);
    rois = id;
    % Take out ROIs with no samples:
    no_match = isnan(me_ROI(1,:));
    rois(no_match) = [];

    % Find indices of matched ROIs in id:
    match_idx = zeros(length(rois),1);
    for i = 1:1:length(rois)
        match_idx(i) = find(id == rois(i));
    end
    
    matched_exp = me_ROI(:,~no_match);
    
    
end