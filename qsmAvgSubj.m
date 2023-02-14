% Calculate average QSM across deep grey nuclei ROIs, for all QSM subjects
% Subject QSMs are located in mydir1, and their corresponding segmentations
% (registered to subject space) are located in mydir2
% First column of ROItable contains the labels corresponding to deep grey
% nuclei ROIs in the segmentation
% Returns:
% avg_qsm = average QSM across 34 deep grey nuclei regions, for all 9
% subjects

function avg_qsm = qsmAvgSubj(ROItable, mydir1, mydir2)
    % Get segmentation labels of ROIs
    id = readmatrix(ROItable,'Range',[1 1 34 1], 'OutputType', 'double');
    
    qsm_subs = dir(fullfile(mydir1,'*padded.nii'));
    seg_subs = dir(fullfile(mydir2,'*.nii.gz'));

    avg_qsm = zeros(length(qsm_subs),length(id));
    for k = 1:length(qsm_subs)
        % Get subject QSM
        qsm_name = qsm_subs(k).name;
        qsm_path = fullfile(mydir1, qsm_name);
        qsm = niftiread(qsm_path);
        % Get segmentation registered to QSM subject
        seg_name = seg_subs(k).name;
        seg_path = fullfile(mydir2, seg_name);
        seg = niftiread(seg_path);
        fprintf(1, 'Now reading %s\n', qsm_name);
        fprintf(1, 'Now reading %s\n', seg_name);
        for i = 1:1:length(id)
            % Calculate average QSM across deep grey nuclei ROIs
            avg_qsm(k,i) = avgqsm(seg,id(i),qsm);
        end
    end
end
