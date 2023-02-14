% Calculate average QSM across region corresponding to id in segmentation
% seg
% qsm is the QSM in the same space as seg
% Returns:
% avg_qsm = average QSM across specified region
function avg_qsm = avgqsm(seg, id, qsm)
    mask = (seg == id); 
    dmask = double(mask);
    dmask(dmask == 0) = NaN;
    roi = qsm.*dmask;
    avg_qsm = mean(roi(~isnan(roi)));
end