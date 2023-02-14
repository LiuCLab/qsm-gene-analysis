Gene processing and regression with QSM:
geneQSMComparison.m - Script to run gene processing functions, regression function, and PLSR
getAllExpressionDG.m - Processing of gene expression data, including matching regions between segmentations, normalizing across ROIs, and averaging across subjects
matchSamplesDG.m - Match samples between segmentations
normalizeExpression.m - Normalize expression across genes
regressionGeneDG.m - perform regression between gene expression and QSM

QSM processing:
zeropad.sh - Zeropad subject QSMs so they are same size as group atlas
inv.sh - Transform segmentation from atlas to subject space using registration
reg.sh - Perform registration between subject QSMs and QSM atlas
qsmAvgSubj.m - Load QSM subject MRIs and registered segmentations, and use avgqsm function to calculate Avg. QSM across all ROIs for each subject
avgqsm.m - Calculate average QSM across an ROI, given mask and image
