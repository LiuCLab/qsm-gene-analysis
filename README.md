Quantitative susceptibility mapping in the brain reflects spatial expression of genes involved in iron homeostasis and myelination

QSM pre-processing:

zeropad.sh - Zeropad subject QSMs so they are same size as group atlas

inv.sh - Transform segmentation from atlas to subject space using registration

reg.sh - Perform registration between subject QSMs and QSM atlas


geneQSMComparison.m - Script to run gene processing functions, regression function, and PLSR. Calls the following functions:

qsmAvgSubj.m - Load QSM subject MRIs and registered segmentations, and use avgqsm function to calculate Avg. QSM across all ROIs for each subject

avgqsm.m - Calculate average QSM across an ROI, given mask and image

getAllExpressionDG.m - Processing of gene expression data, including matching regions between segmentations, normalizing across ROIs, and averaging across subjects

matchSamplesDG.m - Match samples between segmentations

normalizeExpression.m - Normalize expression across genes

regressionGeneDG.m - Perform regression between gene expression and QSM

Other files"

SampleAnnot1.xlsx - SampleAnnot6.xlsx - Lists tissue samples for Allen Human Brain Atlas (AHBA) subjects 1-6

ROI_labels.xlsx - Match table, used to match AHBA tissue samples to regions in QSM segmentation
