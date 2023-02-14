% Normalize gene expression across ROIs
% Return normalized expression of target gene or gene set (specified by
% gene_names)

function gene_expression = normalizeExpression(exp_all_gene, gene_names)
    
    % normalize across all ROIs for current subject:
    mu = mean(exp_all_gene,2);
    sdev = std(exp_all_gene,0,2);
       
    norm_exp_all = (exp_all_gene - mu)./sdev;
    
    % extract unit normalized expression of specific gene set:
    gene_expression = zeros(length(gene_names),size(exp_all_gene,2));
    for i = 1:1:length(gene_names)
        gene_expression(i,:) = norm_exp_all(i,:);
    end
end

