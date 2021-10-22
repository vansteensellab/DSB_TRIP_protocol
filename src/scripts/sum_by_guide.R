library(data.table)
## This function is used to create a list of unique sequences found in all
## experiments using a specific guide RNA to cut.
## Main purpose of this clustering is to reduce computational time for
## pairwise alignment or mutational analysis using inDelphi/SelfTarget.

save(snakemake, file='test.Rdata')

command = paste("cat ",  paste(snakemake@input, collapse = ' '),
                "| awk '{if ($3==\"ins\"||$3==\"del\"){print $0}}'")
dt = fread(cmd=command, col.names=c('count', 'seq', 'call'), key='seq')

sum_dt = dt[,list(count=sum(count)),by=c('seq', 'call')]
fwrite(sum_dt[order(count, decreasing=T),], file=snakemake@output[[1]], sep='\t')
