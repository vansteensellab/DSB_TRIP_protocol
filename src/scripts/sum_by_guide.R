library(data.table)


save(snakemake, file='test.Rdata')

command = paste("cat ",  paste(snakemake@input, collapse = ' '),
                "| awk '{if ($3==\"ins\"||$3==\"del\"){print $0}}'")
dt = fread(cmd=command, col.names=c('count', 'seq', 'call'), key='seq')

sum_dt = dt[,list(count=sum(count)),by=c('seq', 'call')]
fwrite(sum_dt[order(count, decreasing=T),], file=snakemake@output[[1]], sep='\t')
