library(data.table)
library(Biostrings)

save(snakemake, file='test.Rdata')



min_count = snakemake@params[['min_count']]
breaksite = snakemake@params[['breaksite']]
sequence = snakemake@params[['sequence']]

command = paste("cat ",  paste(snakemake@input, collapse = ' '),
                "| awk '{if ($3==\"ins\"||$3==\"del\"){print $0}}'")
print(command)
dt = fread(cmd=command, col.names=c('count', 'seq', 'call'), key='seq')

sum_dt = dt[,list(call=call[1], count=sum(count)),by=c('seq')]


test <- function(start, end, width, aln, breaksite){
    left = T
    right = T
    i = 0
    res = c(start)
    while ((left | right) & i < width){
        i = i + 1
        left_start = start-i
        left_end = end-i+1
        left_pat = substr(pattern(aln), left_start, start-1)
        left_sub = substr(subject(aln), left_end, end)

        right_end = end+i
        right_start = start+i-1
        right_pat = substr(pattern(aln), end+1, right_end)
        right_sub = substr(subject(aln), start, right_start)

        if (left_pat==left_sub){
            res = rbind(res, left_start)
        } else {
            left=F
        }

        if (right_pat==right_sub){
            res = rbind(res, right_start)
        } else {
            right=F
        }
    }

    return(paste(range(res) - breaksite - 1, collapse = ':'))
}


align <- function(wt_seq, mut_seq, type, breaksite){
    aln = pairwiseAlignment(mut_seq, wt_seq)
    if (type == "ins"){
        mut_table = as.data.table(insertion(aln)[[1]])
        if (nrow(mut_table) > 0){
            mut_table[, region:="NA"]
        } else {
            type = 'not_clear'
            mut_table[, region:=character()]
        }
    } else if (type == 'del'){
        mut_table = as.data.table(deletion(aln)[[1]])
        if(nrow(mut_table) > 0){
            mut_table[, region:=test(start,end,width,aln, breaksite),
                      by=c('start','end','width')]
        } else {
            mut_table[, region:=character()]
        }
    }
    else {
        print(aln)
    }
    mut_table[,nuc:=ifelse(type=="ins", substr(mut_seq,start,end) ,"")]
    mut_table[,start:=start - breaksite - 1]
    mut_table[,end:=end - breaksite - 1]
    if (nrow(mut_table) == 1){
        return(as.list(mut_table))
    } else if (nrow(mut_table) == 0){
        return(as.list(mut_table))
    } else {
        return(as.list(mut_table[which.min(abs((start+end)/2)),]))
    }
}


mut.dt = sum_dt[count > min_count,
                align(sequence, seq, call, breaksite), by=c('seq', 'call', 'count')]

fwrite(mut.dt, snakemake@output[[1]], sep='\t')
