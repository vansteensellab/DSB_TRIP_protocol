


awk '{print $4, $2}' $1 | \
	sort -k1 | \
	uniq -c | \
	sort -k1rn | \
	awk -vOFS='\t' '{print $1, $2, $3}'
