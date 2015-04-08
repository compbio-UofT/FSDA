cat | awk -v begin=$1 -v end=$2 'BEGIN {sum=0} {split($2,a,":"); split(a[2],reg,"-"); if (!(reg[1]>end || reg[2]<begin)) {q=reg[1]; p=reg[2]; if (begin>q) q=begin;if (end<p) p=end; if ($1=="deletion") sum-=p-q; if ($1=="duplication") sum+=p-q;}} END {print sum}'

