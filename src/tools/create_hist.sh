frag_thresh=$1
if [[ -z "$frag_thresh" ]]
then
	frag_thresh=500
fi
awk -v frag_thresh=$frag_thresh '{if ($1>0 && $1<frag_thresh) print $1}' | awk -v frag_thresh=$frag_thresh 'BEGIN{for (i=1; i<=frag_thresh; i++) c[i]=0} {c[$1]+=1} END {for (i=1; i<=frag_thresh; i++) print c[i]}' | awk '{printf "%f ",$1} END { printf "\n"}'
#awk '{tmp=$9; if (tmp<0) tmp=tmp*(-1); print tmp}' | awk -v frag_thresh=$frag_thresh '{if ($1>0 && $1<frag_thresh) print $1}' | awk -v frag_thresh=$frag_thresh 'BEGIN{for (i=1; i<=frag_thresh; i++) c[i]=0} {c[$1]+=1} END {for (i=1; i<=frag_thresh; i++) print c[i]}' | awk '{printf "%f ",$1} END { printf "\n"}'
