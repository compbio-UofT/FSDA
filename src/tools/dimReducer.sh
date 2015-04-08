home="/home/aryan/fsda/ver-oct/"
tmp_dir="$home/tmp-dir"
ref_dir="$home/gen-ref"

cat > $tmp_dir/__hist_$$
pypy $home/src/tools/estimate.py $tmp_dir/__hist_$$ $ref_dir/newStats | tail -n 1
rm $tmp_dir/__hist_$$
