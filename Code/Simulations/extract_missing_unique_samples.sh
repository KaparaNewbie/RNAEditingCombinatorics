# #!/usr/bin/env bash
# set -euo pipefail

# DIR="${1:-.}"

# tmp="$(mktemp)"
# trap 'rm -f "$tmp"' EXIT

# # Extract per-file value: col10 from the row with max col5
# # find "$DIR" -type f -name '*.DistinctUniqueProteins.*.Updated.csv' -print0 \

# find "$DIR" -type f -name '*.DistinctUniqueProteins.*.Updated.*csv' -print0 \
# | while IFS= read -r -d '' f; do
#     awk -F',' 'NR==1{next} { if ($5+0 > max) { max=$5+0; val=$10+0 } } END { print val }' "$f"
#   done > "$tmp"

# echo "N_files: $(wc -l < "$tmp")"

# # Summary stats (numeric)
# awk '
# function abs(x){return x<0?-x:x}
# {
#   x=$1+0
#   n++
#   sum+=x
#   sumsq+=x*x
#   if(n==1 || x<min) min=x
#   if(n==1 || x>max) max=x
#   a[n]=x
# }
# END{
#   if(n==0){print "No matching files / no data"; exit 0}
#   mean=sum/n
#   var=(sumsq/n)-(mean*mean)
#   if(var<0) var=0
#   sd=sqrt(var)

#   # median (sort array)
#   asort(a)
#   if(n%2==1) med=a[(n+1)/2]
#   else med=(a[n/2]+a[n/2+1])/2

#   printf "min\t%g\nmax\t%g\nmean\t%g\nmedian\t%g\nsd\t%g\n", min, max, mean, med, sd
# }' "$tmp"





# new version


#!/usr/bin/env bash
set -euo pipefail

DIR="${1:-.}"

tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

find "$DIR" -type f -name '*.DistinctUniqueProteins.*.Updated.*csv' -print0 \
| while IFS= read -r -d '' f; do
    awk '
      NR==1{
        # autodetect delimiter by seeing which split yields >=10 fields on the header
        n_tab = split($0, a, "\t")
        n_com = split($0, b, ",")
        FS = (n_tab >= 10 ? "\t" : ",")
        $0 = $0   # re-split current record using chosen FS
        next      # skip header
      }
      {
        # pick row with maximum col5, store col10 from that row
        if (($5+0) > max) { max = $5+0; val = $10+0 }
      }
      END{ print val+0 }
    ' "$f"
  done > "$tmp"

echo "N_files_total: $(wc -l < "$tmp")"

awk '
{
  x=$1+0

  # ALL
  n++; sum+=x; sumsq+=x*x
  if(n==1 || x<min) min=x
  if(n==1 || x>max) max=x
  a[n]=x

  # SUBSET >=1
  if (x >= 1) { m++; s+=x; ss+=x*x; b[m]=x }
}
END{
  if(n==0){ print "No matching files / no data"; exit 0 }

  # ALL
  mean=sum/n
  var=(sumsq/n)-(mean*mean); if(var<0) var=0
  sd=sqrt(var)
  asort(a)
  med = (n%2 ? a[(n+1)/2] : (a[n/2]+a[n/2+1])/2)

  print "=== ALL FILES ==="
  printf "min\t%g\nmax\t%g\nmean\t%g\nmedian\t%g\nsd\t%g\n\n", min, max, mean, med, sd

  # >=1
  print "=== FILES WITH NumMissingUniqueSamples >= 1 ==="
  printf "count\t%d\n", m
  if(m>0){
    mean2=s/m
    var2=(ss/m)-(mean2*mean2); if(var2<0) var2=0
    sd2=sqrt(var2)
    asort(b)
    med2 = (m%2 ? b[(m+1)/2] : (b[m/2]+b[m/2+1])/2)
    printf "min\t%g\nmax\t%g\nmean\t%g\nmedian\t%g\nsd\t%g\n", b[1], b[m], mean2, med2, sd2
  }
}' "$tmp"