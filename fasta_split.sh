
INPUT=$1
OUTPUT=$2

echo "Started `date`"

echo "Creating output directory"
mkdir -p $OUTPUT

echo "Splitting Fasta..."
while read line ; do
  if [ ${line:0:1} == ">" ] ; then
    filename=$(echo "$line" | cut -d ":" -f1 | tr -d ">")
    touch "$OUTPUT"/"$filename".fasta
    echo "$line" >> "$OUTPUT"/"${filename}".fasta
  else
    echo "$line" >> "$OUTPUT"/"${filename}".fasta
  fi
done < $INPUT

echo "Completed `date`"
