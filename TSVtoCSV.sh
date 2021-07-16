# tsv to csv
sed -i -e 's/,/;/g' $1 # replace commas in drug names with semicolons
sed -i -e 's/   /,/g' $1 # replace tabs between columns with commas