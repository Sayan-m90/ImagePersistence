for f in 256_ObjectCategories/*
do

for f1 in $f/*.jpg
do
echo "1. Keep in a folder with Caltech-256 Dataset."
echo "2. Run using ./caltech.sh "
echo "3. Outputs stat.csv containing topological signature for each image."
f3=${f1:21:3}
f2="stat.csv"
echo $f2
{ ./pers -i $f1 -o $f2 -n 17 -s 15 -v 10 -b $f3; }
done
done