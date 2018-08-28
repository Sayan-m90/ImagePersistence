for f in 256_ObjectCategories/*
do

for f1 in $f/*.jpg
do
echo "Program to take as input Caltech256 and output topological signature."
f3=${f1:21:3}
f2="stat.csv"
echo $f2
{ ./pers -i $f1 -o $f2 -n 17 -s 15 -v 10 -b $f3; }
done
done