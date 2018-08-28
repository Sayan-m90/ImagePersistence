echo "Extension of image files: "
read ext
echo ".csv file to push topological signatures: "
read output
for f2 in input_images/*.$ext
do
	
{ ./pers -i $f2 -o $output -n 15 -s 12 -v 10; }
done
