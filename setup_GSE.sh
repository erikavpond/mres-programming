!/bin/bash
mkdir erika_module4 #this will be the working directory

mv master_script.sbatch erika_module4/ #moves both other downloaded scripts to correct working directory
mv GSE_analysis.R erika_module4/

cd erika_module4 #working directory

#downloading data
set -euo pipefail
mkdir -p GSE163974 results
# Above is creating the directory structure

echo "Downloading data..."
wget -O GSE163974/GSE163974_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE163973&format=file"

echo "Extracting..."
##tar -xvf data/GSE282885_RAW.tar -C data/

echo "Done. File for analysis located in GSE163974/  directory"
tar -xvf GSE163974/GSE163974_RAW.tar -C GSE163974/
cd GSE163974/ 
#creates a directory (unzip, untar) per sample, containing barcodes.tsv/features.tsv/matrix.mtx files.
for f in *.tar.gz; do
	tar -xvzf "$f"
done
cd ..

#now should be in erika_module4 directory, ready to run R analysis!
#sbatch script is next to run, which submits the large analysis
#please see master_script.sbatch