#!/usr/bin/env bash
echo
echo "Burrows-Wheeler Aligner, SAMTools and VarScan are required for the NGS pipeline of MutaNET. Do you wish to install them on your system?"
echo
read -p "Type y for yes or n for no." -n 1
echo
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then 
	sudo apt-get install -y bwa
	echo
	sudo apt-get install -y samtools
	echo
	sudo apt-get install -y varscan
	echo
	sudo apt-get install -y default-jre
fi

echo
echo "Python 3 is required if you want to run MutaNET from source code. It is NOT required when using the executable. Do you wish to install Python 3 and the required Python 3 packages (numpy, matplotlib, scipy, fpdf, pyyaml) on your system?"
echo
echo
read -p  "Type y for yes or n for no." -n 1
echo
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then 
	sudo apt-get install -y python3
	echo
	sudo apt-get install -y python3-tk
	echo
	sudo apt-get install -y python3-pip
	echo
	pip3 install --upgrade pip --user
	echo
	pip3 install numpy --user
	echo
	pip3 install matplotlib --user
	echo
	pip3 install scipy --user
	echo
	pip3 install fpdf --user
	echo
	pip3 install pyyaml --user
	echo
fi
