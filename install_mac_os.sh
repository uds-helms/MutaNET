#!/bin/bash
echo
echo "Burrows-Wheeler Aligner, SAMTools and VarScan are required for the NGS pipeline of MutaNET. Do you wish to install them on your system? This will also install the package manager Homebrew if it is not installed yet."
echo
read -p "Type y for yes or n for no." -n 1
echo
echo

if [[ $REPLY =~ ^[Yy]$ ]]
then 
	which -s brew
	if [[ $? != 0 ]] ; then
		# Install Homebrew
		ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
		echo
	else
		brew update
		echo
	fi
	brew cask install java
	echo
	brew install bwa
	echo
	brew install samtools
	echo
	brew install brewsci/bio/varscan
	echo
fi

echo
echo "Python 3 is required if you want to run MutaNET from source code. It is NOT required when using the executable. Do you wish to install Python 3 and the required Python 3 packages (numpy, matplotlib, scipy, fpdf, pyyaml) on your system? This will also install the package manager Homebrew if it is not installed yet."
echo
echo
read -p  "Type y for yes or n for no." -n 1
echo
echo

if [[ $REPLY =~ ^[Yy]$ ]]
then 
	which -s brew
	if [[ $? != 0 ]] ; then
		ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
		echo
	else
		brew update
		echo
	fi
	brew install python3
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
fi
