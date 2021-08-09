if [ -e v1.3.zip ]; then
	rm v1.3.zip
fi

wget https://github.com/FredHutch/SEACR/archive/refs/tags/v1.3.zip
unzip v1.3.zip
rm v1.3.zip

if [ -e scripts/SEACR-1.3 ]; then
	rm -r scripts/SEACR-1.3
fi

mv -f SEACR-1.3 src/
chmod +x src/SEACR-1.3/SEACR_1.3.sh
