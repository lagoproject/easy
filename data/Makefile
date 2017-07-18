all:
	cd ../src && make
	@echo "../src/EasySim"
	@../src/EasySim || true
	@for i in *.root ; do nm=`basename $$i .root`; [ -f $$nm.dat ] || (echo "../src/Ascii $$i > $$nm.dat" && ../src/Ascii $$i > $$nm.dat && echo "Done"); done
	@for i in *.root ; do nm=`basename $$i .root`; [ $$nm.dat -nt $$i ] || (echo "../src/Ascii $$i > $$nm.dat" && ../src/Ascii $$i > $$nm.dat && echo "Done"); done
	@rm *.root

clean:
	rm -f calib_*.dat calib_*.root
