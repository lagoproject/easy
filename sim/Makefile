all:
	@echo "../bin/EasySim"
	@../bin/EasySim || true
	@for i in *.root; do nm=`basename $$i .root`; [ -f $$nm.dat ] || (echo "../bin/Ascii $$i > $$nm.dat" && ../bin/Ascii $$i > $$nm.dat && echo "Done"); done
	@for i in *.root ; do nm=`basename $$i .root`; [ $$nm.dat -nt $$i ] || (echo "../bin/Ascii $$i > $$nm.dat" && ../bin/Ascii $$i > $$nm.dat && echo "Done"); done

clean:
	rm -f calib_*.dat calib_*.root
