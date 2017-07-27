CC=g++
CPPFLAGS= -DCALIBONLY -Wno-deprecated -I`root-config --cflags` -Wall
LINKFLAGS= `root-config --libs`

SRC=Event.o Station.o Array.o Particle.o EventDict.o ShowerParam.o Utils.o BuildProcessesTables.o FastSimulation.o ManageEvent.o EasySimConfig.o Analyze.o AnaDict.o

TARGETS=EasySim Ascii Trace

all: $(TESTS) $(TARGETS)

EventDict.cc: Event.h Station.h Array.h Particle.h EasySimLinkDef.h
	rootcint -f EventDict.cc -c Event.h Station.h Array.h Particle.h EasySimLinkDef.h

AnaDict.cc: Analyze.h
	rootcint -f AnaDict.cc -c Analyze.h

.o:
	$(CC) -c $<

EasySim: ${SRC} EasySim.o
	$(CC) -o $@ ${SRC} EasySim.o ${LINKFLAGS}

Ascii: ${SRC} Ascii.o
	$(CC) -o $@ ${SRC} Ascii.o ${LINKFLAGS}

Trace: 
	$(CC) -o $@ trace.cc

clean:
	rm -f $(TARGETS) *.o *Dict*

install:
	mkdir -p ../bin/
	cp $(TARGETS) ../bin/

# dependencies added by makedep (and slightly alterated for ROOT dictionnaries)
# DO NOT DELETE

Analyze.o: Analyze.h Event.h Station.h Constants.h Particle.h Array.h Trigger.h
AnalyzeStar.o: Event.h Station.h Constants.h Particle.h Array.h Trigger.h
Array.o: Array.h Constants.h Station.h Particle.h EasySimConfig.h Utils.h
Ascii.o: Analyze.h Event.h Station.h Constants.h Particle.h Array.h
BuildProcessesTables.o: BuildProcessesTables.h Calib.h Constants.h EasySim.h EasySimConfig.h Utils.h
EasySim.o: EasySim.h EasySimConfig.h Utils.h Array.h Constants.h Station.h Particle.h Event.h Calib.h ManageEvent.h BuildProcessesTables.h ShowerParam.h
EasySimConfig.o: EasySimConfig.h Utils.h
Event.o: Event.h Station.h Constants.h Particle.h EasySim.h EasySimConfig.h Utils.h
FastSimulation.o: Constants.h
ManageEvent.o: ManageEvent.h Event.h Station.h Constants.h Particle.h ShowerParam.h Array.h EasySim.h EasySimConfig.h Utils.h Calib.h
Particle.o: Particle.h EasySimConfig.h Utils.h BuildProcessesTables.h Constants.h
ShowerParam.o: EasySim.h EasySimConfig.h Utils.h ShowerParam.h Constants.h
Station.o: Station.h Constants.h Particle.h Array.h Trigger.h ElecConstants.h ShowerParam.h EasySim.h EasySimConfig.h Utils.h
Utils.o: Utils.h
anamuons.o: Analyze.h Event.h Station.h Constants.h Particle.h Array.h
easysim2cdas.o: Event.h Station.h Constants.h Particle.h Trigger.h ElecConstants.h
