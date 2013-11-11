OBJS = main.o BranchHistory.o mcmc.o event.o node.o Settings.o TraitModel.o TraitMCMC.o MbRandom.o Utilities.o TraitBranchHistory.o TraitBranchEvent.o model.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c -m64 -O3 -funroll-loops $(DEBUG)
LFLAGS = -Wall -m64 $(DEBUG)

bamm: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o bamm

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

BranchHistory.o: BranchHistory.cpp BranchHistory.h event.h Utilities.h
	$(CC) $(CFLAGS) BranchHistory.cpp

mcmc.o: mcmc.cpp mcmc.h
	$(CC) $(CFLAGS) mcmc.cpp

event.o: event.cpp event.h
	$(CC) $(CFLAGS) event.cpp

node.o: node.cpp node.h
	$(CC) $(CFLAGS) node.cpp

Settings.o: Settings.cpp Settings.h
	$(CC) $(CFLAGS) Settings.cpp

TraitModel.o: TraitModel.cpp TraitModel.h TraitBranchEvent.h
	$(CC) $(CFLAGS) TraitModel.cpp

TraitMCMC.o: TraitMCMC.cpp TraitMCMC.h
	$(CC) $(CFLAGS) TraitMCMC.cpp

MbRandom.o: MbRandom.cpp MbRandom.h
	$(CC) $(CFLAGS) MbRandom.cpp

Utilities.o: Utilities.cpp Utilities.h event.h TraitBranchEvent.h
	$(CC) $(CFLAGS) Utilities.cpp

TraitBranchHistory.o: TraitBranchHistory.cpp TraitBranchHistory.h TraitBranchEvent.h Utilities.h
	$(CC) $(CFLAGS) TraitBranchHistory.cpp

TraitBranchEvent.o: TraitBranchEvent.cpp TraitBranchEvent.h
	$(CC) $(CFLAGS) TraitBranchEvent.cpp

model.o: model.cpp model.h event.h
	$(CC) $(CFLAGS) model.cpp

clean:
	rm -rf *.o bamm
