OBJS = main.o BranchHistory.o MCMC.o BranchEvent.o node.o Settings.o TraitModel.o TraitMCMC.o MbRandom.o Utilities.o TraitBranchHistory.o TraitBranchEvent.o Model.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c -m64 -O3 -funroll-loops $(DEBUG)
LFLAGS = -Wall -m64 $(DEBUG)

bamm: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o bamm

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

BranchHistory.o: BranchHistory.cpp BranchHistory.h BranchEvent.h Utilities.h
	$(CC) $(CFLAGS) BranchHistory.cpp

MCMC.o: MCMC.cpp MCMC.h
	$(CC) $(CFLAGS) MCMC.cpp

BranchEvent.o: BranchEvent.cpp BranchEvent.h
	$(CC) $(CFLAGS) BranchEvent.cpp

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

Utilities.o: Utilities.cpp Utilities.h BranchEvent.h TraitBranchEvent.h
	$(CC) $(CFLAGS) Utilities.cpp

TraitBranchHistory.o: TraitBranchHistory.cpp TraitBranchHistory.h TraitBranchEvent.h Utilities.h
	$(CC) $(CFLAGS) TraitBranchHistory.cpp

TraitBranchEvent.o: TraitBranchEvent.cpp TraitBranchEvent.h
	$(CC) $(CFLAGS) TraitBranchEvent.cpp

Model.o: Model.cpp Model.h BranchEvent.h
	$(CC) $(CFLAGS) Model.cpp

clean:
	rm -rf *.o bamm
