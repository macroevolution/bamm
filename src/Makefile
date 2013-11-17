OBJS = main.o Autotune.o BranchHistory.o MCMC.o BranchEvent.o Node.o Tree.o Settings.o TraitModel.o TraitMCMC.o MbRandom.o TraitBranchHistory.o TraitBranchEvent.o Model.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c -m64 -O3 -funroll-loops $(DEBUG)
LFLAGS = -Wall -m64 $(DEBUG)

bamm: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o bamm

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

Autotune.o: Autotune.cpp MbRandom.h Model.h TraitModel.h Settings.h
	$(CC) $(CFLAGS) Autotune.cpp
	
BranchHistory.o: BranchHistory.cpp BranchHistory.h BranchEvent.h Utilities.h
	$(CC) $(CFLAGS) BranchHistory.cpp

MCMC.o: MCMC.cpp MCMC.h
	$(CC) $(CFLAGS) MCMC.cpp

BranchEvent.o: BranchEvent.cpp BranchEvent.h
	$(CC) $(CFLAGS) BranchEvent.cpp

Node.o: Node.cpp Node.h
	$(CC) $(CFLAGS) Node.cpp

Tree.o: Tree.cpp Tree.h
	$(CC) $(CFLAGS) Tree.cpp

Settings.o: Settings.cpp Settings.h
	$(CC) $(CFLAGS) Settings.cpp

TraitModel.o: TraitModel.cpp TraitModel.h TraitBranchEvent.h
	$(CC) $(CFLAGS) TraitModel.cpp

TraitMCMC.o: TraitMCMC.cpp TraitMCMC.h
	$(CC) $(CFLAGS) TraitMCMC.cpp

MbRandom.o: MbRandom.cpp MbRandom.h
	$(CC) $(CFLAGS) MbRandom.cpp

TraitBranchHistory.o: TraitBranchHistory.cpp TraitBranchHistory.h TraitBranchEvent.h Utilities.h
	$(CC) $(CFLAGS) TraitBranchHistory.cpp

TraitBranchEvent.o: TraitBranchEvent.cpp TraitBranchEvent.h
	$(CC) $(CFLAGS) TraitBranchEvent.cpp

Model.o: Model.cpp Model.h BranchEvent.h
	$(CC) $(CFLAGS) Model.cpp

clean:
	rm -rf *.o bamm
