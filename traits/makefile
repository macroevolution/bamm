OBJS = main.o branchHistory.o mcmc.o event.o node.o Settings.o TraitModel.o TraitMCMC.o MbRandom.o Utilities.o TraitBranchHistory.o TraitEvent.o model.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c -m64 -O3 -funroll-loops $(DEBUG)
LFLAGS = -Wall -m64 $(DEBUG)

bammtrait: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o bammtrait

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

branchHistory.o: branchHistory.cpp branchHistory.h event.h Utilities.h
	$(CC) $(CFLAGS) branchHistory.cpp

mcmc.o: mcmc.cpp mcmc.h
	$(CC) $(CFLAGS) mcmc.cpp

event.o: event.cpp event.h
	$(CC) $(CFLAGS) event.cpp

node.o: node.cpp node.h
	$(CC) $(CFLAGS) node.cpp

Settings.o: Settings.cpp Settings.h
	$(CC) $(CFLAGS) Settings.cpp

TraitModel.o: TraitModel.cpp TraitModel.h TraitEvent.h
	$(CC) $(CFLAGS) TraitModel.cpp

TraitMCMC.o: TraitMCMC.cpp TraitMCMC.h
	$(CC) $(CFLAGS) TraitMCMC.cpp

MbRandom.o: MbRandom.cpp MbRandom.h
	$(CC) $(CFLAGS) MbRandom.cpp

Utilities.o: Utilities.cpp Utilities.h event.h TraitEvent.h
	$(CC) $(CFLAGS) Utilities.cpp

TraitBranchHistory.o: TraitBranchHistory.cpp TraitBranchHistory.h TraitEvent.h Utilities.h
	$(CC) $(CFLAGS) TraitBranchHistory.cpp

TraitEvent.o: TraitEvent.cpp TraitEvent.h
	$(CC) $(CFLAGS) TraitEvent.cpp

model.o: model.cpp model.h event.h
	$(CC) $(CFLAGS) model.cpp

clean:
	rm -rf *.o bammtrait