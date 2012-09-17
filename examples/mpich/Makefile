BASE = ../..
include $(BASE)/Makefile.inc

CC = mpicc_shmem
#CC = /home/eddemain/Packages/MPICH/bin/mpicc
CFLAGS = $(INCLUDES) $(OPTFLAGS)

all: fanin_smp fanin_p4 fanout_smp fanout_p4 sfanin_smp sfanin_p4 sfanout_smp sfanout_p4 syseval_smp syseval_p4 ssyseval_smp ssyseval_p4

clean:
	rm -f *.o
distclean: clean
	rm -f *~ fanin_smp fanin_p4 fanout_smp fanout_p4 sfanin_smp sfanin_p4 sfanout_smp sfanout_p4 syseval_smp syseval_p4 ssyseval_smp ssyseval_p4

fanin.o: fan.c
	$(CC) $(CFLAGS) -DISOUT=0 -c fan.c -o fanin.o
fanin_smp: fanin.o $(LIBDEP)
	$(CC) -o fanin_smp fanin.o -lm
fanin_p4: fanin.o $(LIBDEP)
	mpicc_p4 -o fanin_p4 fanin.o -lm
fanout.o: fan.c
	$(CC) $(CFLAGS) -DISOUT=1 -c fan.c -o fanout.o
fanout_smp: fanout.o $(LIBDEP)
	$(CC) -o fanout_smp fanout.o -lm
sfanin.o: sfan.c
	$(CC) $(CFLAGS) -DISOUT=0 -c sfan.c -o sfanin.o
sfanin_smp: sfanin.o $(LIBDEP)
	$(CC) -o sfanin_smp sfanin.o -lm
sfanout.o: sfan.c
	$(CC) $(CFLAGS) -DISOUT=1 -c sfan.c -o sfanout.o
sfanout_smp: sfanout.o $(LIBDEP)
	$(CC) -o sfanout_smp sfanout.o -lm

syseval_p4: syseval.o
	mpicc_p4 -o syseval_p4 syseval.o

syseval_smp: syseval.o
	mpicc_shmem -o syseval_smp syseval.o

ssyseval_p4: ssyseval.o
	mpicc_p4 -o ssyseval_p4 ssyseval.o

ssyseval_smp: ssyseval.o
	mpicc_shmem -o ssyseval_smp ssyseval.o

