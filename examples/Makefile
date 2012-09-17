BASE = ..
include $(BASE)/Makefile.inc

all: hello fanin mergesort syseval

clean:
	rm -f *.o */*.o
	cd mpich; make clean
distclean: clean
	rm -f *~ hello hello_2p fanin fanin_2p fanout fanin_info fanout_info mergesort syseval syseval_1p syseval_2p testqueue testglobal test_tids snooze snoozepp
	cd mpich; make distclean

hello: hello.o $(LIBDEP)
	$(LINK) -o hello hello.o $(LIBMPI)
hello_2p: hello_2p.o $(LIBDEP)
	$(LINK) -o hello_2p hello_2p.o $(LIBMPI)

snooze: snooze.o $(LIBDEP)
	$(LINK) -o snooze snooze.o $(LIBMPI)
snoozepp: snooze.c $(PLIBDEP)
	../mpicc/mpicc -profile -echo -o snoozepp snooze.c ../../PUPPET/src/logger.o

fanin.o: fan.c
	$(CC) $(CFLAGS) -DISOUT=0 -c fan.c -o fanin.o
fanin: fanin.o $(LIBDEP)
	$(LINK) -o fanin fanin.o $(LIBMPI) -lm
fanpp: fan.c $(PLIBDEP)
	../mpicc/mpicc -profile -echo -o fanpp -DISOUT=0 fan.c ../../PUPPET/src/logger.o -lm
fanin_2p.o: fan_2p.c fan.c
	$(CC) $(CFLAGS) -DISOUT=0 -c fan_2p.c -o fanin_2p.o
fanin_2p: fanin_2p.o $(LIBDEP)
	$(LINK) -o fanin_2p fanin_2p.o $(LIBMPI) -lm
fanin_info.o: faninfo.c
	$(CC) $(CFLAGS) -DISOUT=0 -c faninfo.c -o fanin_info.o
fanin_info: fanin_info.o $(LIBDEP)
	$(LINK) -o fanin_info fanin_info.o $(LIBMPI) -lm
fanout.o: fan.c
	$(CC) $(CFLAGS) -DISOUT=1 -c fan.c -o fanout.o
fanout: fanout.o $(LIBDEP)
	$(LINK) -o fanout fanout.o $(LIBMPI) -lm
fanout_2p.o: fan_2p.c fan.c
	$(CC) $(CFLAGS) -DISOUT=1 -c fan_2p.c -o fanout_2p.o
fanout_2p: fanout_2p.o $(LIBDEP)
	$(LINK) -o fanout_2p fanout_2p.o $(LIBMPI) -lm
fanout_info.o: faninfo.c
	$(CC) $(CFLAGS) -DISOUT=1 -c faninfo.c -o fanout_info.o
fanout_info: fanout_info.o $(LIBDEP)
	$(LINK) -o fanout_info fanout_info.o $(LIBMPI) -lm

mergesort: mergesort.o $(LIBDEP)
	$(LINK) -o mergesort mergesort.o $(LIBMPI)
mergesortpp: mergesort.c $(PLIBDEP)
	../mpicc/mpicc -profile -echo -o mergesortpp -DISOUT=0 mergesort.c ../../PUPPET/src/logger.o -lm

syseval: syseval.o $(LIBDEP)
	$(LINK) -o syseval syseval.o $(LIBMPI)

syseval_1p: syseval_1p.o $(LIBDEP)
	$(LINK) -o syseval_1p syseval_1p.o $(LIBMPI)

syseval_2p: syseval_2p.o $(LIBDEP)
	$(LINK) -o syseval_2p syseval_2p.o $(LIBMPI)

testqueue: testqueue.o $(LIBDEP)
	$(LINK) -o testqueue testqueue.o $(LIBMPI)

testglobal: testglobal.c $(LIBDEP)
	../mpicc/mpicc $(CFLAGS) -o testglobal testglobal.c
test_tids: test_tids.c $(LIBDEP)
	$(CC) $(CFLAGS) -o test_tids test_tids.c $(LIBMPI)
