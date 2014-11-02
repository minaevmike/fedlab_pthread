DEBUG=-Wall -g -DDEBUG
LDFLAGS=-g -lpthread
BINS=calc
CC=g++
FILE=main.cpp 
all: ${BINS}

calc: clean
		${CC} ${FILE} ${LDFLAGS} -o ${BINS}

debug: clean
		${CC} ${FILE} ${DEBUG} ${LDFLAGS} -o ${BINS}

clean:
			/bin/rm -rf ${BINS} *.o core *.core
