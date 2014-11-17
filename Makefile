DEBUG=-Wall -g -DDEBUG
LDFLAGS=-g -lpthread -O3
BINS=calc
CC=g++
NOGNUPLOT=-DNOGNUPLOT
FILE=main.cpp 
all: ${BINS}

calc: clean
		${CC} ${FILE} ${LDFLAGS} -o ${BINS}

debug: clean
		${CC} ${FILE} ${DEBUG} ${LDFLAGS} -o ${BINS}

nognuplot: clean
		${CC} ${FILE} ${LDFLAGS}  ${NOGNUPLOT} -o ${BINS}


clean:
			/bin/rm -rf ${BINS} *.o core *.core
