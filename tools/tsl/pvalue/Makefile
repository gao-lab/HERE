# Two Sample Logo
#
# Vladimir Vacic, University of California, Riverside
# Lilia M. Iakoucheva, Rockefeller University, New York
# Predrag Radivojac, Indiana University, Bloomington

CC=gcc
LIBS=-lm
FLAGS=
OBJS=bdtr.o const.o exp.o gamma.o incbet.o incbi.o isnan.o log.o mtherr.o ndtri.o polevl.o pow.o powi.o fabs.o

all: pvalue


bdtr.o: bdtr.c
	$(CC) $(FLAGS) bdtr.c -c

const.o: const.c
	$(CC) $(FLAGS) const.c -c

exp.o: exp.c
	$(CC) $(FLAGS) exp.c -c

gamma.o: gamma.c
	$(CC) $(FLAGS) gamma.c -c

incbet.o: incbet.c
	$(CC) $(FLAGS) incbet.c -c

incbi.o: incbi.c
	$(CC) $(FLAGS) incbi.c -c

isnan.o: isnan.c
	$(CC) $(FLAGS) isnan.c -c

log.o: log.c
	$(CC) $(FLAGS) log.c -c

ndtri.o: ndtri.c
	$(CC) $(FLAGS) ndtri.c -c

mtherr.o: mtherr.c
	$(CC) $(FLAGS) mtherr.c -c

polevl.o: polevl.c
	$(CC) $(FLAGS) polevl.c -c

pow.o: pow.c
	$(CC) $(FLAGS) pow.c -c

powi.o: powi.c
	$(CC) $(FLAGS) powi.c -c

fabs.o: fabs.c
	$(CC) $(FLAGS) fabs.c -c

pvalue: ${OBJS} pvalue.c
	$(CC) $(FLAGS) -o pvalue pvalue.c $(LIBS) ${OBJS}

clean:
	rm -f pvalue ${OBJS}

