opt: main.c dijkstra.c fib.c
	gcc -g -Wall -std=c99 -lm -lnlopt -fopenmp -o bin/oopt main.c dijkstra.c fib.c

calc: calculate.c dijkstra.c fib.c
	gcc -g -Wall -std=c23 -lm -fopenmp -o bin/calculate calculate.c dijkstra.c fib.c

server: server.c
	gcc -g -Wall -std=c99 -lm -o bin/server server.c
