#OPTIONS
FLAGS = #-DPRINTBELIEF
#COMPILER
CC = g++
CFLAGS = -O3 -std=c++11 #-fsanitize=address -Wall -O0 -g #DEBUG
CFLAGS += -I/opt/local/include  -L/opt/local/lib -lm

#FILES
INCLUDES =

default: search post

search: ../src/main.cc
	${CC} -o multiagent_search ../src/main.cc ${CFLAGS} ${INCLUDES} ${FLAGS}

clean:
	rm *.o 
	@echo "Cleaned."

post:
	@echo "compiled"
	@echo "output is: multiagent_search"

