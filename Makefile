BIN=HybridAttack
SRCDIR=src
SRCS=$(wildcard $(SRCDIR)/*.cpp)
OBJS=$(SRCS:$(SRCDIR)/%.cpp=$(SRCDIR)/%.o)
DEPS=$(SRCS:$(SRCDIR)/%.cpp=$(SRCDIR)/%.d)

CC = mpic++
#CC = g++ 
#g++

CPPFLAGS += -DUSING_TBB -DUSING_MPI
#
CFLAGS=  -fopenmp -g -Wfatal-errors -march=native -std=c++11  -Ofast
#-Ofast

LFLAGS= -fopenmp -g -lntl -lgmp -ltbb -ltbbmalloc_proxy -ltbbmalloc

INCLUDE=-I./include 

all: $(BIN)

$(BIN): $(OBJS)

	$(CC) $(INCLUDE) $^ -o $@ $(LFLAGS)

%.o: %.cpp

	$(CC) $(INCLUDE) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@ 

clean:

	rm $(SRCDIR)/*.d $(SRCDIR)/*.o $(BIN)

echo:

	echo $(OBJS)

-include $(DEPS)

.PHONY: all clean

