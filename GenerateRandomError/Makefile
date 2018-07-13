BIN=GenerateInputHybridAttackRandomError
SRCDIR=src
SRCS=$(wildcard $(SRCDIR)/*.cpp)
OBJS=$(SRCS:$(SRCDIR)/%.cpp=$(SRCDIR)/%.o)
DEPS=$(SRCS:$(SRCDIR)/%.cpp=$(SRCDIR)/%.d)

CC = g++

CFLAGS= -fopenmp -g -Wfatal-errors -march=native -std=c++11  -Ofast
#-Ofast

LFLAGS= -fopenmp -g -lntl -lgmp

INCLUDE=-I./include 

all: $(BIN)

$(BIN): $(OBJS)

	$(CC) $(INCLUDE) $^ -o $@ $(LFLAGS)

%.o: %.cpp

	$(CC) $(INCLUDE) $(CFLAGS) -MMD -MP -c $< -o $@ 

clean:

	rm $(SRCDIR)/*.d $(SRCDIR)/*.o $(BIN)

echo:

	echo $(OBJS)

-include $(DEPS)

.PHONY: all clean

