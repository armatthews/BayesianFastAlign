CC=g++
CPYP_DIR = /home/austinma/git/cpyp
INCS=-I$(CPYP_DIR)
LIBS=-L$(PREFIX)/lib
FINAL=-lboost_regex -lboost_program_options
#FINAL=-lcnn -lcnncuda -lboost_regex -lboost_serialization -lboost_program_options -lcuda -lcudart -lcublas
#CFLAGS=-std=c++11 -O3 -g -march=native -pipe # WE MUST NEVER USE -Ofast with CPYP!
CFLAGS=-std=c++11 -Wall -pedantic -O0 -g -pipe -DDEBUG
BINDIR=bin
OBJDIR=obj
SRCDIR=src

.PHONY: clean
all: make_dirs $(BINDIR)/mode $(BINDIR)/bayesianfastalign

make_dirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)

include $(wildcard $(OBJDIR)/*.d)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@
	$(CC) -MM -MP -MT "$@" $(CFLAGS) $(INCS) $< > $(OBJDIR)/$*.d

$(BINDIR)/bayesianfastalign: $(addprefix $(OBJDIR)/, bayesianfastalign.o)
	$(CC) $(CFLAGS) $(LIBS) $(INCS) $^ -o $@ $(FINAL)

$(BINDIR)/mode: $(addprefix $(OBJDIR)/, mode.o)
	$(CC) $(CFLAGS) $(LIBS) $(INCS) $^ -o $@ $(FINAL)

clean:
	rm -f $(BINDIR)/*
	rm -f $(OBJDIR)/*
