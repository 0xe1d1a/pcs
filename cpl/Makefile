
CHPL = chpl
FLAGS = --fast

MUL_FLAGS = CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv

HEAT_PAR = heat
HEAT_MUL = heat_multi
heat_SOURCES = main.chpl common.chpl compute.chpl

EXAMPLE_PAR = example
EXAMPLE_MUL = example_multi
example_SOURCES = example.chpl

HELLO_LOC = hello
HELLO_DIST = hello_multi
hello_SOURCES = hello.chpl

all: $(CONF) $(heat_SOURCES) $(example_sources) $(hello_SOURCES) $(HEAT_SEQ) $(HEAT_PAR) $(HEAT_MUL) $(EXAMPLE_SEQ) $(EXAMPLE_PAR) $(EXAMPLE_MUL) $(HELLO_LOC) $(HELLO_DIST)

$(HEAT_PAR): $(heat_SOURCES)
	$(CHPL) -o $(HEAT_PAR) $(FLAGS) $(PAR_FLAGS) $(heat_SOURCES)

$(HEAT_MUL): $(heat_SOURCES)
	$(MUL_FLAGS) $(CHPL) -o $(HEAT_MUL) $(FLAGS) $(heat_SOURCES)

$(EXAMPLE_PAR): $(example_SOURCES)
	$(CHPL) -o $(EXAMPLE_PAR) $(FLAGS) $(PAR_FLAGS) $(example_SOURCES)

$(EXAMPLE_MUL): $(example_SOURCES)
	$(MUL_FLAGS) $(CHPL) -o $(EXAMPLE_MUL) $(FLAGS) $(example_SOURCES)

$(HELLO_LOC): $(hello_SOURCES)
	$(CHPL) -o $(HELLO_LOC) $(FLAGS) $(hello_SOURCES)

$(HELLO_DIST): $(hello_SOURCES)
	$(MUL_FLAGS) $(CHPL) -o $(HELLO_DIST) $(FLAGS) $(hello_SOURCES)

clean: 
	rm -f $(EXAMPLE_PAR) $(EXAMPLE_MUL) $(HEAT_PAR) $(HEAT_MUL) $(HEAT_MUL)_real $(EXAMPLE_MUL)_real $(HELLO_LOC){,_real} $(HELLO_DIST){,_real}

distcheck: guard-distdir
	tar -c --transform 's,^,$(distdir)/,' -zvf $(distdir).tar.gz common.chpl  compute.chpl  dist.job  example.chpl  hello.chpl  main.chpl  Makefile  set_env

guard-distdir:
	@ if [ "${distdir}" == "" ]; then \
        echo "Environment variable $* not set"; \
        echo "Run with: make distcheck distdir=STUDID"; \
        exit 1; \
    fi

