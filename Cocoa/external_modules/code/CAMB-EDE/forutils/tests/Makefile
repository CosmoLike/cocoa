
all: default

default: Release Debug

include ../Makefile_compiler

export

Debug:
	@$(MAKE) -C '..' $(MAKEOPT) Debug
	@$(MAKE) tester_Debug OUTPUT=Debug F90FLAGS="$(F90DEBUGFLAGS) $(F90COMMONFLAGS)" $(MAKEOPT)

Release:
	@$(MAKE) -C '..' $(MAKEOPT) Release
	@$(MAKE) tester_Release OUTPUT=Release F90FLAGS="$(F90RELEASEFLAGS) $(F90COMMONFLAGS)" $(MAKEOPT)

OUTPUT_DIR=../$(OUTPUT)

FILES = file_tests array_tests ObjectLists_tests interpolation_tests inifile_tests \
           string_tests ranges_tests run_tests

OBJS = $(patsubst %,$(OUTPUT_DIR)/%.o,$(FILES))

tester_$(OUTPUT): $(OUTPUT_DIR)/libforutils.a $(OBJS)
	$(F90C) $(F90FLAGS) $(OBJS) -o $@ -L$(OUTPUT_DIR) -lforutils

$(OUTPUT_DIR)/%.o: %.f90 $(OUTPUT_DIR)/compiler.ver
	$(F90C) $(F90FLAGS) $(MODLINK)$(OUTPUT_DIR) -c $< -o $(OUTPUT_DIR)/$*.o

run: Release Debug
	./tester_Debug
	./tester_Release
	@echo 'All Ok'

clean:
	@$(MAKE) -C '..' $(MAKEOPT) clean
	rm -f *.tmp

delete: clean
	rm -f tester_*

.PHONY: Debug Release clean test clean
