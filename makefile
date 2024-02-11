#####################################################################
# Makefile for the parallel NUFFT project
# E Boström 240209
#####################################################################

# Define shell
SHELL = /bin/sh

# Directorys
SRC_DIR := ./src
BUILD_DIR := ./build
TEST_DIR := ./test
INCL_DIR := ./include

# Compiler, flags and libraries
CC = g++
CFLAGS = -I$(INCL_DIR) -g -Wall
DEPFLAGS = -MMD -MP
LIBS = -fopenmp -lm -lfftw3

# Specify the maximum number of OpenMP threads you want to use
MAX_OMP_THREADS = 12

# Get time. Used to generate a time stamp
TIMENOW = $(shell date +%y%m%d%H%M%S)

# Targets
_TEST_TARGETS = test_naive_nufft3d test_gridding test_FGG test_parallel
TEST_TARGETS = $(addprefix $(TEST_DIR)/, $(_TEST_TARGETS))

# Files
SRC_FILES = gridding3d.cpp sorting.cpp dirsum3d.cpp deconv3d.cpp utils.cpp
OBJ_FILES := $(addprefix $(BUILD_DIR)/, $(SRC_FILES:.cpp=.o))
DEP_FILES := $(addprefix $(BUILD_DIR)/, $(SRC_FILES:.cpp=.d))

.PHONY: clean

#--------------------------------------------------------------------
# Build rules
#--------------------------------------------------------------------

all: build_object_files build_tests

build_object_files : $(OBJ_FILES)

# Build object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -c -o $@ $< $(DEPFLAGS) $(LIBS)

-include $(DEP_FILES)

# Test programs
build_tests : $(TEST_TARGETS)
$(TEST_TARGETS) : $(TEST_TARGETS:=.cpp) $(OBJ_FILES)
	$(CC) $(CFLAGS) $@.cpp -o $@ $(OBJ_FILES) $(LIBS)


#--------------------------------------------------------------------
# Run test scripts, data is stored in .dat files
#--------------------------------------------------------------------
run_serial_tests: $(TEST_DIR)/test_FGG
	echo "Runing all serial tests"
	(cd $(TEST_DIR) ; sh run_serial_tests.sh)

run_parallel_tests: $(TEST_DIR)/test_parallel
	echo "Runing parallel speedup tests"
	(cd $(TEST_DIR) ; sh run_parallel_tests.sh)


#--------------------------------------------------------------------
# Create plots from generated data
#--------------------------------------------------------------------
serial_plots:
	gnuplot $(TEST_DIR)/plot_fgg_time.gnu
	gnuplot $(TEST_DIR)/plot_fgg_speedup.gnu

parallel_plots:
	gnuplot $(TEST_DIR)/plot_parallel_time.gnu
	gnuplot $(TEST_DIR)/plot_parallel_speedup.gnu			

# Backup data into a separate data folder
backup_data:
	@mkdir -p $(TEST_DIR)/data_$(TIMENOW)
	@cp $(TEST_DIR)/*.dat $(TEST_DIR)/data_$(TIMENOW)/.


#--------------------------------------------------------------------
# Cleanup rules
#--------------------------------------------------------------------

# Clean data that is not backed-up
clean_data:
	$(RM) $(TEST_DIR)/*.dat


# Clean images
clean_imgs:
	$(RM) $(TEST_DIR)/*.eps

# Clean all build files
clean_build:
	$(RM) -r $(BUILD_DIR)

clean_tests:
	$(RM) $(TEST_TARGETS)

clean: clean_build clean_tests

# Clean all (except backed up data)
clean_all: clean_data clean_imgs clean_build clean_tests