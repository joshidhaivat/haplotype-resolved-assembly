# Compiler settings
CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -fopenmp -pthread -Wall -Wextra -Iinclude
LDFLAGS = -fopenmp -pthread -lz

# Directories
SRC_DIR = src
INC_DIR = include
OBJ_DIR = obj
BIN_DIR = bin

# Target executables
TARGET = $(BIN_DIR)/hmm_decoder
JFTOOL = $(BIN_DIR)/JFToOl

# Source files for main program
MAIN_SOURCES = $(SRC_DIR)/main.cpp \
               $(INC_DIR)/ol_worker.cpp \
               $(INC_DIR)/hmm_functions.cpp \
               $(INC_DIR)/pipeline_workers.cpp \
               $(INC_DIR)/hmm_pipeline.cpp \
               $(INC_DIR)/output_writer.cpp

# Source file for JFToOl converter
JFTOOL_SOURCES = $(INC_DIR)/JFToOl.cpp

# Header files
HEADERS = $(INC_DIR)/thread_safe_queue.h \
          $(INC_DIR)/math_utils.h \
          $(INC_DIR)/matrix.h \
          $(INC_DIR)/data_structures.h \
          $(INC_DIR)/hmm_functions.h \
          $(INC_DIR)/fasta_reader.h \
          $(INC_DIR)/ol_worker.h \
          $(INC_DIR)/pipeline_workers.h \
          $(INC_DIR)/hmm_pipeline.h \
          $(INC_DIR)/output_writer.h

# Object files for main program
MAIN_OBJECTS = $(OBJ_DIR)/main.o \
               $(OBJ_DIR)/ol_worker.o \
               $(OBJ_DIR)/hmm_functions.o \
               $(OBJ_DIR)/pipeline_workers.o \
               $(OBJ_DIR)/hmm_pipeline.o \
               $(OBJ_DIR)/output_writer.o

# Object files for JFToOl converter
JFTOOL_OBJECTS = $(OBJ_DIR)/JFToOl.o

# Default target - build both executables
all: directories $(TARGET) $(JFTOOL)

# Create necessary directories
directories:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(BIN_DIR)

# Link object files to create main executable
$(TARGET): $(MAIN_OBJECTS)
	@echo "Linking $(TARGET)..."
	$(CXX) $(MAIN_OBJECTS) -o $(TARGET) $(LDFLAGS)
	@echo "Build complete: $(TARGET)"

# Link object files to create JFToOl converter
$(JFTOOL): $(JFTOOL_OBJECTS)
	@echo "Linking $(JFTOOL)..."
	$(CXX) $(JFTOOL_OBJECTS) -o $(JFTOOL) $(LDFLAGS)
	@echo "Build complete: $(JFTOOL)"

# Compile main.cpp from src directory
$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp $(HEADERS)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile ol_worker.cpp from src directory
$(OBJ_DIR)/ol_worker.o: $(INC_DIR)/ol_worker.cpp $(INC_DIR)/ol_worker.h $(HEADERS)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile cpp files from include directory (except JFToOl.cpp)
$(OBJ_DIR)/hmm_functions.o: $(INC_DIR)/hmm_functions.cpp $(HEADERS)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/pipeline_workers.o: $(INC_DIR)/pipeline_workers.cpp $(HEADERS)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/hmm_pipeline.o: $(INC_DIR)/hmm_pipeline.cpp $(HEADERS)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/output_writer.o: $(INC_DIR)/output_writer.cpp $(HEADERS)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile JFToOl.cpp separately (has its own main())
$(OBJ_DIR)/JFToOl.o: $(INC_DIR)/JFToOl.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -rf $(OBJ_DIR) $(BIN_DIR)
	@echo "Clean complete"

# Clean and rebuild
rebuild: clean all

# Install both executables (optional - copies to /usr/local/bin)
install: $(TARGET) $(JFTOOL)
	@echo "Installing $(TARGET) and $(JFTOOL) to /usr/local/bin..."
	@sudo cp $(TARGET) /usr/local/bin/
	@sudo cp $(JFTOOL) /usr/local/bin/
	@echo "Installation complete"

# Uninstall (optional)
uninstall:
	@echo "Uninstalling hmm_decoder and JFToOl from /usr/local/bin..."
	@sudo rm -f /usr/local/bin/hmm_decoder
	@sudo rm -f /usr/local/bin/JFToOl
	@echo "Uninstall complete"

# Run with default test parameters (modify as needed)
test: $(TARGET) $(JFTOOL)
	@echo "Running test..."
	./$(TARGET) test/params.txt test/sequences.fasta test/counts.jf

# Display help
help:
	@echo "HMM Decoder Makefile"
	@echo "===================="
	@echo "Available targets:"
	@echo "  all       - Build the project (default) - creates hmm_decoder and JFToOl"
	@echo "  clean     - Remove build artifacts"
	@echo "  rebuild   - Clean and rebuild"
	@echo "  install   - Install to /usr/local/bin (requires sudo)"
	@echo "  uninstall - Remove from /usr/local/bin (requires sudo)"
	@echo "  test      - Run with test parameters"
	@echo "  help      - Display this help message"
	@echo ""
	@echo "Usage:"
	@echo "  make [target]"
	@echo ""
	@echo "After building, run with:"
	@echo "  ./bin/hmm_decoder <params> <fasta> <jellyfish_db.jf> [output] [threads] [batch_size]"
	@echo ""
	@echo "The program will automatically:"
	@echo "  1. Check if OL database exists (jellyfish_db.jf.prefix, .suffix, .counts)"
	@echo "  2. Build it if missing using ./bin/JFToOl"
	@echo "  3. Run the HMM pipeline"
	@echo ""
	@echo "Examples:"
	@echo "  ./bin/hmm_decoder params.txt reads.fasta counts.jf"
	@echo "  ./bin/hmm_decoder params.txt reads.fasta counts.jf output.txt 32 5000"
	@echo "  ./bin/hmm_decoder params.txt reads.fasta counts.jf --rebuild"

.PHONY: all directories clean rebuild install uninstall test help