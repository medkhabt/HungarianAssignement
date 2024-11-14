SRC := main.cpp
BUILD_DIR := build
DEBUG_DIR := $(BUILD_DIR)/debug
RELEASE_DIR := $(BUILD_DIR)/release

all: release debug 

release: $(RELEASE_DIR)/main

$(RELEASE_DIR)/main: $(SRC)
	mkdir -p $(RELEASE_DIR)
	g++ -std=c++14 -O2 -o $@ $< 

debug: $(DEBUG_DIR)/main

$(DEBUG_DIR)/main: $(SRC)
	mkdir -p $(DEBUG_DIR)
	g++ -std=c++14 -g -o $@ $< 

clean:
	rm -rf $(BUILD_DIR)
