SEL_DIR		:= src/jec4prompt/selections
SEL_FILES	:= $(wildcard $(SEL_DIR)/*/*.cpp)

$(info SEL_FILES: $(SEL_FILES))

.PHONY: all
.PHONY: $(SEL_FILES)

all: $(SEL_FILES)

install: $(SEL_FILES)
	python3 -m pip install -e .

$(SEL_FILES):
	root -e ".L $@++" -q


uninstall: clean
	python3 -m pip uninstall jec4prompt

clean:
	rm -fv $(SEL_DIR)/*/*_cpp*
