CC=chpl
CC_FLAGS=

.PHONEY: all
all: pi

SRC=$(wildcard *.chpl)

pi: $(SRC)
	@$(CC) $(CC_FLAGS) $(SRC) -o pi

.PHONEY: clean
clean:
	@rm pi

.PHONEY: run
run: all
	@echo "Sequential:"
	@./pi
	@echo
