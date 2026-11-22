CC      = clang
CFLAGS  = -O3 -Wall -Wextra -std=c11
LDFLAGS = -lm
TARGET  = born_microscopic_sim
SRC     = born_microscopic_sim.c

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.csv
