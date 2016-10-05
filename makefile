all: raycast.c
	gcc -o raycast raycast.c -lm

clean:
	rm -rf raycast.c *~
