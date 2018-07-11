all:
	g++ -o main main.c `root-config --cflags --glibs`
