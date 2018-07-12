all:
	g++ -o main main.c `root-config --cflags --glibs`
	g++ -o convert create_tree.cpp `root-config --cflags --glibs`
