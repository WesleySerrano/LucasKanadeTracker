COMPILER = g++

main: main.cpp utils.h
	$(COMPILER) -g -o main main.cpp utils.h -lX11 -lpthread -lm