COMPILER = g++

main: main.cpp
	$(COMPILER) -g -o main main.cpp -lX11 -lpthread -lm