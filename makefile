COMPILER = g++

main: main.cpp
	$(COMPILER) -o main main.cpp -lX11 -lpthread -lm