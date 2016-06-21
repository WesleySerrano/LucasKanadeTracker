#!/bin/bash
rm -f *.o;

make;

if [ "$#" -lt 3 ]; then
	./main 20 jpg
elif [ "$#" -ge 4 ]; then
	./main "$1" "$2" "$3" "$4"
elif [ "$#" -ge 3 ]; then
	./main "$1" "$2" "$3"
else
	./main
fi