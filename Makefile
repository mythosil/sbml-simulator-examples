build:
	@g++ -std=c++11 -Wall simulator.cpp -o simulator -lsbml

sample1: build run/sample1 plot/sample1

mapk: build run/mapk plot/mapk

run/sample1:
	@./simulator models/sample1.xml 10.0 0.1 | tee result.sample1.dat

run/mapk:
	@./simulator models/MAPK.xml 1000.0 0.1 | tee result.mapk.dat

plot/sample1:
	@gnuplot -c gnuplot/sample1.plt

plot/mapk:
	@gnuplot -c gnuplot/mapk.plt

.PHONY: clean
clean:
	@rm -f simulator result.sample1.dat result.mapk.dat

