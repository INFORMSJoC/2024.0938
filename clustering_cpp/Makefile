gurobi_path = /home/antonio/gurobi1100/linux64

sources = main.cpp cg_reduction.cpp cg_stabilization.cpp cg_util.cpp cg_pricing.cpp CGStruct.cpp Kmeans2D.cpp cg_pricing_multi.cpp cg_branching.cpp branch_and_price.cpp JobQueue.cpp ThreadPool.cpp
objects = $(sources:.cpp=.o)

all: bb

bb: $(objects)
	g++ -L $(gurobi_path)/lib/ -o $@ $^ -pthread -lgurobi_c++ -lgurobi110 -Wl,-rpath,$(gurobi_path)/lib/

%.o: %.cpp
		g++ -g -Wall -O2 -o $@ -c $< -I $(gurobi_path)/include/


.PHONY: clean
		
clean:
		rm -f bb *.o
