all: 
	nvcc --expt-relaxed-constexpr parallel.cu -o parallel

debug: 
	nvcc -g -G parallel.cu -o parallel

clean:
	rm -f parallel
