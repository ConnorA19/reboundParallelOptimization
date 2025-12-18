export OPENMP=1
export OPENGL=0# Set this to 1 to enable OpenGL
export SERVER=1# Set this to 1 to enable the visualization web server
include ../../src/Makefile.defs

# CCPROBLEM is defined in Makefile.defs to allow for
# a compact cross platform Makefile
.PHONY: all librebound
all: problem.c mockCollision.c librebound
	@echo "Compiling $< ..."
	$(CCPROBLEM) mockCollision.c speedupCollisionAttempts.c
# 	cc -I../../src/ -Wl,-rpath,./ -std=c99 -Wpointer-arith -D_GNU_SOURCE -O3 -fPIC -Wall -g  -fopenmp -DSERVER -DOPENMP -DGITHASH=16cd699dc0bc3aa32c879d3ae155a250ea279085 problem.c  -L. -lrebound -lm -lrt -fopenmp -o rebound mockCollision.c speedupCollisionAttempts.c
	@echo ""
	@echo "Compilation successful. To run REBOUND, execute the file '$(EXEREBOUND)'."
	@echo ""

test4: problem.c mockCollision.c librebound
	@echo "Compiling $< ..."
	$(CCPROBLEM) mockCollision.c speedupCollisionAttempts.c
	echo " "
	echo "4 Threads"
	OMP_NUM_THREADS=4 ./rebound


test: problem.c mockCollision.c librebound
	@echo "Compiling $< ..."
	$(CCPROBLEM) mockCollision.c speedupCollisionAttempts.c
	echo " "
	echo "1 Thread"
	OMP_NUM_THREADS=1 ./rebound
	echo " "
	echo "2 Threads"
	OMP_NUM_THREADS=2 ./rebound
	echo " "
	echo "4 Threads"
	OMP_NUM_THREADS=4 ./rebound
	echo " "
	echo "8 Threads"
	OMP_NUM_THREADS=8 ./rebound

librebound:
	@echo "Compiling shared library $(LIBREBOUND) ..."
	$(MAKE) -C ../../src/
	@-$(RM) $(LIBREBOUND)
	@$(LINKORCOPYLIBREBOUND)
	@echo ""


clean:
	@echo "Cleaning up shared library $(LIBREBOUND) ..."
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-$(RM) $(LIBREBOUND)
	@-$(RM) $(EXEREBOUND)

rebound_webgl.html: problem.c
	@echo "Compiling problem.c with emscripten (WebGL enabled)..."
	emcc -O3 -I../../src/ ../../src/*.c problem.c -DSERVERHIDEWARNING -DOPENGL=1 -sSTACK_SIZE=655360 -s USE_GLFW=3 -s FULL_ES3=1 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../web_client/shell_rebound_webgl.html -o rebound_webgl.html

rebound_console.html: problem.c
	@echo "Compiling problem.c with emscripten (WebGL disabled)..."
	emcc -O3 -I../../src/ ../../src/*.c problem.c -DSERVERHIDEWARNING -sSTACK_SIZE=655360 -sASYNCIFY -sALLOW_MEMORY_GROWTH -sSINGLE_FILE -sEXPORTED_RUNTIME_METHODS="callMain" --shell-file ../../web_client/shell_rebound_console.html -o rebound_console.html
