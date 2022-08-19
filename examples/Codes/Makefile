ARMA_INCLUDE = -I ~/libraries/armadillo/include
ARMA_LIB = -L ~/libraries/armadillo/lib64
QIC_INCLUDE= -I ~libraries/QIClib/include
BEC_INCLUDE= -I ~/libraries/BEClib/include

APP = test
name = $(APP).cpp
obj = $(APP).o
bin = $(APP)
CXX = mpicxx
CXXFLAGS = -Wall -O3 -std=c++17 -Wfatal-errors ${ARMA_INCLUDE} ${QIC_INCLUDE} ${BEC_INCLUDE}
LIBFLAGS = -larmadillo ${ARMA_LIB}

${bin}: ${obj}
	${CXX} ${obj} -o ${bin} ${CXXFLAGS} ${LIBFLAGS} 

${obj}: ${name}
	${CXX} -c ${name} ${CXXFLAGS}

.PHONY: clean
clean:
	-rm -f ${bin} ${obj}

#.PHONY: run
#run:
#	./${bin}
