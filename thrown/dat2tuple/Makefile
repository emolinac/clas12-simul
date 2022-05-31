## GENERAL 
GXX := g++

SRC  := ./src
INC  := ./include
BIN  := ./bin

NAME := dat2tuple

## ROOT related stuff
ROOT_CONFIG  := root-config
# ROOT dirs
ROOT_INCDIR  := $(shell $(ROOT_CONFIG) --incdir)
# ROOT flags
ROOT_CFLAGS  := $(shell $(ROOT_CONFIG) --cflags)
# ROOT libs
ROOT_LIBS    := $(shell $(ROOT_CONFIG) --libs)

## SHOWTIME
${BIN}/${NAME}: ${SRC}/${NAME}.cpp
	${GXX} ${SRC}/${NAME}.cpp -o ${BIN}/${NAME} -I${INC} ${ROOT_CFLAGS} ${ROOT_LIBS}

clean:
	rm ${BIN}/${NAME}
