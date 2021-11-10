# Object files to either reference or create                          
OBJ = $(wildcard ./C++/*.cpp)

# The executable file that will be created at the end                 
EXEC = ./run/UDASF

# The flags to use for compilation     
FLAGS = -fexceptions               \
        -O2                        \
        -fexpensive-optimizations  \
        -Os                        \
        -O3                        \
        -O2                        \
        -O1                        \
        -O                         \
        -Wformat                   \
        -Werror=format-security    \
        -D_FORTIFY_SOURCE=2        \
        -g                         \
        -fconcepts                 \
        -Wall

# The code compiler to use for compilation                            
CC = g++ -std=c++2a

# Libraries
LIBS = -Wl,--no-export-dynamic     \
       -fopenmp

# Includes
INC = -I"./C++/"

# make
all: info command                                        
	
info:
	@echo "Compilando codigo - makefile linux. Aguarde ..."

command:
	$(CC) $(INC) $(FLAGS) $(OBJ) $(LIBS) -o $(EXEC)
	@echo "Executavel construido com sucesso"

