.SUFFIXES: .C .o
CPP = g++
flags = -O3


lflags = -lm
TNT_PATH = /Users/sanghamitradeb/TNT/
JAMA_PATH = /Users/sanghamitradeb/JAMA/
INC = -I $(TNT_PATH) -I $(JAMA_PATH)
BIN_DIR = ../../bin

all: recon_pbl

clean:
	-rm *.o
	-rm *~

.C.o: $?
	$(CPP) $(INC) $(flags) -c $*.C

recon_pbl: recon_pbl.o binnedgrid.o grid_field.o invert.o modelpbl.o datapbl.o
	$(CPP) recon_pbl.o binnedgrid.o grid_field.o invert.o modelpbl.o \
	datapbl.o  $(lflags) -o $(BIN_DIR)/recon_pbl

recon_pbl_bullet: recon_pbl_bullet.o binnedgrid.o grid_field.o invert.o \
	modelpbl.o datapbl.o
	$(CPP) recon_pbl_bullet.o binnedgrid.o grid_field.o invert.o \
	modelpbl.o datapbl.o  $(lflags) -o $(BIN_DIR)/recon_pbl_bullet