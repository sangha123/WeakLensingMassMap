.SUFFIXES: .C .o
CPP = g++
flags = -O3


lflags = -lm
path = /users/astro/sdeb
TNT = -I $(path)/TNT -I $(path)/JAMA
BIN_DIR = ../../bin

all: recon_pbl

clean:
	-rm *.o
	-rm *~

.C.o: $?
	$(CPP) $(TNT) $(flags) -c $*.C

recon_pbl: recon_pbl.o binnedgrid.o grid_field.o invert.o modelpbl.o datapbl.o
	$(CPP) recon_pbl.o binnedgrid.o grid_field.o invert.o modelpbl.o \
	datapbl.o  $(lflags) -o $(BIN_DIR)/recon_pbl

recon_pbl_bullet: recon_pbl_bullet.o binnedgrid.o grid_field.o invert.o \
	modelpbl.o datapbl.o
	$(CPP) recon_pbl_bullet.o binnedgrid.o grid_field.o invert.o \
	modelpbl.o datapbl.o  $(lflags) -o $(BIN_DIR)/recon_pbl_bullet