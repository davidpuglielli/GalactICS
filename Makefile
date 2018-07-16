CC=g++
CCC=nvcc
FLAGS=-ffast-math -O4 -funroll-loops -Wall -fbounds-check
CUDAFLAGS=--compiler-bindir=/home/david/gcc-bin -arch=sm_21
LIBS=-lm -lgsl -lgslcblas -L/usr/local/cuda/lib64 -lcudart
BIN=../bin-G
ALLFILES=GenerateGalaxy MCMC GetOneModel getPDFs
CONFIG=-DDEBUG

#.cu.o:
#	$(CCC) $(CUDAFLAGS) -c $*.cu

.cpp.o:
	$(CC) $(FLAGS) $(LIBS) $(CONFIG) -c $*.cpp -I/usr/local/cuda/include/

all: $(ALLFILES)
	
GenerateGalaxy: GenerateGalaxy.o cuda_core.o adiabaticcontraction.o allocatevectors.o appdiskforce.o bulgedenspsi.o bulgedf.o bulgepotential.o chisq.o dbh.o densrpsi.o diskdens.o diskdf.o diskdfez.o diskpotentialestimate.o diskprofile.o dofileio.o dpolardens.o force.o galaxy.o gasdiskdens.o gasdiskprofile.o genblackhole.o gendenspsibulge.o gendenspsihalo.o gendf.o getchisquare.o getfreqs.o getnbody.o getomegakappa.o getp.o getpsi.o gettrunc.o halodenspsi.o halodf.o haloforce.o halopotential.o halopotentialestimate.o haloprofiles.o mcmcfunctions.o nbody.o observations.o parameters.o pot.o rotation.o sersicprofiles.o sigrz.o splinetd.o totaldens.o
	$(CC) GenerateGalaxy.o cuda_core.o adiabaticcontraction.o allocatevectors.o appdiskforce.o bulgedenspsi.o bulgedf.o bulgepotential.o chisq.o dbh.o densrpsi.o diskdens.o diskdf.o diskdfez.o diskpotentialestimate.o diskprofile.o dofileio.o dpolardens.o force.o galaxy.o gasdiskdens.o gasdiskprofile.o genblackhole.o gendenspsibulge.o gendenspsihalo.o gendf.o getchisquare.o getfreqs.o getnbody.o getomegakappa.o getp.o getpsi.o gettrunc.o halodenspsi.o halodf.o haloforce.o halopotential.o halopotentialestimate.o haloprofiles.o mcmcfunctions.o nbody.o observations.o parameters.o pot.o rotation.o sersicprofiles.o sigrz.o splinetd.o totaldens.o -o GenerateGalaxy $(LIBS) $(FLAGS)

MCMC: mcmc.o cuda_core.o adiabaticcontraction.o allocatevectors.o appdiskforce.o bulgedenspsi.o bulgedf.o bulgepotential.o chisq.o dbh.o densrpsi.o diskdens.o diskdf.o diskdfez.o diskpotentialestimate.o diskprofile.o dofileio.o dpolardens.o force.o galaxy.o gasdiskdens.o gasdiskprofile.o genblackhole.o gendenspsibulge.o gendenspsihalo.o gendf.o getchisquare.o getfreqs.o getnbody.o getomegakappa.o getp.o getpsi.o gettrunc.o halodenspsi.o halodf.o haloforce.o halopotential.o halopotentialestimate.o haloprofiles.o mcmcfunctions.o nbody.o observations.o parameters.o pot.o rotation.o sersicprofiles.o sigrz.o splinetd.o totaldens.o
	$(CC) mcmc.o cuda_core.o adiabaticcontraction.o allocatevectors.o appdiskforce.o bulgedenspsi.o bulgedf.o bulgepotential.o chisq.o dbh.o densrpsi.o diskdens.o diskdf.o diskdfez.o diskpotentialestimate.o diskprofile.o dofileio.o dpolardens.o force.o galaxy.o gasdiskdens.o gasdiskprofile.o genblackhole.o gendenspsibulge.o gendenspsihalo.o gendf.o getchisquare.o getfreqs.o getnbody.o getomegakappa.o getp.o getpsi.o gettrunc.o halodenspsi.o halodf.o haloforce.o halopotential.o halopotentialestimate.o haloprofiles.o mcmcfunctions.o nbody.o observations.o parameters.o pot.o rotation.o sersicprofiles.o sigrz.o splinetd.o totaldens.o -o MCMC $(LIBS) $(FLAGS)

GetOneModel: GetOneModel.o cuda_core.o adiabaticcontraction.o allocatevectors.o appdiskforce.o bulgedenspsi.o bulgedf.o bulgepotential.o chisq.o dbh.o densrpsi.o diskdens.o diskdf.o diskdfez.o diskpotentialestimate.o diskprofile.o dofileio.o dpolardens.o force.o galaxy.o gasdiskdens.o gasdiskprofile.o genblackhole.o gendenspsibulge.o gendenspsihalo.o gendf.o getchisquare.o getfreqs.o getnbody.o getomegakappa.o getp.o getpsi.o gettrunc.o halodenspsi.o halodf.o haloforce.o halopotential.o halopotentialestimate.o haloprofiles.o mcmcfunctions.o nbody.o observations.o parameters.o pot.o rotation.o sersicprofiles.o sigrz.o splinetd.o totaldens.o
	$(CC) GetOneModel.o cuda_core.o adiabaticcontraction.o allocatevectors.o appdiskforce.o bulgedenspsi.o bulgedf.o bulgepotential.o chisq.o dbh.o densrpsi.o diskdens.o diskdf.o diskdfez.o diskpotentialestimate.o diskprofile.o dofileio.o dpolardens.o force.o galaxy.o gasdiskdens.o gasdiskprofile.o genblackhole.o gendenspsibulge.o gendenspsihalo.o gendf.o getchisquare.o getfreqs.o getnbody.o getomegakappa.o getp.o getpsi.o gettrunc.o halodenspsi.o halodf.o haloforce.o halopotential.o halopotentialestimate.o haloprofiles.o mcmcfunctions.o nbody.o observations.o parameters.o pot.o rotation.o sersicprofiles.o sigrz.o splinetd.o totaldens.o -o GetOneModel $(LIBS) $(FLAGS)

getPDFs: getPDFs.o cuda_core.o adiabaticcontraction.o allocatevectors.o appdiskforce.o bulgedenspsi.o bulgedf.o bulgepotential.o chisq.o dbh.o densrpsi.o diskdens.o diskdf.o diskdfez.o diskpotentialestimate.o diskprofile.o dofileio.o dpolardens.o force.o galaxy.o gasdiskdens.o gasdiskprofile.o genblackhole.o gendenspsibulge.o gendenspsihalo.o gendf.o getchisquare.o getfreqs.o getnbody.o getomegakappa.o getp.o getpsi.o gettrunc.o halodenspsi.o halodf.o haloforce.o halopotential.o halopotentialestimate.o haloprofiles.o mcmcfunctions.o nbody.o observations.o parameters.o pot.o rotation.o sersicprofiles.o sigrz.o splinetd.o totaldens.o
	$(CC) getPDFs.o cuda_core.o adiabaticcontraction.o allocatevectors.o appdiskforce.o bulgedenspsi.o bulgedf.o bulgepotential.o chisq.o dbh.o densrpsi.o diskdens.o diskdf.o diskdfez.o diskpotentialestimate.o diskprofile.o dofileio.o dpolardens.o force.o galaxy.o gasdiskdens.o gasdiskprofile.o genblackhole.o gendenspsibulge.o gendenspsihalo.o gendf.o getchisquare.o getfreqs.o getnbody.o getomegakappa.o getp.o getpsi.o gettrunc.o halodenspsi.o halodf.o haloforce.o halopotential.o halopotentialestimate.o haloprofiles.o mcmcfunctions.o nbody.o observations.o parameters.o pot.o rotation.o sersicprofiles.o sigrz.o splinetd.o totaldens.o -o getPDFs $(LIBS) $(FLAGS)

#allocatevectors.o: allocatevectors.cpp galaxy.h
#	$(CC) $(LIBS) -c allocatevectors.cpp $(FLAGS)
#
#appdiskforce.o: appdiskforce.cpp galaxy.h
#	$(CC) $(LIBS) -c appdiskforce.cpp $(FLAGS)
#        
#bulgedenspsi.o: bulgedenspsi.cpp galaxy.h
#	$(CC) $(LIBS) -c bulgedenspsi.cpp $(FLAGS)
#
#bulgedf.o: bulgedf.cpp galaxy.h
#	$(CC) $(LIBS) -c bulgedf.cpp $(FLAGS)
#
#bulgepotential.o: bulgepotential.cpp galaxy.h
#	$(CC) $(LIBS) -c bulgepotential.cpp $(FLAGS)
#
#chisq.o: chisq.cpp galaxy.h chisq.h
#	$(CC) $(LIBS) -c chisq.cpp $(FLAGS)
#        
#coef.o: coef.cpp galaxy.h
#	$(CC) $(LIBS) -c coef.cpp $(FLAGS)
#        
#dbh.o: dbh.cpp galaxy.h 
#	$(CC) $(LIBS) -c dbh.cpp $(FLAGS)
#        
#densrpsi.o: densrpsi.cpp galaxy.h
#	$(CC) $(LIBS) -c densrpsi.cpp $(FLAGS)
#
#diskdens.o: diskdens.cpp galaxy.h
#	$(CC) $(LIBS) -c diskdens.cpp $(FLAGS)
#
#diskdf.o: diskdf.cpp galaxy.h
#	$(CC) $(LIBS) -c diskdf.cpp $(FLAGS)
#
#diskdfez.o: diskdfez.cpp galaxy.h
#	$(CC) $(LIBS) -c diskdfez.cpp $(FLAGS)
#
#diskpotentialestimate.o: diskpotentialestimate.cpp galaxy.h
#	$(CC) $(LIBS) -c diskpotentialestimate.cpp $(FLAGS)
#
#diskprofile.o: diskprofile.cpp galaxy.h
#	$(CC) $(LIBS) -c diskprofile.cpp $(FLAGS)
#
#dofileio.o: dofileio.cpp galaxy.h
#	$(CC) $(LIBS) -c dofileio.cpp $(FLAGS)
#
#dpolardens.o: dpolardens.cpp galaxy.h
#	$(CC) $(LIBS) -c dpolardens.cpp $(FLAGS)
#
#galaxy.o: galaxy.cpp galaxy.h
#	$(CC) $(LIBS) -c galaxy.cpp $(FLAGS)
#
#genblackhole.o: genblackhole.cpp galaxy.h
#	$(CC) $(LIBS) -c genblackhole.cpp $(FLAGS)
#
#gendenspsibulge.o: gendenspsibulge.cpp galaxy.h
#	$(CC) $(LIBS) -c gendenspsibulge.cpp $(FLAGS)
#
#gendenspsihalo.o: gendenspsihalo.cpp galaxy.h
#	$(CC) $(LIBS) -c gendenspsihalo.cpp $(FLAGS)
#
#gendf.o: gendf.cpp galaxy.h
#	$(CC) $(LIBS) -c gendf.cpp $(FLAGS)
#
#GenerateGalaxy.o: GenerateGalaxy.cpp galaxy.h chisq.h
#	$(CC) $(LIBS) -c GenerateGalaxy.cpp $(FLAGS)
#
#getchisquare.o: getchisquare.cpp galaxy.h chisq.h
#	$(CC) $(LIBS) -c getchisquare.cpp $(FLAGS)
#
#getfreqs.o: getfreqs.cpp galaxy.h
#	$(CC) $(LIBS) -c getfreqs.cpp $(FLAGS)
#
#getnbody.o: getnbody.cpp galaxy.h
#	$(CC) $(LIBS) -c getnbody.cpp $(FLAGS)
#
#getomegakappa.o: getomegakappa.cpp galaxy.h
#	$(CC) $(LIBS) -c getomegakappa.cpp $(FLAGS)
#
#getpsi.o: getpsi.cpp galaxy.h
#	$(CC) $(LIBS) -c getpsi.cpp $(FLAGS)
#
#gettrunc.o: gettrunc.cpp galaxy.h
#	$(CC) $(LIBS) -c gettrunc.cpp $(FLAGS)
#
#halodenspsi.o: halodenspsi.cpp galaxy.h
#	$(CC) $(LIBS) -c halodenspsi.cpp $(FLAGS)
#
#halodf.o: halodf.cpp galaxy.h
#	$(CC) $(LIBS) -c halodf.cpp $(FLAGS)
#
#haloforce.o: haloforce.cpp galaxy.h
#	$(CC) $(LIBS) -c haloforce.cpp $(FLAGS)
#
#halopotential.o: halopotential.cpp galaxy.h
#	$(CC) $(LIBS) -c halopotential.cpp $(FLAGS)
#
#halopotentialestimate.o: halopotentialestimate.cpp galaxy.h
#	$(CC) $(LIBS) -c halopotentialestimate.cpp $(FLAGS)
#
#haloprofiles.o: haloprofiles.cpp galaxy.h
#	$(CC) $(LIBS) -c haloprofiles.cpp $(FLAGS)
#        
#mcmc.o: mcmc.cpp galaxy.h chisq.h
#	$(CC) $(LIBS) -c mcmc.cpp $(FLAGS)
#
#pot.o: pot.cpp galaxy.h
#	$(CC) $(LIBS) -c pot.cpp $(FLAGS)
#
#sersicprofiles.o: sersicprofiles.cpp galaxy.h
#	$(CC) $(LIBS) -c sersicprofiles.cpp $(FLAGS)
#
#sigrz.o: sigrz.cpp galaxy.h
#	$(CC) $(LIBS) -c sigrz.cpp $(FLAGS)
#
#splinetd.o: splinetd.cpp galaxy.h
#	$(CC) $(LIBS) -c splinetd.cpp $(FLAGS)
#
#totaldens.o: totaldens.cpp galaxy.h
#	$(CC) $(LIBS) -c totaldens.cpp $(FLAGS)

clean:
	rm *.o $(ALLFILES)

install:
	cp $(ALLFILES) ../bin-G
