#set default architecture, can be overridden from the compile line
#ARCH = $(VLASIATOR_ARCH)
ARCH = meteo
include MAKE/Makefile.${ARCH}

#set FP precision to SP (single) or DP (double)
FP_PRECISION = DP

#set a default archive utility, can also be set in Makefile.arch
AR ?= ar

#leveque (no other options)
TRANSSOLVER ?= leveque

#leveque or semilag
ACCSOLVER ?= leveque

#londrillo_delzanna (no other options)
FIELDSOLVER ?= londrillo_delzanna

CXXFLAGS += -DPROFILE

#Add -DNDEBUG to turn debugging off. If debugging is enabled performance will degrade significantly
CXXFLAGS += -DNDEBUG

#define USE_AGNER_VECTORCLASS to use an external vector class that is used in some of the solvers
#If not defined a slower but portable implementation is used, as the external one only supports 
#Linux & x86 processors  
CXXFLAGS += -DUSE_AGNER_VECTORCLASS
#Add -DCATCH_FPE to catch floating point exceptions and stop execution
#May cause problems
# CXXFLAGS += -DCATCH_FPE


#Add -DFS_1ST_ORDER_SPACE or -DFS_1ST_ORDER_TIME to make the field solver first-order in space or time
# CXXFLAGS += -DFS_1ST_ORDER_SPACE
# CXXFLAGS += -DFS_1ST_ORDER_TIME


#//////////////////////////////////////////////////////
# The rest of this file users shouldn't need to change
#//////////////////////////////////////////////////////


ifeq ($(strip $(ACCSOLVER)),semilag)
# If defined accelerates velocity space using cpu_acc_semilag.hpp.
# Otherwise the same solver is used in real and velocity spaces.
CXXFLAGS += -DSEMILAG
endif

ifeq ($(strip $(TRANSSOLVER)),leveque)
CXXFLAGS += -DSOLVER_LEVEQUE
endif

#will need profiler in most places..
CXXFLAGS += ${INC_PROFILE} 

#define precision
CXXFLAGS += -D${FP_PRECISION} 

CXXEXTRAFLAGS = ${CXXFLAGS} -DTOOL_NOT_PARALLEL


default: vlasiator

tools: parallel_tools not_parallel_tools

parallel_tools: vlsv2vtk vlsv2silo vlsv2bzt vlsvextract

FORCE:
# On FERMI one has to use the front-end compiler (e.g. g++) to compile this tool.
# This target here defines a flag which removes the mpi headers from the code with 
# #ifdef pragmas such that one can compile this tool to be used on the login nodes.
# To ensure this works one also needs to change the compiler at the top of Makefile.fermi*.
not_parallel_tools: vlsvdiff

all: vlasiator tools

# Compile directory:
INSTALL = $(CURDIR)

# Executable:
EXE = vlasiator


# Collect libraries into single variable:
LIBS = ${LIB_BOOST}
LIBS += ${LIB_ZOLTAN}
LIBS += ${LIB_MPI}
LIBS += ${LIB_CUDA}
LIBS += ${LIB_PROFILE}

# Define common dependencies
DEPS_COMMON = common.h definitions.h mpiconversion.h logger.h

# Define dependencies on all project files
DEPS_PROJECTS =	projects/project.h projects/project.cpp \
		projects/projectIsotropicMaxwellian.h projects/projectIsotropicMaxwellian.cpp \
		projects/Alfven/Alfven.h projects/Alfven/Alfven.cpp \
		projects/Diffusion/Diffusion.h projects/Diffusion/Diffusion.cpp \
		projects/Dispersion/Dispersion.h projects/Dispersion/Dispersion.cpp \
		projects/Firehose/Firehose.h projects/Firehose/Firehose.cpp \
		projects/Flowthrough/Flowthrough.h projects/Flowthrough/Flowthrough.cpp \
		projects/Fluctuations/Fluctuations.h projects/Fluctuations/Fluctuations.cpp \
		projects/KHB/KHB.h projects/KHB/KHB.cpp \
		projects/Larmor/Larmor.h projects/Larmor/Larmor.cpp \
		projects/Magnetosphere/Magnetosphere.h projects/Magnetosphere/Magnetosphere.cpp\
		projects/MultiPeak/MultiPeak.h projects/MultiPeak/MultiPeak.cpp \
		projects/Riemann1/Riemann1.h projects/Riemann1/Riemann1.cpp \
		projects/Shock/Shock.h projects/Shock/Shock.cpp \
		projects/test_fp/test_fp.h projects/test_fp/test_fp.cpp \
		projects/test_trans/test_trans.h projects/test_trans/test_trans.cpp \
		projects/verificationLarmor/verificationLarmor.h projects/verificationLarmor/verificationLarmor.cpp
#all objects for vlasiator

OBJS = 	version.o backgroundfield.o ode.o quadr.o dipole.o constantfield.o integratefunction.o \
	datareducer.o datareductionoperator.o \
	donotcompute.o ionosphere.o outflow.o setbyuser.o setmaxwellian.o \
	sysboundary.o sysboundarycondition.o \
	project.o projectIsotropicMaxwellian.o \
	Alfven.o Diffusion.o Dispersion.o Firehose.o Flowthrough.o Fluctuations.o  KHB.o Larmor.o Magnetosphere.o MultiPeak.o Riemann1.o Shock.o test_fp.o test_trans.o verificationLarmor.o\
	grid.o ioread.o iowrite.o vlasiator.o logger.o muxml.o \
	parameters.o readparameters.o spatial_cell.o \
	vlscommon.o vlsvreader2.o vlsvwriter2.o vlasovmover_$(TRANSSOLVER).o $(FIELDSOLVER).o


help:
	@echo ''
	@echo 'make c(lean)             delete all generated files'
	@echo 'make dist                make tar file of the source code'
	@echo 'make ARCH=arch Compile vlasiator '
	@echo '                           ARCH:  Set machine specific Makefile Makefile.arch'

# remove data generated by simulation
d: data
data:
	rm -rf phiprof*txt restart*vlsv grid*vlsv diagnostic.txt logfile.txt

c: clean
clean: data
	rm -rf *.o *~ */*~ */*/*~ ${EXE} vlsv2silo_${FP_PRECISION} vlsvextract_${FP_PRECISION} vlsv2vtk_${FP_PRECISION} vlsvdiff_${FP_PRECISION} vlsv2bzt_${FP_PRECISION} check_projects_compil_logs/ check_projects_cfg_logs/


# Rules for making each object file needed by the executable

version.cpp: FORCE
	./generate_version.sh "${CMP}" "${CXXFLAGS}" "${FLAGS}" "${INC_MPI}" "${INC_DCCRG}" "${INC_ZOLTAN}" "${INC_BOOST}"

version.o: version.cpp 
	 ${CMP} ${CXXFLAGS} ${FLAGS} -c version.cpp

dipole.o: backgroundfield/dipole.cpp backgroundfield/dipole.hpp backgroundfield/fieldfunction.hpp backgroundfield/functions.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/dipole.cpp 

constantfield.o: backgroundfield/constantfield.cpp backgroundfield/constantfield.hpp backgroundfield/fieldfunction.hpp backgroundfield/functions.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/constantfield.cpp 

ode.o: backgroundfield/ode.cpp backgroundfield/ode.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/ode.cpp

quadr.o: backgroundfield/quadr.cpp backgroundfield/quadr.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/quadr.cpp

backgroundfield.o: ${DEPS_COMMON} backgroundfield/backgroundfield.cpp backgroundfield/backgroundfield.h backgroundfield/fieldfunction.hpp backgroundfield/functions.hpp backgroundfield/integratefunction.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/backgroundfield.cpp 

integratefunction.o: ${DEPS_COMMON} backgroundfield/integratefunction.cpp backgroundfield/integratefunction.hpp backgroundfield/functions.hpp  backgroundfield/quadr.cpp backgroundfield/quadr.hpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c backgroundfield/integratefunction.cpp 

datareducer.o: ${DEPS_COMMON} spatial_cell.hpp datareduction/datareducer.h datareduction/datareductionoperator.h datareduction/datareducer.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c datareduction/datareducer.cpp ${INC_MPI} ${INC_BOOST}

datareductionoperator.o:  ${DEPS_COMMON} spatial_cell.hpp datareduction/datareductionoperator.h datareduction/datareductionoperator.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c datareduction/datareductionoperator.cpp ${INC_MPI} ${INC_BOOST}

donotcompute.o: ${DEPS_COMMON} sysboundary/donotcompute.h sysboundary/donotcompute.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c sysboundary/donotcompute.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

ionosphere.o: ${DEPS_COMMON} sysboundary/ionosphere.h sysboundary/ionosphere.cpp backgroundfield/backgroundfield.cpp backgroundfield/backgroundfield.h projects/project.h projects/project.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c sysboundary/ionosphere.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

outflow.o: ${DEPS_COMMON} sysboundary/outflow.h sysboundary/outflow.cpp projects/project.h projects/project.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c sysboundary/outflow.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

setmaxwellian.o: ${DEPS_COMMON} sysboundary/setmaxwellian.h sysboundary/setmaxwellian.cpp sysboundary/setbyuser.h sysboundary/setbyuser.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c sysboundary/setmaxwellian.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

setbyuser.o: ${DEPS_COMMON} sysboundary/setbyuser.h sysboundary/setbyuser.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c sysboundary/setbyuser.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

sysboundary.o: ${DEPS_COMMON} sysboundary/sysboundary.h sysboundary/sysboundary.cpp sysboundary/sysboundarycondition.h sysboundary/sysboundarycondition.cpp sysboundary/donotcompute.h sysboundary/donotcompute.cpp sysboundary/ionosphere.h sysboundary/ionosphere.cpp sysboundary/outflow.h sysboundary/outflow.cpp sysboundary/setmaxwellian.h sysboundary/setmaxwellian.cpp sysboundary/setbyuser.h sysboundary/setbyuser.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c sysboundary/sysboundary.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

sysboundarycondition.o: ${DEPS_COMMON} sysboundary/sysboundarycondition.h sysboundary/sysboundarycondition.cpp sysboundary/donotcompute.h sysboundary/donotcompute.cpp sysboundary/ionosphere.h sysboundary/ionosphere.cpp sysboundary/outflow.h sysboundary/outflow.cpp sysboundary/setmaxwellian.h sysboundary/setmaxwellian.cpp sysboundary/setbyuser.h sysboundary/setbyuser.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c sysboundary/sysboundarycondition.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

Alfven.o: ${DEPS_COMMON} projects/Alfven/Alfven.h projects/Alfven/Alfven.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/Alfven/Alfven.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

Diffusion.o: ${DEPS_COMMON} projects/Diffusion/Diffusion.h projects/Diffusion/Diffusion.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/Diffusion/Diffusion.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

Dispersion.o: ${DEPS_COMMON} projects/Dispersion/Dispersion.h projects/Dispersion/Dispersion.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/Dispersion/Dispersion.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

Firehose.o: ${DEPS_COMMON} projects/Firehose/Firehose.h projects/Firehose/Firehose.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/Firehose/Firehose.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

Flowthrough.o: ${DEPS_COMMON} projects/Flowthrough/Flowthrough.h projects/Flowthrough/Flowthrough.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/Flowthrough/Flowthrough.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

Fluctuations.o: ${DEPS_COMMON} projects/Fluctuations/Fluctuations.h projects/Fluctuations/Fluctuations.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/Fluctuations/Fluctuations.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

KHB.o: ${DEPS_COMMON} projects/KHB/KHB.h projects/KHB/KHB.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/KHB/KHB.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

Larmor.o: ${DEPS_COMMON} projects/Larmor/Larmor.h projects/Larmor/Larmor.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/Larmor/Larmor.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

Magnetosphere.o: ${DEPS_COMMON} projects/Magnetosphere/Magnetosphere.h projects/Magnetosphere/Magnetosphere.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/Magnetosphere/Magnetosphere.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

MultiPeak.o: ${DEPS_COMMON} projects/MultiPeak/MultiPeak.h projects/MultiPeak/MultiPeak.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/MultiPeak/MultiPeak.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

Riemann1.o: ${DEPS_COMMON} projects/Riemann1/Riemann1.h projects/Riemann1/Riemann1.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/Riemann1/Riemann1.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

Shock.o: ${DEPS_COMMON} projects/Shock/Shock.h projects/Shock/Shock.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/Shock/Shock.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

test_fp.o: ${DEPS_COMMON} projects/test_fp/test_fp.h projects/test_fp/test_fp.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/test_fp/test_fp.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

test_trans.o: ${DEPS_COMMON} projects/test_trans/test_trans.h projects/test_trans/test_trans.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/test_trans/test_trans.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

verificationLarmor.o: ${DEPS_COMMON} projects/verificationLarmor/verificationLarmor.h projects/verificationLarmor/verificationLarmor.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/verificationLarmor/verificationLarmor.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

project.o: ${DEPS_COMMON} $(DEPS_PROJECTS)
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/project.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

projectIsotropicMaxwellian.o: ${DEPS_COMMON} $(DEPS_PROJECTS)
	${CMP} ${CXXFLAGS} ${FLAGS} -c projects/projectIsotropicMaxwellian.cpp ${INC_DCCRG} ${INC_ZOLTAN} ${INC_BOOST}

spatial_cell.o: spatial_cell.cpp spatial_cell.hpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c spatial_cell.cpp $(INC_BOOST)

vlasovmover_leveque.o: spatial_cell.hpp transferstencil.h  vlasovsolver/cpu_acc_$(ACCSOLVER).hpp vlasovsolver/cpu_trans_leveque.h vlasovsolver/cpu_lorentz.hpp vlasovsolver/limiters.h  vlasovsolver/limiters_vec4.h vlasovsolver/leveque_common.h vlasovsolver/vlasovmover_leveque.cpp		
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${MATHFLAGS} ${FLAGS}  -DMOVER_VLASOV_ORDER=2  -c vlasovsolver/vlasovmover_leveque.cpp -I$(CURDIR) ${INC_ZOLTAN} ${INC_BOOST} ${INC_DCCRG}  ${INC_PROFILE}  ${INC_VECTORCLASS}

londrillo_delzanna.o: spatial_cell.hpp transferstencil.h   parameters.h common.h fieldsolver/londrillo_delzanna.cpp
	 ${CMP} ${CXXFLAGS} ${FLAGS} -c fieldsolver/londrillo_delzanna.cpp -I$(CURDIR)  ${INC_BOOST} ${INC_DCCRG}  ${INC_PROFILE}  ${INC_ZOLTAN}

vlasiator.o:  ${DEPS_COMMON} readparameters.h parameters.h ${DEPS_PROJECTS} grid.h spatial_cell.hpp vlasiator.cpp
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c vlasiator.cpp ${INC_MPI} ${INC_DCCRG} ${INC_BOOST} ${INC_ZOLTAN} ${INC_PROFILE}

grid.o:  ${DEPS_COMMON} parameters.h ${DEPS_PROJECTS} spatial_cell.hpp grid.cpp grid.h  sysboundary/sysboundary.h
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c grid.cpp ${INC_MPI} ${INC_DCCRG} ${INC_BOOST} ${INC_ZOLTAN} ${INC_PROFILE}

ioread.o:  ${DEPS_COMMON} parameters.h  spatial_cell.hpp ioread.cpp ioread.h  vlsvreader2.cpp
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c ioread.cpp ${INC_MPI} ${INC_DCCRG} ${INC_BOOST} ${INC_ZOLTAN} ${INC_PROFILE}


iowrite.o:  ${DEPS_COMMON} parameters.h  spatial_cell.hpp iowrite.cpp iowrite.h  vlsvwriter2.cpp
	${CMP} ${CXXFLAGS} ${FLAG_OPENMP} ${FLAGS} -c iowrite.cpp ${INC_MPI} ${INC_DCCRG} ${INC_BOOST} ${INC_ZOLTAN} ${INC_PROFILE}

logger.o: logger.h logger.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c logger.cpp ${INC_MPI}

muxml.o: muxml.h muxml.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c muxml.cpp

parameters.o: parameters.h parameters.cpp readparameters.h
	$(CMP) $(CXXFLAGS) $(FLAGS) -c parameters.cpp ${INC_BOOST}

readparameters.o: readparameters.h readparameters.cpp version.h version.cpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c readparameters.cpp ${INC_BOOST}

vlscommon.o:  $(DEPS_COMMON)  vlscommon.h vlscommon.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlscommon.cpp


vlsvreader2.o:  muxml.h muxml.cpp vlscommon.h vlsvreader2.h vlsvreader2.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsvreader2.cpp

vlsvreader2extra.o:  muxml.h muxml.cpp vlscommon.h vlsvreader2.h vlsvreader2.cpp
	${CMP} ${CXXEXTRAFLAGS} ${FLAGS} -c vlsvreader2.cpp

vlsvwriter2.o: mpiconversion.h muxml.h muxml.cpp vlscommon.h vlsvwriter2.h vlsvwriter2.cpp 
	${CMP} ${CXXFLAGS} ${FLAGS} -c vlsvwriter2.cpp ${INC_MPI}



# Make executable
vlasiator: $(OBJS)
	$(LNK) ${LDFLAGS} -o ${EXE} $(OBJS) $(LIBS)


#/// TOOLS section/////

#common reader filter
DEPS_VLSVREADER = muxml.h muxml.cpp vlscommon.h vlsvreader2.h vlsvreader2.cpp 

#common reader objects for tools
OBJS_VLSVREADER = muxml.o vlscommon.o vlsvreader2.o
OBJS_VLSVREADEREXTRA = muxml.o vlscommon.o vlsvreader2extra.o


vlsvextract: ${DEPS_VLSVREADER} tools/vlsvextract.cpp ${OBJS_VLSVREADER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsvextract.cpp ${INC_BOOST} ${INC_DCCRG} ${INC_SILO} ${INC_EIGEN} -I$(CURDIR) 
	${LNK} -o vlsvextract_${FP_PRECISION} vlsvextract.o ${OBJS_VLSVREADER} ${LIB_BOOST} ${LIB_DCCRG} ${LIB_SILO} ${LDFLAGS}


vlsv2vtk: ${DEPS_VLSVREADER} ${OBJS_VLSVREADER} tools/vlsv2vtk.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsv2vtk.cpp ${INC_BOOST} -I$(CURDIR) 
	${LNK} -o vlsv2vtk_${FP_PRECISION} vlsv2vtk.o ${OBJS_VLSVREADER} ${INC_BOOST} ${LDFLAGS}


vlsv2silo: ${DEPS_VLSVREADER} ${OBJS_VLSVREADER} tools/vlsv2silo.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsv2silo.cpp ${INC_SILO} -I$(CURDIR) 
	${LNK} -o vlsv2silo_${FP_PRECISION} vlsv2silo.o ${OBJS_VLSVREADER} ${LIB_SILO} ${LDFLAGS}

vlsv2bzt: ${DEPS_VLSVREADER} ${OBJS_VLSVREADER} tools/vlsv2bzt.cpp
	${CMP} ${CXXFLAGS} ${FLAGS} -c tools/vlsv2bzt.cpp -I$(CURDIR) 
	${LNK} -o vlsv2bzt_${FP_PRECISION} vlsv2bzt.o ${OBJS_VLSVREADER} ${LDFLAGS}

vlsvdiff: ${DEPS_VLSVREADER} ${OBJS_VLSVREADEREXTRA} tools/vlsvdiff.cpp
	${CMP} ${CXXEXTRAFLAGS} ${FLAGS} -c tools/vlsvdiff.cpp -I$(CURDIR)
	${LNK} -o vlsvdiff_${FP_PRECISION} vlsvdiff.o ${OBJS_VLSVREADER} ${LDFLAGS}


# DO NOT DELETE
