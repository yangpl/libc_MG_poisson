SXAMGLIB        = sxamg
LIB_SUFFIX      = .a

CPP	            = /usr/bin/gcc -E
CPPFLAGS        = -I/home/pyang/Documents/libc_EM_iterative/sxamg-master/include  -I/usr/local/include

CC              = /usr/bin/gcc 
CFLAGS          = -g -O2 -Wall -Wextra -Wmissing-declarations -Wmissing-prototypes -Wimplicit-function-declaration -Wunused -Wuninitialized -Wvla -Wformat -Wformat-security -Wmain -Wsequence-point -Wshadow -Wbad-function-cast -Wcomments -Wdeclaration-after-statement -Waddress -Wpointer-sign -Wempty-body -Wnested-externs -Wold-style-declaration -Wold-style-definition -Wmissing-field-initializers -Wmaybe-uninitialized -Wno-return-local-addr -Wreturn-type -Wsequence-point -Wundef -Wunknown-pragmas -Wunused-label

CXX             = @CXX@
CXXFLAGS        = @CXXFLAGS@

FC              = @FC@
FCFLAGS         = @FCFLAGS@
FCLIBS          = @FCLIBS@

LINKER          = /usr/bin/gcc 
LDFLAGS         = -L/home/pyang/Documents/libc_EM_iterative/sxamg-master/src  -L/usr/local/lib
LIBS            = -lsxamg -L/usr/local/lib  -lm -lm

AR              = ar
RANLIB          = ranlib

prefix		= /usr/local/sxamg
libdir      = /usr/local/sxamg/lib
includedir  = /usr/local/sxamg/include
RPATH_FLAG	= -Wl,-rpath,

%.o : %.c
	${CC} ${CFLAGS} ${CPPFLAGS} -c $*.c
	@echo

%.o : %.cxx
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -c $*.cxx
	@echo

%.o : %.f
	${FC} ${FCFLAGS} -c $*.f
	@echo

% : %.o
	${LINKER} ${LDFLAGS} -o $@ $*.o ${LIBS}
	@echo

