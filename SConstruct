
import os

# source files and libraries
src_files = ['dpss.c','mtm.cc','components.cc']
hdr_files = ['blitz_io.hh','mtm.hh']
include_path = ['/opt/local/include']
lib_path = ['/lib','/usr/lib','/usr/lib64','/usr/local/lib','/usr/local/lib64','/opt/local/lib']

env = Environment(LIBPATH = lib_path,
                  CCFLAGS = ['-O2','-Wall'],
                  LIBS=['m','blitz','sndfile','fftw3','lapack',])
env.Append(CPPPATH = include_path)

env.Program('znote_label',src_files + ['znote_label.cc'])
