import os

threads = 2

common_src = ['components.cc']
include_path = ['/opt/local/include']
lib_path = ['/lib','/usr/lib','/usr/lib64','/usr/local/lib','/usr/local/lib64','/opt/local/lib']

env = Environment(LIBPATH = lib_path,
                  CCFLAGS = ['-O2','-Wall'],
                  LIBS=['m','blitz','sndfile','fftw3','lapack',])
env.Append(CPPPATH = include_path)
if threads > 1:
    env.Append(CCFLAGS = '-DTHREADS=%d' % threads,
               LIBS = ['fftw3_threads','pthread'])

env.Program('znote_label',['znote_label.cc','mtm.cc','dpss.c'] + common_src)
env.Program('znote_extract',['znote_extract.cc','spect.cc'] + common_src)
