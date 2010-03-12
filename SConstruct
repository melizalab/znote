import os

threads = 1
if hasattr(os,'uname'):
    system = os.uname()[0]
else:
    system = 'Windows'

common_src = ['components.cc']
lib_path = ['/lib','/usr/lib','/usr/lib64','/usr/local/lib','/usr/local/lib64']

env = Environment(LIBPATH = lib_path,
                  CCFLAGS = ['-O2','-Wall'],
                  LIBS=['m','blitz','sndfile','fftw3'])
if threads > 1:
    env.Append(CCFLAGS = '-DTHREADS=%d' % threads,
               LIBS = ['fftw3_threads','pthread'])

# system-specific settings
if system=='Darwin':
    env.Append(LIBS = ['lapack'],
               CPPPATH = ['/opt/local/include'],
               LIBPATH = ['/opt/local/lib'])
elif system=='Linux':
    env.Append(LIBS = ['atlas','cblas','f77blas','lapack'])
elif system=='Windows':
    env = Environment(platform = Platform('win32'), tools=['mingw'])
    env.Append(CCFLAGS = ['-Wall'],
               CPPPATH = ['./blitz-0.9','./libsndfile/include','./fftw/include'],
               LIBPATH = ['./blitz-0.9/lib','./libsndfile/lib','./fftw/lib','./lapack/'],
               LIBS=['m','blitz','sndfile','fftw3','lapack','blas','g2c'])

env.Program('znote_label',['znote_label.cc','mtm.cc','dpss.c'] + common_src)
env.Program('znote_extract',['znote_extract.cc','spect.cc'] + common_src)
