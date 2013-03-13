import os

if hasattr(os,'uname'):
    system = os.uname()[0]
else:
    system = 'Windows'

debug = ARGUMENTS.get('debug',0)
threads = int(ARGUMENTS.get('thread',1))

common_src = ['components.cc']

env = Environment(CCFLAGS = ['-Wall'],
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

if int(debug):
    env.Append(CCFLAGS=['-g2', '-DDEBUG=1'])
else:
    env.Append(CCFLAGS=['-O2'])


env.Program('znote_label',['znote_label.cc','mtm.cc','dpss.c'] + common_src)
env.Program('znote_extract',['znote_extract.cc','spect.cc'] + common_src)
