from distutils.core import setup
from distutils.extension import Extension

module1 = Extension('msAlignForPy',
                    #include_dirs = ['/usr/local/Cellar/boost/1.65.1/include'],
                    libraries = ['boost_python3'],
                    #library_dirs = ['/usr/local/Cellar/boost/1.65.1/lib'],
                    extra_compile_args=['-std=c++11'],
                    sources = ['msAlignForPy.cpp', 'ActiveSequence.cpp', 'Heap.cpp'])

setup (name = 'msAlign',
       version = '1.0',
       description = 'This is package for ms alignment',
       author = 'Mario Marchand',
       ext_modules = [module1])