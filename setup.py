from Cython.Build import cythonize
from setuptools import setup, Extension
import numpy

setup(
    version='1.0',
    author='Oliver Muthmann, Matthias H Hennig, Albert Puente Encinas',
    license='GPL3',
    description='Efficient spike detection for extracellular recordings.',
    url='http://github.com/',
    ext_modules=cythonize(Extension(
           'interpDetect',
           sources=['interpDetect.pyx', 'SpkDslowFilter.cpp'],
           language='c++',
           extra_compile_args=['-DBUILD_PLATFORM_SPIR',
                               '-I/afs/inf.ed.ac.uk/user/s15/s1575609/ComputeCpp-16.05-Linux//include',
                               '-I./temp_files',
                               '-D_GLIBCXX_USE_CXX11_ABI=0',
                               '-include', './temp_files/SpkDslowFilter.cpp.sycl',
                               '-std=c++11', 
                               '-pthread'],
           extra_link_args=['-rdynamic',                            
                            '-L/disk/scratch/apuente/drivers/AMDAPPSDK-3.0/lib/x86_64/','-lOpenCL',
                            '-L/afs/inf.ed.ac.uk/user/s15/s1575609/ComputeCpp-16.05-Linux//lib/', '-lSYCL']
           )),
    include_dirs=[numpy.get_include()], requires=['h5py']
)
