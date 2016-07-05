from Cython.Build import cythonize
from setuptools import setup, Extension
import numpy

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

setup(
    version='0.1',
    author='...',
    license='...',
    description='...',
    long_description='...',
    url='...',
    ext_modules=cythonize(Extension(
           "interpDetect",
           sources=["interpDetect.pyx", "SpkDslowFilter.cpp"],
           language="c++",
           extra_compile_args=['-std=c++11', '-Wunused-function'],
           )),
    include_dirs=[numpy.get_include()]
)
