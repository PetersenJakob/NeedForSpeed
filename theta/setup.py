from setuptools import setup, Extension
import pybind11

# cpp_args = ['-std=c++11', '-stdlib=libc++', '-mmacosx-version-min=10.7']
cpp_args = ['-std=c++11', '-stdlib=libc++']

sfc_module = Extension(
    'nfs',
    sources=['module.cpp'],
    include_dirs=[pybind11.get_include()],
    language='c++',
    extra_compile_args=cpp_args,
    )

setup(
    name='nfs',
    version='1.0',
    description='C++ accelerator module for FinPy',
    ext_modules=[sfc_module],
)