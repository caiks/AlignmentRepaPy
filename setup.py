from distutils.core import setup, Extension

setup (name = 'AlignmentForeignPy',
    version = '1.0',
    ext_modules = [Extension('AlignmentForeignPy', ['AlignmentForeignPy.c'])])

