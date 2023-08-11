from setuptools import Extension, find_packages, setup
from Cython.Build import cythonize
import numpy

extentions = [
    Extension("solver.utils", ["solver/utils.pyx"], include_dirs=[numpy.get_include()]),
    Extension("solver.packing_constraints", ["solver/packing_constraints.pyx"], include_dirs=[numpy.get_include()]),
    Extension("solver.clique_cover", ["solver/clique_cover.pyx"], include_dirs=[numpy.get_include()]),
]

setup(
    name="VCSolver",
    ext_modules=cythonize(extentions),
    # include_dirs=[numpy.get_include()],
    packages=find_packages(),
    zip_safe=False,
)
