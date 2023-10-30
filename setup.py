from setuptools import Extension, find_packages, setup
from Cython.Build import cythonize
import numpy

extentions = [
    Extension("src.vcsolver.utils", ["src/vcsolver/utils.pyx"], include_dirs=[numpy.get_include()]),
    Extension("src.vcsolver.packing_constraints", ["src/vcsolver/packing_constraints.pyx"], include_dirs=[numpy.get_include()]),
    Extension("src.vcsolver.clique_cover", ["src/vcsolver/clique_cover.pyx"], include_dirs=[numpy.get_include()]),
]

setup(
    name="vcsolver",
    ext_modules=cythonize(extentions),
    # include_dirs=[numpy.get_include()],
    packages=find_packages("src/vcsolver"),
    zip_safe=False,
)
