import setuptools
import os.path
import re
from numpy.distutils.core import setup, Extension
from numpy.distutils.system_info import get_info

ext = Extension(name='fraysum',
                sources = ['pyraysum/src/buildmodel.f',
                           'pyraysum/src/eigenvec.f',
                           'pyraysum/src/eispack-cg.f',
                           'pyraysum/src/matrixops.f',
                           'pyraysum/src/phaselist.f',
                           'pyraysum/src/raysum.f',
                           'pyraysum/src/readwrite.f',
                           'pyraysum/src/call-seis-spread.f',
                           'pyraysum/src/trace.f'],
                extra_compile_args=['-O3'])

def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname, encoding='utf-8') as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name='pyraysum',
    version=find_version('pyraysum', '__init__.py'),
    description='Python wrapper around Fortran code Raysum',
    author='Andrew Frederiksen, Pascal Audet',
    maintainer='Pascal Audet',
    author_email='pascal.audet@uottawa.ca',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Fortran',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8'],
    install_requires=['numpy>=1.15', 'obspy>=1.0.0', 'matplotlib'],
    python_requires='>=3.7',
    ext_modules = [ext],
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={
        'pyraysum': [
            'examples/models/*.txt',
            'examples/Notebooks/*.ipynb']},
)
