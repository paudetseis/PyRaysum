import os.path
import re
from numpy.distutils.core import setup, Extension
from numpy.distutils.system_info import get_info

def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname, encoding='utf-8') as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


# ext = [Extension(name='pyraysum.raysum_f',
#                  sources=['src/raysum.f90', 'src/raysum_sub.f90'],
#                  libraries=['lapack'],
#                  library_dirs=get_info('lapack_opt', 1).get('library_dirs'))]

ext = [Extension(name='pyraysum.raysum_f',
                 sources=['src/raysum.f90'])]

setup(
    name='pyraysum',
    version=find_version('pyraysum', '__init__.py'),
    description='Python version of Fortran code Raysum',
    author='Andrew Frederiksen, Pascal Audet',
    maintainer='Pascal Audet',
    author_email='pascal.audet@uottawa.ca',
    # url='https://github.com/paudetseis/Telewavesim',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Fortran',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8'],
    install_requires=['numpy>=1.15', 'obspy>=1.0.0', 'matplotlib'],
    python_requires='>=3.7',
    tests_require=['pytest'],
    ext_modules=ext,
    packages=['pyraysum', 'pyraysum.tests'],
    package_data={
        'pyraysum': [
            # 'examples/*.ipynb',
            'examples/models/*.txt']
            # 'examples/Notebooks/*.ipynb']
    }
)
