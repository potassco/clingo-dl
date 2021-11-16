import sys
import site
from os.path import dirname, abspath
from textwrap import dedent
from skbuild import setup
import clingo


if not site.ENABLE_USER_SITE and "--user" in sys.argv[1:]:
    site.ENABLE_USER_SITE = True

clingopath = abspath(dirname(clingo.__file__))

setup(
    version = '1.3.0',
    name = 'clingo-dl',
    description = 'CFFI-based bindings to the clingo-dl solver.',
    long_description = dedent('''\
        This package allows for adding the clingo-dl propagator as a
        theory to clingo.

        It can also be used as a clingo-dl solver running:

            python -m clingodl CLINGODL_ARGUMENTS
        '''),
    long_description_content_type='text/markdown',
    author = 'Roland Kaminski',
    author_email = 'kaminski@cs.uni-potsdam.de',
    license = 'MIT',
    url = 'https://github.com/potassco/clingo-dl',
    install_requires=[ 'cffi', 'clingo' ],
    cmake_args=[ '-DCLINGODL_MANAGE_RPATH=OFF',
                 '-DPYCLINGODL_ENABLE=pip',
                 '-DPYCLINGODL_INSTALL_DIR=libpyclingo-dl',
                 f'-DPYCLINGODL_PIP_PATH={clingopath}' ],
    packages=[ 'clingodl' ],
    package_data={ 'clingodl': [ 'py.typed', 'import__clingo-dl.lib', 'clingo-dl.h' ] },
    package_dir={ '': 'libpyclingo-dl' },
    python_requires=">=3.6"
)
