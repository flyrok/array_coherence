"""

See:
https://github.com/flyrok/array_coherence
"""

from setuptools import setup, find_packages
from pathlib import Path

here = Path(__file__).resolve().parent

# Get the long description from the README file
readme=here / 'README.md'
with open(readme, encoding='utf-8') as f:
    long_description = f.read()

PROJECT_NAME="array_coherence"
exec(open(here / "array_coherence/__init__.py").read())

VERSION=__version__
DESCRIPTION="Search and pull station metadata from FDSN server"
URL="https://github.com/flyrok/array_coherence"
AUTHOR="A Ferris"
EMAIL="aferris@flyrok.org"
CLASSIFIERS=['Development Status :: 3 - Alpha',
    'Intended Audience :: Seismic Researchers',
    'Topic :: Obspy/FDSN :: Helper Scripts',
    'License :: OSI Approved :: GPL-3 License',
     'Programming Language :: Python :: 3']
KEYWORDS="seismology obspy earthquakes fdsn"     

INSTALL_REQUIRES = [
    'numpy>=1.6.1',
    'matplotlib>=1.1.0',
    'obspy',
    ]

setup(
    name=PROJECT_NAME,  # Required
    version=VERSION,  # Required
    description=DESCRIPTION,  # Optional
    long_description=long_description,  # Optional
    long_description_content_type='text/markdown',  # Optional (see note above)
    url=URL,  # Optional
    author=AUTHOR,  # Optional
    author_email=EMAIL,  # Optional
    classifiers=CLASSIFIERS ,
    keywords=KEYWORDS,  # Optional
    python_requires='>=3.6',
    include_package_data=True,
    packages=['array_coherence'],
    entry_points={  # Optional
        'console_scripts': [
            'array_coherence=array_coherence.__main__:main',
        ],
    },
    install_requires=INSTALL_REQUIRES,
    extras_require={},
    package_data={  
    },
    project_urls={  # Optional
        'Source': URL,
    },
)

