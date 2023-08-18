from setuptools import setup

setup(
    name='TrajMaker',
    version='0.1',
    packages=[''],
    url='https://github.com/pablogalaviz/trajmaker.git',
    license='GPLv3',
    author='Pablo Galaviz',
    install_requires=[
        'numpy',
        'PyYAML',
        'pandas',
        'ase',
        'jsonschema'
      ],
    scripts=['scripts/trajmaker.py'],
    author_email='galavizp@ansto.gov.au',
    description='Takes an input structure in POSCAR format and creates an XDATCAR trajectory with a prescribed '
                'oscillatory dynamics.'
)
