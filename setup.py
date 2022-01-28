from setuptools import setup,find_packages

setup(
    name='lobsterpy',
    version='0.1.0',
    description='Package for autmatic bonding analysis with Lobster/VASP',
    url='https://github.com/jageo/lobsterpy',
    author='Janine George',
    author_email='janine.george@bam.de',
    license='BSD 3-clause',
    packages=find_packages(),
    install_requires=['pymatgen>=2022.1.9',
                      'numpy',
                      'typing'],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],

    entry_points={
        "console_scripts": [
            "lobsterpy = lobsterpy.cli:main",

        ]
    },
)
