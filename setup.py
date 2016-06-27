from distutils.core import setup

setup(
    name='Pyspeqtr',
    version='0.2.0',
    author='Nathan Brewer',
    author_email='brewer.nathant@gmail.com',
    packages=['Pyspeqtr'],
    url=['https://github.com/ntbrewer/Pyspeqtr'],
    scripts=['bin/py_grow_decay.py',
             'bin/spectrum_fitter.py'],
    license='LICENSE.txt',
    description='Useful spectroscopic tools and display in qt',
    long_description=open('README.txt').read(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    requires=['matplotlib', 'numpy', 'lmfit','pyqtgraph'],
)
