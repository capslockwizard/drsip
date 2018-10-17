import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="drsip",
    version="0.5",
    author="Justin Chan",
    author_email="capslockwizard@gmail.com",
    description="DRSIP docking package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    install_requires=['zdock-parser','numba', 'numpy', 'scipy', 'pandas', 'u-msgpack-python', 'MDAnalysis', 'biopython', 'drsip-common', 'docking-eval'],
    classifiers=[
        "Environment :: Console",
        "Environment :: Plugins",
        "Intended Audience :: Science/Research  ",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
