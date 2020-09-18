import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

setuptools.setup(
	name="AIRgb",
	version="0.0.1",
	author="Tilman Mehl",
	author_email="tilmanmehl@zedat.fu-berlin.de",
    description="A package to automatically access the inverted repeats of plastid genomes stored on NCBI GenBank",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/michaelgruenstaeudl/AIRgb',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: UNIX/Linux",
    ],
    python_requires='>=3.6',
)
