import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SageNP",
    version="0.1",
    author="Tolga Birkandan, Emir Baysazan, Pelin Ozturk",
    author_email="birkandant@itu.edu.tr",
    description="Newman-Penrose calculations for SageMath",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tbirkandan/SageNP",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    keywords="SageMath, Newman-Penrose formalism, Petrov classification"
)