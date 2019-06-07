import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="IR_JaGeo",
    version="1.0.4",
    author="Janine George",
    author_email="janine.george@uclouvain.be",
    description="Package to calculate infrared ntensities with the dipole approximation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JaGeo/IR",
    packages=setuptools.find_packages(),
    setup_requires=['setuptools>=18.0'],
    install_requires=["phonopy>=2.1.1"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
)
