import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PCIst",
    version="0.1.15",
    author="Renzo Comolatti",
    author_email="renzo.com@gmail.com",
    description="Short library to calculate PCIst.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/renzocom/PCIst",
    packages=setuptools.find_packages(include=['PCIst']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL License",
        "Operating System :: OS Independent",
    ],
)
