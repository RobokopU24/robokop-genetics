import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="robokop-genetics",
    version="0.1.0",
    author="Evan Morris",
    author_email="evandietzmorris@gmail.com",
    description="A package for Robokop genetics tools and services.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ObesityHub/robokop-genetics",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=["requests", "redis"]
)