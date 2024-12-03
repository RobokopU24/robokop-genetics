import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="robokop-genetics",
    version="0.6.0",
    author="Evan Morris",
    author_email="evandietzmorris@gmail.com",
    maintainer="Evan Morris",
    maintainer_email="evandietzmorris@gmail.com",
    description="A package for Robokop genetics tools and services.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/RobokopU24/robokop-genetics",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    license="CC-0",
    install_requires=[
        "requests>=2.32.3",
        "redis>=5.0.4"
    ]
)
