import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pypm3D",
    version="0.0.1",
    author="Lucas Ribeiro Alves",
    author_email="lucas.ra.mast@gmail.com",
    description="Python package.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lucasralves/pypm3D",
    project_urls={
        "Bug Tracker": "https://github.com/lucasralves/pypm3D/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)