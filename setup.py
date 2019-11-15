from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

requirements = [
    "numpy>=1.14.0",
    "pandas>=0.21.0",
    "tensorflow>=1.13.1",
    "h5py>=2.9.0"
]

setup(
    name="disc",
    version="0.0.0.4",
    author="iyhaoo",
    author_email="904469382@qq.com",
    description="An accurate and scalable imputation algorithm based on semi-supervised deep learning for single-cell transcriptome",
    install_requires=requirements,
    long_description=readme,
    url="https://github.com/iyhaoo/DISC",
    packages=find_packages(),
    classifiers=[
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires='>=3.6',
)




