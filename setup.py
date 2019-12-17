from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

requirements = [
    "numpy>=1.14.0",
    "pandas>=0.21.0",
    "tensorflow>=1.13.1,<2.0.0",
    "h5py>=2.9.0"
]

setup(
    name="disc",
    version="0.0.0.29",
    author="Zhongshan Ophthalmic Centre (ZOC), Sun Yat-sen University (SYSU)",
    author_email="904469382@qq.com",
    description="A deep learning method with semi-supervised learning, aiming to accurately recover true gene expression and expression structures for large scRNA-seq data.",
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
        "Topic :: Scientific/Engineering :: Artificial Intelligence"
    ],
    entry_points={
        'console_scripts': [
            'disc = disc.__main__:main'
        ]},
    python_requires='>=3.6',
)




