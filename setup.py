from setuptools import setup, find_packages

requirements = [
    "numpy>=1.14.0",
    "pandas>=0.21.0",
    "tensorflow>=1.13.1,<2.0.0",
    "h5py>=2.9.0"
    "matplotlib>=3.0.0"
]

setup(
    name="disc",
    version="1.1.6",
    author="Zhongshan Ophthalmic Centre (ZOC), Sun Yat-sen University (SYSU)",
    author_email="904469382@qq.com",
    description="A highly scalable and accurate inference of gene expression and structure for single-cell transcriptomes using semi-supervised deep learning.",
    install_requires=requirements,
    url="https://github.com/iyhaoo/DISC",
    packages=find_packages(),
    license='Apache License 2.0',
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: Apache Software License",
        "Topic :: Scientific/Engineering :: Artificial Intelligence"
    ],
    entry_points={
        'console_scripts': [
            'disc = disc.__main__:main'
        ]},
    python_requires='>=3.6',
)




