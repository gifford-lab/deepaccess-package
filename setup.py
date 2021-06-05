import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="deepaccess", 
    version="0.1.0",
    author="Jennifer Hammelman",
    author_email="jhammelm@mit.edu",
    description="A package for training and interpreting an "
    +"ensemble of neural networks for chromatin accessibility",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    project_urls={
        "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"deepaccess": "deepaccess"},
    python_requires=">=3.6",
    entry_points={
        "console_scripts": ["deepaccess=deepaccess.deepaccess:main"]
    },
    include_package_data=True,
    package_data={'': ['train/homer_matrix.npy'],
                  '':['interpret/data/*']},
)
