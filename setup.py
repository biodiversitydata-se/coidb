from setuptools import setup, find_namespace_packages

with open("README.md") as f:
    long_description = f.read()

setup(
    name="coidb",
    version="0.2",
    author="John Sundh",
    url="https://github.com/NBISweden/coidb/",
    description="Workflow for downloading and formatting COI database",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",

    python_requires=">=3.6",
    install_requires=[
        "snakemake",
        "biopython",
        "tqdm",
        "pandas",
        "vsearch",
        "importlib_resources"
    ],
    package_dir={"": "src"},
    packages=find_namespace_packages("src"),
    scripts=['src/coidb/scripts/cluster_bold.py'],
    package_data={
        "coidb": [
            "Snakefile",
            "config.schema.yaml",
            "config.yaml",
        ]
    },
    entry_points={
        "console_scripts": [
            "coidb = coidb.__main__:main"
        ]
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ]
)