from setuptools import setup, find_packages
import nimporter

setup(
    name="librtd",
    version="0.0.1",
    author="Benjamin D. Lee",
    author_email="benjamindlee@me.com",
    packages=find_packages(),
    description="Generalized k-mer return time distribution calculation",
    # long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    install_requires=["docopt"],
    url="https://github.com/Lab41/librtd",
    ext_modules=nimporter.build_nim_extensions(danger=True),
    entry_points={"console_scripts": ["rtd = cli_wrapper:cli_wrapper"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
