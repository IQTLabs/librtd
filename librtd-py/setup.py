from setuptools import setup
import nimporter

setup(
    name="librtd",
    version="0.0.2.2",
    author="Benjamin D. Lee",
    author_email="benjamindlee@me.com",
    packages=["librtd"],
    description="Generalized k-mer return time distribution calculation",
    long_description="For the README, please look [here](https://github.com/IQTLabs/librtd)",
    long_description_content_type="text/markdown",
    install_requires=["docopt"],
    url="https://github.com/IQTLabs/librtd",
    ext_modules=nimporter.build_nim_extensions(danger=True),
    entry_points={"console_scripts": ["rtd = librtd.librtd:_cli_wrapper"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
