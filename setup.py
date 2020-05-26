from setuptools import setup
import nimporter

setup(
    name="librtd",
    version="0.1dev",
    packages=["librtd"],
    license="Creative Commons Attribution-Noncommercial-Share Alike license",
    long_description="hi",
    install_requires=["nimporter", "docopt"],
    ext_modules=nimporter.build_nim_extensions(),
    entry_points={"console_scripts": ["librtd = librtd.cli_wrapper:cli_wrapper"]},
)
