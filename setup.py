from setuptools import find_packages, setup
import nimporter

setup(
    name='librtd',
    version='0.1dev',
    packages=find_packages(),
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description="hi",
    install_requires=["nimporter"],
    ext_modules=nimporter.build_nim_extensions(),
)
