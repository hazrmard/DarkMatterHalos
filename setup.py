from setuptools import setup, find_packages

setup(
    name="halos",
    version="0.7.3",
    author="Ibrahim Ahmed",
    description="Find half mass radii of halos.",
    packages=find_packages(),
    install_requires=['numpy', 'matplotlib']
)
