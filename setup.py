from setuptools import setup
import versioneer

setup(
    name = "pyrnet",
    version = versioneer.get_version(),
    cmdclass= versioneer.get_cmdclass(),
    author = "Jonas Witthuhn",
    author_email= "witthuhn@tropos.de",
    license= "OSI Approved :: GNU General Public License v3 (GPLv3)",
    packages= ["pyrnet"],
    package_dir={"":"src"},
    package_data={"pyrnet": ["share/*.json"]},
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'pyrnet = pyrnet.click:cli',
        ],
    },
    python_requires=">=3.10",
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "xarray",
        "netcdf4",
        "pyproj",
        "pip",
        "jstyleson",
        "Click",
        "toolz",
        "addict",
        "openpyxl",
        "xlrd",
        "parse>=1.20",
        "limepy",
        "trosat-base @ git+https://github.com/hdeneke/trosat-base.git#egg=trosat-base",
    ],
    extras_require={
        "nbs": [
            "jupyter",
            "nbdev",
            "nbformat",
            "cfchecker",
            "udunits2>=2.2.25"
        ],
        "docs": [
            "sphinx",
            "myst-parser",
            "myst-nb",
        ],
    },



)