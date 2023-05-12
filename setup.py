from setuptools import setup

setup(
    name = "pyrnet",
    version = "0.0.1alpha",
    author = "Jonas Witthuhn",
    author_email= "witthuhn@tropos.de",
    license= "OSI Approved :: GNU General Public License v3 (GPLv3)",
    packages= ["pyrnet"],
    package_dir={"":"src"},
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'pyrnet = pyrnet.click:cli',
        ],
    },
    python_requires= ">=3.7, <3.10",
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
        "trosat-base @ git+https://github.com/hdeneke/trosat-base.git#egg=trosat-base",
    ]

)