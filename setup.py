from setuptools import setup, find_packages

setup(
    name='assets',
    version='0.0.1',
    packages=find_packages(include=['assets']),
    install_requires=[
        "deepdiff == 5.6.0",
        "natsort == 7.1.1",
        "numpy == 1.21.2",
        "ordered-set == 4.0.2",
        "pandas == 1.3.3",
        "pysam == 0.17.0",
        "pyvcf == 0.6.8",
        "pyyaml == 5.4.1",
    ],
    setup_requires=['pytest-reqs==0.2.1', 'pytest==6.2.5'],
    tests_require=['pytest-reqs==0.2.1', 'pytest==6.2.5'],
)
