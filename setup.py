from setuptools import setup, find_packages

with open("./assets/requirements.txt") as requirements_file:
    requirements = requirements_file.readlines()

setup(
    name='assets',
    version='0.0.1',
    packages=find_packages(include=['assets']),
    install_requires=requirements,
    setup_requires=['pytest-reqs>=0.2.1', 'pytest>=6.2.5', 'pytest-flake8>=1.0.7', 'flake8<=4.0.0'],
    tests_require=['pytest-reqs>=0.2.1', 'pytest>=6.2.5', 'pytest-flake8>=1.0.7', 'flake8<=4.0.0'],
)
