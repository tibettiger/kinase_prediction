from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='KBPRNA',
    version='0.0.1',
    packages=['KBPRNA'],
    desription='Machine Learning model of predicting kinase activity',
    license='MIT',
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'scanpy',
        'anndata',
        'matplotlib',
        'xgboost',
        'sklearn',
        'random'
    ]
)
