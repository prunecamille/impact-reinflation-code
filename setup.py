"""
Setup script for the impact-reinflation package.
"""

from setuptools import setup, find_packages

setup(
    name="impact-reinflation",
    version="1.0.0",
    description="Model for planetary atmosphere evolution through impacts, escape, and outgassing",
    author="Prune Camille August",
    author_email="pcaugust@seas.harvard.edu",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "matplotlib>=3.4.0",
    ],
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    entry_points={
        'console_scripts': [
            'impact-run=run_single:main',
            'impact-montecarlo=run_montecarlo:main',
            'impact-plot=plot_results:main',
        ],
    },
)
