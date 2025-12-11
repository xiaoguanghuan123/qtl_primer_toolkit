#!/usr/bin/env python3
"""
打包配置文件 - 简化版本，避免依赖问题
"""
from setuptools import setup, find_packages

# 读取版本号
with open('qtlprimer/__init__.py', 'r') as f:
    for line in f:
        if line.startswith('__version__'):
            version = line.split('=')[1].strip().strip("'\"")
            break
    else:
        version = '1.0.0'

setup(
    name='qtl-primer-toolkit',
    version=version,
    author='Your Name',
    author_email='your.email@unl.edu',
    description='QTL精细定位引物设计工具包',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/qtl-primer-toolkit',
    packages=find_packages(),
    install_requires=[
        'primer3-py>=2.0.0',
        'biopython>=1.79',
        'pandas>=1.3.0',
        'matplotlib>=3.4.0',
        'numpy>=1.21.0',
        'jinja2>=3.0.0',
        'openpyxl>=3.0.0',
        'pyyaml>=6.0',
    ],
    # 可选依赖
    extras_require={
        'full': [
            'pyfaidx>=0.6.0',
            'seaborn>=0.11.0',
            'plotly>=5.0.0',
        ],
        'dev': [
            'pytest>=6.0.0',
            'flake8>=3.9.0',
        ],
    },
    entry_points={
        'console_scripts': [
            'qtlprimer=qtlprimer.cli:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    python_requires='>=3.8',
)
