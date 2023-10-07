from setuptools import find_packages, setup

setup(
    name='lmbio', # 应用名
    packages=find_packages(),
    version='1.0.0.1', # 版本号
    description='lmbio analysis tools',
    author='lujiawei',
    license='MIT',
    install_requires=["numpy","pandas",
                      "matplotlib",
                      "rpy2","argparse"], # 依赖列表
    python_requires='>=3.7',
    package_data={'lmbio': ['spatialmetabolism/spatial/spatialMeta/R/*',]}
)
