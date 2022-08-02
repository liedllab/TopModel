from setuptools import setup, find_packages

setup(
        name='pdbear', 
        packages=find_packages(),
        include_package_data=True,
        install_requires=[
            'Click',
            ],
        entry_points={
            'console_scripts' : [
                'pdbear = pdbear.app:main',
                ]
            },
        )
