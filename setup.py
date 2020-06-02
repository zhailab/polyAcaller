from setuptools import setup, find_packages

def get_version(string):
    """ Parse the version number variable __version__ from a script. """
    import re
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    version_str = re.search(version_re, string, re.M).group(1)
    return version_str

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')


install_requires_py = ["numpy >=1.16",
                       "matplotlib ==3.1.1",
                       "pandas",
                       "ont_fast5_api"]

setup(
    name='polyAcaller',
    version=get_version(open('polyacaller/__init__.py').read()),
    author='Jinbu Jia',
    author_email='jiajb@sustech.edu.cn',
    packages=find_packages(exclude=['tests']),
    scripts=['bin/polyAcaller'],
    include_package_data=True,
    package_dir={'polyacaller': 'polyacaller'},
    url='https://github.com/zhailab/polyAcaller',
    license='LICENSE.txt',
    description='Command-line tool and api for find polyA from nanopore basecalled fast5 file',
    long_description=open('Readme.md').read(),
    long_description_content_type="text/markdown",
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=install_requires_py,
    zip_safe=False,
    python_requires='>=3.6.*, <4'
)
