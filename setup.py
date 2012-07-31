from distutils.core import setup

setup(
        name='srtools',
        version='0.1.0',
        author='Alex J. Hammel',
        author_email='ahammel87@gmail.com',
        packages=['srtools', 'srstats', 'sr_postgres'],
        url='http://pypi.python.org/pypi/srtools',
        license='LICENSE.txt',
        description='Tools for short-read genetic data.',
        long_description=open('README.txt').read(),
        install_requires=[
        "Python >= 3.2.3"
        ],
     )
