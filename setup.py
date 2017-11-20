from setuptools import setup
def readme():
    with open('README.md') as f:
        return f.read()
setup(name='chem_G5',
      version='1.0',
      description='Elementray reaction and reversible reaction',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Text Processing :: Linguistic',
      ],
      url='https://github.com/CS207-G5/cs207-FinalProject',
      author='Brianna, Thomas, Yujiao, Yi',
      author_email='yi_zhai@g.harvard.edu',
      license='MIT',
      packages=['chem_G5'],
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose'])