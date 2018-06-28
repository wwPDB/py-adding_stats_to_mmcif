from setuptools import setup, find_packages



setup(
    name='adding_stats_to_mmcif',
    version='0.1',
    url='https://github.com/berrisfordjohn/adding_stats_to_mmcif',
    author='John Berrisford',
    test_suite = 'tests',
    #test_suite='tests.my_test_suite',
    #test_suite='nose.collector',
    #tests_require=['nose'],
    zip_safe=False,
    packages=find_packages('src'),
    #dependency_links=['https://github.com/project-gemmi/gemmi.git']
    #install_requires=['unittest', 'logging', 'xml', 'os'], 
)