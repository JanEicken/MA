from distutils.core import setup, Extension

module1 = Extension('trajcomp',
                    define_macros = [('MAJOR_VERSION', '0'),
                                     ('MINOR_VERSION', '1')],
#                    include_dirs = ['/usr/local/include','/usr/include/python3.2mu'],
#                    libraries = ['boost'],
#                    library_dirs = ['/usr/include/boost'],
                    include_dirs = ['/usr/include'],
					extra_compile_args = ['-std=c++11', '-lboost_filesystem', '-lboost_system'],
                    sources = ['trajcomp_module.cpp'])

setup (name = 'trajcomp',
       version = '0.1',
       description = 'This is a testing package',
       author = 'Martin Werner',
       author_email = 'martin@trajectorycomputing.com',
       url = 'http://trajectorycomputing.com',
       long_description = '''
This is just to see, whether my trajcomp can be used in python efficiently
''',
       ext_modules = [module1])
