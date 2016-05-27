from distutils.core import setup, Extension

module1 = Extension('RandomAngle', sources = ['Core.c'])

setup (name = 'Adv_Calc',
        version = '1.0',
        description = 'Random Or. a 3d vector in an angle alpha',
        ext_modules = [module1])