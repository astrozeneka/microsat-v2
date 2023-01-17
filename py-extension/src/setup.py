from distutils.core import Extension, setup
module1 = Extension("microsatellite", sources=["microsatellite.c"])
setup(name="microsatellite",
      version="1.0",
      description="Search for microsatellites in a genome",
      ext_modules=[module1])