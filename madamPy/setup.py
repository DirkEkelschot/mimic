from distutils.core import setup, Extension

setup(name='voicePkg', version='1.0',  \
      ext_modules=[Extension('voice', ['voiceWrapper.cpp','Voice.cpp'])])
