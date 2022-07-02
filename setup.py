import os
import platform
import re
import subprocess
import sys

from setuptools import setup, Extension, find_namespace_packages
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir='', namespace=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)
        self.namespace = namespace


class CMakeBuild(build_ext):
    def run(self):
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.namespace+ext.name)))
        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      '-DBINDINGS=ON'
                      ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            cmake_args += ['-T', 'ClangCl']
            if sys.maxsize > 2 ** 32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            cpus = os.cpu_count()
            if cpus is None:
                cpus = 2
            build_args += ['--', '-j{}'.format(cpus)]

        # cross-compile support for macOS - respect ARCHFLAGS if set
        if sys.platform.startswith("darwin"):
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                arch_argument = "-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))
                print('macOS building with: ', arch_argument, flush=True)
                cmake_args += [arch_argument]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.', '--target', ext.name] + build_args, cwd=self.build_temp)


README_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                           'README.md')
with open(README_PATH) as readme_file:
    README = readme_file.read()

setup(
    name='mqt.qecc',
    version='0.1.0',
    author='Lucas Berent',
    author_email='lucas.berent@tum.de',
    description='MQT QECC - QECC - A Tool for Quantum Error Correcting Codes Written in C++',
    long_description=README,
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/lucasberent/qecc",
    ext_modules=[CMakeExtension('pyqecc', namespace='mqt.qecc.')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    packages=find_namespace_packages(include=['mqt.*']),
    classifiers=[
        'Development Status :: 4 - Beta',
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Electronic Design Automation (EDA)",
    ],
    keywords="mqt quantum error-correcting codes",
    project_urls={
        'Source': 'https://github.com/lucasberent/qecc',
        'Tracker': 'https://github.com/lucasberent/qecc/issues'
    }
)
