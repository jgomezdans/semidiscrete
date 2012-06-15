import imp
import os
import os.path
import shutil
import sys
import tempfile
if os.name == 'nt': import win32api


def _tmp_pkg(dir="/tmp/"):
    """
    Create a temporary package.

    Returns (name, path)
    """
    while True:
        path = tempfile.mkdtemp(dir=dir)
        name = os.path.basename(path)
        try:
            modinfo = imp.find_module(name)
            # if name is found, delete and try again
            os.rmdir(path)
        except:
            break   
    init = file(os.path.join(path, '__init__.py'), 'w')
    init.close()
    return name, path


def mext(name):
    """
    Load and return a unique copy of a (possibly) already loaded module.

    This makes it possible to treat a module as a class and load multiple
    instances of it.
    """

    # first find the "real" module on the "real" syspath
    srcfile, srcpath, srcdesc = imp.find_module(name)
    # now create a temp directory for the bogus package
    pkgname, pkgdir = _tmp_pkg('.')
    # copy the original module to the new package
    shutil.copy(srcpath, pkgdir)
    # add the directory containing the new package to the search path
    #sys.path.append(tmpdir)
    # import the module
    # __import__ returns the package, not the sub-module
    pkg = __import__(pkgname, globals(), locals(), [name])
    # remove the bogus directory from sys.path ???
    #sys.path.remove(tmpdir)
    # return the module object
    print pkgdir
    return getattr(pkg, name)



class MExt(object):
    """
    Load a unique copy of a module that can be treated as a "class instance".
    """

    def __init__(self, name):
        self.name = name
        # first find the "real" module on the "real" syspath
        srcfile, srcpath, srcdesc = imp.find_module(name)
        # now create a temp directory for the bogus package
        self._pkgname, self._pkgdir = _tmp_pkg('.')
        # copy the original module to the new package
        shutil.copy(srcpath, self._pkgdir)
        # add the directory containing the new package to the search path
        #sys.path.append(tmpdir)
        # import the module
        # __import__ returns the package, not the sub-module
        self._pkg = __import__(self._pkgname, globals(), locals(), [self.name])
        # remove the bogus directory from sys.path ???
        #sys.path.remove(tmpdir)
        # return the module object
        self._module = getattr(self._pkg, self.name)
        # now add the module's stuff to this class
        self.__dict__.update(self._module.__dict__)

    def __del__(self):
        # remove module
        del sys.modules[self._module.__name__]
        del sys.modules[self._pkg.__name__]
        # on win32, the DLL must be unloaded forcefully in order to delete it.
        # on Darwin (other unix???) this doesn't appear to be necessary
        # try to unload the dll
        if os.name == 'nt':
            hModule = win32api.GetModuleHandle(self._module.__file__)
            win32api.FreeLibrary(hModule)
        # now try to delete the files and directory
        shutil.rmtree(self._pkgdir)
        # make sure the original module is loaded -
        # otherwise python crashes on exit
        # if MExt objects have not been explicitly 'del'd,
        # and __del__ is occurring on python shutdown, the import will fail
        # and the exception is caught here
        try:
            __import__(self.name)
        except:
            pass

class SemiDiscrete(object):
    def __init__(self, *args, **kwargs):
        semidiscrete_mod = MExt('rtmodel_ad_trans2')
        self.semidiscrete = semidiscrete_mod._module
