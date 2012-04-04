# packages as git submodules

Although the clawpack repositories have been set up as independent modules, occasionally changes in one module will need accompanying changes in other modules.  In this situation, it is a good idea to make a "meta-commit" to a super-repository that contains the other repositories as submodules.

This allows developers to work in one repository while being assured that all other repositories will interoperate correctly.  

## Procedure for initializing the Clawpack superrepository:

```
# change this before the pull request is accepted to clawpack
git clone git://github.com/ahmadia/clawpack.git
cd clawpack
git submodule init
git submodule update
```

This brings down all of the clawpack modules as subrepositories checked out at specific commits (as opposed to the tip of a branch).  The head for each subrepository is detached, so you will need to check out master, then create a branch from there to make changes:

```
cd pyclaw
git checkout master
git checkout -b new_change
```

[Submodules in the Git Community Book](http://book.git-scm.com/5_submodules.html)
[Submodules in the Pro Git book](http://progit.org/book/ch6-6.html)

## Procedure for a Clawpack cross-repository commit:

* make local changes in each module
* commit and publish (push) modular changes
* change to top-level clawpack directory
* commit and push 

```
cd pyclaw
# make local changes
git commit -a -m "Updated pyclaw *cross-changes* with riemann."
git push
cd ../riemann
# make local changes
git commit -a -m "Updated riemann *cross-changes* with pyclaw."
git push
cd ..
git add pyclaw riemann
git commit -m "Cross-update pyclaw and riemann."
git push
```

# distribute installer 

The Python Distribution Utilities ([distutils](http://docs.python.org/distutils/index.html)) are the official Python standard for package installation.  A well-known enhancement to distutils is [setuptools](http://pypi.python.org/pypi/setuptools), which has been unmaintained since 2009.  A fork of setuptools called [distribute](http://pypi.python.org/pypi/distribute) attempts to replace setuptools, and seems to be commonly used in several larger projects.  A new project named [distutils2](http://pypi.python.org/pypi/Distutils2) will eventually supercede distutils, and is on the path for incorporation into Python 3.3 (and possibly backported into Python 2.*).

We are using the [numpy distutils](http://docs.scipy.org/doc/numpy/reference/distutils.html) in clawpack for now because we are using extended f2py extension module installers.

Some bootstrap instructions are available from [The Hitchhiker's Guide to Packaging](http://guide.python-distribute.org/).

[Some useful comments from Tarek Ziade](http://ziade.org/2010/03/03/the-fate-of-distutils-pycon-summit-packaging-sprint-detailed-report/)

# test refactor
