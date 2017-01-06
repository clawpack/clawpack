# Installing the Python Clawpack tools

The following command is sufficient to get up and running with PyClaw, VisClaw, and the Python interfaces to other packages:

    pip install clawpack

If you want to install clawpack in a manner that allows also using the
Fortran codes (Classic, AMRClaw, and GeoClaw), see
http://www.clawpack.org/installing.html

# Installation Tools

## Git

If you are just interested in working with development repositories, those are available from a checked out repository with the following command:

    python setup.py git-dev

Scroll down to "Working with a Development Version of Clawpack" for more details.  

## Python

We use the [numpy distutils](http://docs.scipy.org/doc/numpy/reference/distutils.html) as our fundamental install tool because we are using extended f2py extension module installers.

We really like the pip installer.  If you don't have pip installed and will be working with the full Clawpack stack, we recommend that you use [Anaconda](http://www.continuum.io/downloads) to quickly get up and running with a user-space Scientific Python stack.  

# Working with a Development Version of Clawpack

## Clawpack Packages as Git Submodules

Although the Clawpack packages have been set up as independent repositories, occasionally changes in one module will need accompanying changes in other modules.  In this situation, it is a good idea for maintainers to make a "meta-commit" to this top-level repository that contains the other repositories as submodules.  This allows maintainers to coordinate changes across multiple Clawpack packages.  An example of how to do this is shown later.

## Creating a read-only development version of Clawpack

```
git clone git://github.com/clawpack/clawpack.git
cd clawpack
python setup.py git-dev
# optionally, for installation of Python components
pip install -e .
```

This downloads all of the clawpack modules as subrepositories checked out at specific commits (as opposed to the tip of a branch).  

## Contributing to Clawpack!

If you'd like to contribute to Clawpack, you do not need to worry about the submodules.  Just checkout a branch in the package you are updating and push your changes as a pull request.  Just be sure to post your pull requests at the same time (and refer to the other pull requests) when submitting coordinated changes.

### Steps for contributors
* make local changes in one or more packages
* test your changes, if they're not testable, please add a test!
* commit and publish (push) your changes
* submit pull requests for each package you've modified

Optionally, you may also create a top-level clawpack branch and submit a pull request showing how your changes work together.  This is especially recommended if you are unable to create working tests otherwise.  Remember that you need to push commits in your submodules to GitHub before you can open a pull request that refers to those commits.

## Creating a development version of Clawpack for modification on GitHub/submitting pull requests

First, head to http://github.com/clawpack and fork the repositories that you will be working with.  You will probably need to fork one or more of the following:

* amrclaw
* clawpack
* clawutil
* geoclaw
* pyclaw
* riemann
* visclaw

Then, for each submodule you'd like to modify, add your own repository as a remote.  Here's an example for publishing changes to geoclaw (replace `username` below with your own username).

```
cd geoclaw
git remote add username git@github.com:username
# you may want to start from upstream
git checkout -b new_feature origin/master
# make some changes
git commit
# push your changes to your GitHub repository
git push username 
# open pull request on GitHub
```

## Working with a Maintenance Version of Clawpack

(This section is a work-in-progress)

# Maintaining Clawpack

The Clawpack maintainers are responsible for reviewing and accepting pull requests, although community support is always welcome in reviewing incoming pull requests.  Once a pull request has been accepted and tested, the maintainers need to decide if they'd like to update the top-level clawpack repository, which should be treated as the stable development branch coordinating the different repositories.  Changes in the individual repositories are not visible to somebody using `pip install` or the `python setup.py git-dev` commands unless they manually fetch/checkout changes or the top-level clawpack master branch has been updated.

## Complete Submodule Workflow Example

```
cd pyclaw
# make local changes
git commit -a -m "Updated pyclaw *cross-changes* with riemann."
# push to your GitHub fork of pyclaw
git push
cd ../riemann
# make local changes
git commit -a -m "Updated riemann *cross-changes* with pyclaw."
# push to your GitHub fork of riemann
git push
cd ..
# now at top-level clawpack directory
git add pyclaw riemann
git commit -m "Cross-update pyclaw and riemann."
# push to your GitHub fork of clawpack
git push

# open PRs on GitHub for all branches that need to be merged

```

More reading about submodules.
* [Submodules in the Git Community Book](http://book.git-scm.com/5_submodules.html)
* [Submodules in the Pro Git book](http://progit.org/book/ch6-6.html)
