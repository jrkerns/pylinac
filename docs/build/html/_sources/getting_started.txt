

===============
Getting Started
===============

Getting started with pylinac is easy! Once installed, you can write your own script in a matter of minutes. Each module of pylinac
addresses the topic of its name (e.g. ''VMAT'' performs VMAT-related QA tests), and all modules are importable from pylinac. Furthermore,
each module is designed as similarly as possible to one another. So once you start using one module,
it's easy to use another. Each module also has its own demo to show off what it can do.

Let's get started by running a demo with the ''VMAT'' module. First, import pylinac and its VMAT module::

    from pylinac import VMAT

This imports the VMAT class, which has all the capabilities of loading and analyzing VMAT EPID images. Let's 1) create an instance of that
class and then 2) run its demonstration method::

    myvmat = VMAT()
    myvmat.run_demo_drgs()

Running this should result in a printing of information to the console and an image showing the Open field and MLC field, like so:

#TODO: embed print results and matplotlib image.

Congratulations! In 3 lines you've successfully used a pylinac module. Of course there's more to it than that; you'll want to analyze your
own images.