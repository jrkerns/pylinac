
========================
Pylinac General Overview
========================

.. image:: /images/analyzed_starshot.png
   :height: 340px
   :width: 300px

.. image:: /images/vmat_analyzed.png
   :height: 340px
   :width: 800px


Intended Use
------------

Pylinac is intended to be used by two types of physicists: ones who know at least a bit of programming and those who know nothing about
programming.

For the first group, pylinac can be used within a Python environment to automate analysis of tests, either as a script or used
within their own flavor of graphical user interface (GUI).

For the second group, pylinac will come with its own GUI, allowing those who don't know, or care to know, programming. The pylinac GUI will
function like any other "normal" end-user software, with a display and buttons, etc. The pylinac GUI component of a module will be built on
the underlying pylinac substructure, and will be built after the API component.

Philosophy
----------

Pylinac runs on a few philosophical principles:

* A given module should only address 1 overarching task
* Using pylinac should require a minimal amount of code
* The user should have to supply as little information as necessary to run an analysis
* The underlying code of pylinac should be easy to understand from the first look

Algorithm Design Overview
-------------------------

Generally speaking, the design of algorithms should all follow the same guidelines and appear as similar as possible. Each module will
outline its own specific algorithm in its documentation.

* Descriptions of algorithms are sorted into steps of the following:

  * **Allowances** -- These describe what the pylinac algorithm *can* account for.
  * **Restrictions** -- These are the things pylinac *cannot* do and must be addressed before the module can be properly used.
  * **Pre-Analysis** -- Algorithm steps that prepare for the main algorithm sequence.
  * **Analysis** -- The steps pylinac takes to analyze the image or data.
  * **Post-Analysis** -- What pylinac does or can do after analysis, like showing the data or checking against tolerances.

* Algorithm steps should be expressible in a word or short phrase.
* Algorithm method names should be as similar as possible from module to module.

.. _module_design:

Module Design
-------------

Pylinac has a handful of modules, but they all work somewhat the same, so here we describe the general patterns you'll see when using
pylinac.

* **Each module has its own demonstration method** -- If you don't yet have an image or data and want to see how a module works
  you can run and inspect the code of the demo to get an idea. Most demo methods have a name similar to ``.run_demo*()``.
* **Each module has similar load, analyze, and show methods** -- The normal flow of a pylinac module script is to 1) Load the data in,
  2) Analyze the data, and 3) Show the results
* **Most modules can be fully utilized in a few lines** -- The whole point of pylinac is to automate and simplify the process of
  analyzing routine QA data. Thus, most routines can be written in a few lines.


GUI Design
----------

The design of the graphical user interface (GUI) follows a general philosophy: *There is often too much information to show all at once,
so the GUI should only show the information relevant for a given step.* This is accomplished by sorting information into "layers". Layers
consist of showing only information relevant at that moment, and hiding the rest, until needed or asked for.

For example, in a given situation, there are usually several steps to get to the results:

* Information must be imported into the system. There is no need to take up space with an image or image holder or quantitative analysis at
  this time. Analysis settings are not yet relevant.
* Settings should be reviewed or changed. There is still no need to view the information, nor space reserved for loading images.
* Viewing information qualitatively. When viewing information qualitatively (like an image or plot of results),
  settings are only relevant as a backdrop to the situation.
* Viewing information quantitatively. Quantitative information is often viewed separate (at least in attention focus) from qualitative
  information.

In general, each of these situations requires little simultaneous information from the other. Visually speaking,
these "layers" can be hidden when not needed to save screen space.