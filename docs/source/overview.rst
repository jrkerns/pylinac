
========================
Pylinac General Overview
========================

.. image:: images/starshot_analyzed.png

.. image:: images/vmat_analyzed.png

.. image:: images/PF_tight_tolerance.png

Intended Use
------------

Pylinac is intended to be used by two types of physicists: ones who know at least a bit of programming and those who know nothing about
programming.

For the first group, pylinac can be used within a Python environment to automate analysis of QA images.

For the second group, pylinac is also implemented as a `web app <www.assuranceqa.herokuapp.com>`_ and a desktop GUI.

Philosophy
----------

Pylinac runs on a few philosophical principles:

* A given module should only address 1 overarching task.
* Using pylinac should require a minimal amount of code.
* The user should have to supply as little information as necessary to run an analysis.
* The underlying code of pylinac should be easy to understand.

.. epigraph::
    The joy of coding Python should be in seeing short, concise, readable classes that express
    a lot of action in a small amount of clear code -- not in reams of trivial code that bores
    the reader to death.

    -- Guido van Rossum

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

Pylinac has a handful of modules, but most of them work somewhat the same, so here we describe the general patterns you'll see when using
pylinac.

* **Each module has its own demonstration method(s)** -- If you don't yet have an image or data and want to see how a module works
  you can run and inspect the code of the demo to get an idea. Most demo methods have a name like or starts with ``.run_demo()``.
* **Each module has its own demo image/dataset(s)** -- Want to test the analysis but are having trouble with your image? Use the provided
  demo images. All major classes have a demo image or dataset and are usually similar to ``.from_demo_image()``.
* **Each module has similar load, analyze, and show methods and behavior** -- The normal flow of a pylinac module use is to 1) Load the data in,
  2) Analyze the data, and 3) Show the results.
* **Most modules can be fully utilized in a few lines** -- The whole point of pylinac is to automate and simplify the process of
  analyzing routine QA data. Thus, most routines can be written in a few lines. Each module gives a starting script
  in its documentation.