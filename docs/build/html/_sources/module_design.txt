
=====================
General Module Design
=====================

Pylinac has a handful of modules, but they all work somewhat the same, so here we describe the general patterns you'll see when using
pylinac.

* **Each module has its own demonstration method** -- If you don't yet have an image or data and want to see how a module works
  you can run and inspect the code of the demo to get an idea. Most demo methods have a name similar to ``.run_demo*()``.
* **Each module has similar load, analyze, and show methods** -- The normal flow of a pylinac module script is to 1) Load the data in,
  2) Analyze the data, and 3) Show the results
* **Most modules can be fully utilized in a few lines** -- The whole point of pylinac is to automate and simplify the process of
  analyzing routine QA data. Thus, most routines can be written in a few lines.