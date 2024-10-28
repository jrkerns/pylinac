=================
Memory Efficiency
=================

Generally speaking, Python is not the most memory-efficient language.
Most of the time this doesn't make a significant difference, but in
certain context with constraints it can be a problem.
It's unlikely that you'll actually need to worry about this, but this
provides some context. Within the RadMachine environment, memory can be
a precious commodity, so some work has been done to optimize memory usage.
It is explained here for clarity and completeness.

.. note::

    Memory optimization is only relevant here for DICOM stacks (CT, MR, etc).

.. note::

    "Space" here refers to both memory **and** disk space. In the RadMachine cloud environment,
    disk space is the same as memory so for the following discussion there is no benefit to using memory vs disk.

There are 3 uses of space when analyzing a CT dataset:

* The original dataset as passed (e.g. a zip archive or a list of files)
* The extracted dataset (if applicable)
* The array(s) in memory

In the worst case scenario all three uses are employed.
E.g. a 5MB zip archive is passed.
It is then extracted as a 30MB directory,
and then each image is loaded into numpy as an array (usually uint16 or float).
Depending on the CT scanner this can be ~2.5MB/image.
So in total the we might be using 5 + 30 + 60 * 2.5 = 185MB of space,
assuming 60 images.
This doesn’t sound like that much space but this can be much larger depending on the
scanner and number of slices.

We can’t ever get rid of the first dataset, but we can in theory get rid of the latter two.

The loaded arrays in memory are the easiest to get rid of. To do so, we can simply
load the given image requested on demand, or, "lazily". This will result in only 1
image in memory at any given time. This is referred to as "lazy".

The extracted dataset can be more difficult. Assuming the original dataset is also
counted toward the space used, we want to minimize what we can. This means
we ideally want to be passed a compressed dataset (ZIP archive) and also load one
image lazily as requested. This can be done in vanilla Python via the ``zipfile`` module but suffers from
two problems:

- It can modify the data in-place. I.e. performing a modification actually modifies
  the original data. This is undesirable.
- Changing or removing data from the archive is difficult and slow. Python can easily
  read from a ZIP archive but writing to it or removing data requires rewriting the
  entire archive. This is extremely slow.

The best solution so far is to create a "shadow" object that holds the data in memory,
but is also compressed. This is a compromise between the two. Generally speaking,
the total space usage for this method is around 2x the size of the compressed dataset.
This is referred to as "shadow object".

See the table below for a comparison:

.. note::

    This assumes a ZIP archive is passed to the constructor. Numbers will vary depending on
    the dataset and the computer.

+-----------------------------------+-------------------------+----------------+---------------+
|                                   | non-optimized (default) | lazy           | shadow object |
+===================================+=========================+================+===============+
| Slowdown (%)                      |                         | 5-15%          | 10-20%        |
+-----------------------------------+-------------------------+----------------+---------------+
| zip archive (MB)                  | 13                      | 13             | 13            |
+-----------------------------------+-------------------------+----------------+---------------+
| extracted zip size (MB)           | 49                      | 49             | 0             |
+-----------------------------------+-------------------------+----------------+---------------+
| Shadow object (MB)                | 0                       | 0              | 13            |
+-----------------------------------+-------------------------+----------------+---------------+
| DICOM stack size at any time (MB) | 195                     | 2.5            | 2.5           |
+-----------------------------------+-------------------------+----------------+---------------+
| Max usage (MB)                    | 257                     | 64.5           | 28.5          |
+-----------------------------------+-------------------------+----------------+---------------+
| Memory usage vs dataset (x)       | 19.8                    | 5.0            | 2.2           |
+-----------------------------------+-------------------------+----------------+---------------+
| Memory optimizations (x)          |                         | 4.0            | 9.0           |
+-----------------------------------+-------------------------+----------------+---------------+

Which one should you use? It depends on your environment, but if memory is no concern,
the default is fine. If you are memory-constrained but not disk-constrained, extracted lazy
is the best. If you are both memory and disk constrained, shadow object is the best.

How to enable memory optimizations
----------------------------------

Memory efficiency is not enabled by default for historical reasons. You generally
don't have to worry about implementation details. If you want to enable it,
you can do so in the constructor of any of the DICOM stack classes (CatPhan, ACR MR/CR, Cheese, etc)
by passing ``memory_efficient_mode=True``.

* When ``memory_efficient_mode=False`` (default), the non-optimized method is used and everything is loaded into memory.
* When ``memory_efficient_mode=True`` and the dataset is not a ZIP archive (a list of files already on disk), the lazy method is used.
* When ``memory_efficient_mode=True`` and the dataset is a ZIP archive, the shadow object method is used.
