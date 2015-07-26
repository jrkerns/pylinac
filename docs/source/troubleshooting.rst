
.. _general_troubleshooting:

===============
Troubleshooting
===============

General
-------

Things always go wrong in real life. If you tried analyzing your data in pylinac and it threw an
error, you can try a few simple things to fix it.

* First, **See if the demo works** - If not, pylinac may not have installed correctly or you may not
  have a dependency or the minimum version of a dependency.
* Second, **Check the error** - If it's an error that makes sense, maybe you just forgot something; e.g.
  analyzing an image before it's loaded will raise an error. Asking for a gamma result of a fluence before
  calculating the fluence will raise an error. Such things are easy to fix.
* Third, **Check the Troubleshooting section of the specific module** - Each module may fail in different
  ways, and also have different methods of resolution.
* And if none of those work, **Email me or file an issue on Github** - jkerns100 at gmail.com and I'll do
  my best to figure out the problem; you may have found a bug and it needs fixing! You may also file an issue
  on the Github `Issue tracker <https://github.com/jrkerns/pylinac/issues>`_.

Loading TIFF Files
------------------

Loading TIFF files can be tricky since there are many variations of the TIFF image format.
`Pillow <https://python-pillow.github.io/>`_ is the package for image I/O in Python and is what
pylinac uses. But sometimes even Pillow has trouble. If you've tried loading a TIFF file and it
just doesn't seem to be working you should try resaving the image with another program. While I
can't tell you exactly what will work, one solution that's worked for me is using
`GIMP <http://www.gimp.org/>`_. It's free; just open up your TIFF files and then export them.
It may not seem like that should change anything, but my anecdotal evidence is that every TIFF
image that didn't work that I threw into GIMP allowed me to read it in.
