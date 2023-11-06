.. _contributor_guide:

============
Contributing
============

There are several ways you can contribute to the pylinac project no matter your skill level. Read on for more info.

Submitting bugs
---------------

The easiest way to contribute is to report bugs. Submit bugs via a Github issue `here <https://github.com/jrkerns/pylinac/issues>`__.

Submitting files
----------------

Another easy way to improve pylinac is to submit QA files for the testing repository. Files are treated anonymously and are
added to the test suite so that the package will become more robust. You can submit files `here <https://forms.gle/sfrDXL3XhHsyiKeJ7>`__.

Suggesting ideas
----------------

Ideas are always welcome (though they might not get implemented). You can submit new ideas `here <https://github.com/jrkerns/pylinac/issues>`_.


Commit changes
--------------

Now you're serious about contributing. Awesome! Pylinac has mostly had one maintainer but we are looking for newcomers to contribute.
There are a few things to know before contributing:

* Make an issue first, whether it be for a bug fix or feature request. This helps track the progress and put it on the roadmap appropriately.
* See if there are any other tools out there that can solve the problem. We don't want to write code just to write code. Pylinac should solve a problem
  that no one else has solved, do it better than existing solutions, or do so openly (vs closed). If another library already
  does this, the justification takes more work.
* Evaluate the work involved. This includes reviewing any existing modules or 3rd party tools that can help solve the problem.
* At this point, it looks like you're going to write some code. Make sure you have a
  :ref:`development environment <setup_dev_env>` set up.
* Propose the framework. This includes making the files and boilerplate for the new module/functions. This will allow others to evaluate and make
  suggestions to the framework before the work actually starts.
* Once the framework is agreed upon the code can start flowing. Get to it!
* If you're unsure, just ask.

.. _setup_dev_env:

Setting up a development environment
------------------------------------

If you want to contribute code, you'll need to set up a development environment. This is easy to do with the following steps.
If you're new to git, see the `git guide <https://git-scm.com/book/en/v2/Getting-Started-First-Time-Git-Setup>`__.
Be sure to ask for help if you need it!

* `Fork <https://github.com/jrkerns/pylinac/fork>`__ the pylinac repository.
* Clone your `forked repository to your local machine <https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository>`__: ``git clone`` followed by the URL of your forked repository.
  Some IDEs also have a GUI for this. See also :ref:`distro_stack`.
* Create a virtual environment. This is optional but highly recommended. See the `Python guide <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/>`__ to set up a new venv.
* Install the requirements and developer requirements: ``pip install -r requirements.txt -r requirements-dev.txt``
* Create a new branch for your work: ``git checkout -b my_new_branch``
* Make your changes
* Write tests for your changes. Most test modules in pylinac have a 1:1 correlation with the library modules and thus the
  expected location should be straightforward.
* Run the tests locally: ``pytest``
* If the tests pass, commit your changes: ``git commit -m "my new feature!"``
* Push your changes to your forked repository: ``git push origin my_new_branch``
* Make a new pull request to the main pylinac repository.
