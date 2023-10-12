Profiles
========

Profiles, in the context of pylinac, are 1D arrays of data that
contain a single radiation field. Colloquially, these are what
physicists might call an "inplane profile" or "crossplane profile"
although the usage here is not restricted to those orientations.

Use Cases
---------

Typical use cases for profiles are:

* Calculating a metric such as flatness, symmetry, penumbra, etc.
* Finding the center of the field.
* Finding the field edges.

Assumptions & Constraints
-------------------------

* There is one single, large radiation field in the profile.
* The radiation field should not be at the edge of the array.
* The radiation field should have a higher pixel value than the background. I.e. it should not be inverted.

* The field does not have to be normalized.
* The field can be "horned" (i.e. have a dip in the middle) and/or contain a peak from an FFF field.
* The field can be off-center, but the penumbra should be fully contained in the array.
* The field can be skewed (e.g. a wedge field).
* The field can have any resolution or no resolution.


Legacy vs New Classes
---------------------

The legacy class for analyzing profiles is :class:`~pylinac.core.profile.SingleProfile`.
This class is frozen and will not receive updates.

The modern classes for analyzing profiles are :class:`~pylinac.core.profile.FWXMProfile`,
:class:`~pylinac.core.profile.InflectionDerivativeProfile`, :class:`~pylinac.core.profile.HillProfile`.
These classes deal with **arrays only**.

For physical profiles, i.e. something where the values have a physical size or location like an EPID profile
or water tank scan, use the following classes: :class:`~pylinac.core.profile.FWXMProfilePhysical`,
:class:`~pylinac.core.profile.InflectionDerivativeProfilePhysical`, :class:`~pylinac.core.profile.HillProfilePhysical`.

The difference between ``SingleProfile`` and the other classes
is that ``SingleProfile`` is a swiss army knife. It can
do almost everything the other classes can do and considerably more. However,
the other classes are more specialized and thus more robust
as well as a lot clearer and focused.
Internally, pylinac uses the specialized classes.
A plugin system is being developed for adding additional functionality to the specialized classes
to mimic the ``SingleProfile`` functionality.

The ``SingleProfile`` class is more complicated to both read and use
than the specialized classes. It's also harder to test and maintain.
Thus, the specialized classes have come about as a response.



API
---

.. autoclass:: pylinac.core.profile.SingleProfile
    :inherited-members:
    :members:

.. autoclass:: pylinac.core.profile.FWXMProfile
    :inherited-members:
    :members:

.. autoclass:: pylinac.core.profile.FWXMProfilePhysical
    :inherited-members:
    :members:

.. autoclass:: pylinac.core.profile.InflectionDerivativeProfile
    :inherited-members:
    :members:

.. autoclass:: pylinac.core.profile.InflectionDerivativeProfilePhysical
    :inherited-members:
    :members:

.. autoclass:: pylinac.core.profile.HillProfile
    :inherited-members:
    :members:

.. autoclass:: pylinac.core.profile.HillProfilePhysical
    :inherited-members:
    :members:
