.. _rigid_transformations:

Rigid transformations
=====================

In mathematics, a rigid transformation (also called Euclidean transformation or Euclidean isometry)
is a geometric transformation of a Euclidean space that preserves the Euclidean distance between every pair of points.
In the context of pylinac rigid transformations are used to place ROIs and include rotations and translations while reflections are intentionally left out.


In 2D
-----

Translations
~~~~~~~~~~~~

In Euclidean geometry, a translation is a geometric transformation that moves every point of a figure,
shape or space by the same distance in a given direction.

.. math::

    p^{*} = p + t = \begin{bmatrix}p_x\\p_y\end{bmatrix} + \begin{bmatrix}t_x\\t_y\end{bmatrix}

.. plot::
    :include-source: False

    import matplotlib.pyplot as plt
    import numpy as np

    p0 = np.array([[2, 4],[2, 3],[3,3]]).T
    t = np.array([[1, -2]]).T
    p1 = p0 + t

    ax1 = plt.subplot(1, 1, 1)
    plt.plot(p0[0], p0[1], 'ko-')
    plt.plot(p1[0], p1[1], 'ro-')
    ax1.annotate('', xy=(float(p1[0,1]), float(p1[1,1])), xytext=(float(p0[0,1]), float(p0[1,1])),
                 arrowprops=dict(facecolor='black', shrink=0.1, width=0.1, headlength=5, headwidth=5), )
    ax1.annotate('T', xy=(float(0.5 * p0[0,1] + 0.5 * p1[0,1]) + 0.1, float(0.5 * p0[1,1] + 0.5 * p1[1,1]) + 0.1))
    ax1.annotate('p', xy=(float(p0[0,1] + 0.15), float(p0[1,1] + 0.15)))
    ax1.annotate('p*', xy=(float(p1[0,1] + 0.15), float(p1[1,1] + 0.15)))
    plt.axis((0, 5, 0, 5))
    ax1.set_aspect('equal')
    plt.show()


Rotations
~~~~~~~~~

In Euclidean geometry, a rotation is a geometric transformation that turns every point of a figure,
shape or space around a fixed point through a specified angle and in a specified direction

.. math::

    p^* = R * p = \begin{bmatrix}cos(\alpha)&&-sin(\alpha)\\sin(\alpha)&&cos(\alpha)\end{bmatrix} * \begin{bmatrix}p_x\\p_y\end{bmatrix}

.. plot::
    :include-source: False

    import matplotlib.pyplot as plt
    import numpy as np
    from skimage import transform

    p0 = np.array([[3, 2],[3, 1],[4,1]])
    a = 30
    tf = transform.EuclideanTransform(rotation=np.deg2rad(a))
    p1 = transform.matrix_transform(p0, tf.params)
    t = np.linspace(np.arctan2(p0[1,1],p0[1,0]), np.arctan2(p1[1,1],p1[1,0]), 100)
    r = np.linalg.norm(p0[1,:])

    ax1 = plt.subplot(1, 1, 1)
    plt.plot(p0[:,0], p0[:,1], 'ko-')
    plt.plot(p1[:,0], p1[:,1], 'ro-')
    plt.plot(r*np.cos(t), r*np.sin(t), 'k--')
    plt.plot([0, p0[1,0]], [0, p0[1,1]], 'k:')
    plt.plot([0, p1[1,0]], [0, p1[1,1]], 'k:')

    ax1.annotate(r'$\alpha$', xy=(float(0.5 * p0[0,1] + 0.5 * p1[0,1]) - 0.15, float(0.5 * p0[1,1] + 0.5 * p1[1,1]) - 0.15))
    ax1.annotate('p', xy=(float(p0[1,0] + 0.15), float(p0[1,1] + 0.15)))
    ax1.annotate('p*', xy=(float(p1[1,0]), float(p1[1,1] + 0.15)))
    plt.axis((0, 5, 0, 5))
    ax1.set_aspect('equal')

    plt.show()


Pose
~~~~

In 2D geometry, a pose (or transform) is a combination of a translation and a rotation:

.. math::
    p^* = R * p + T

The way to define this transformation is by using 3x3 matrices:

.. math::
    Translation = \left[\begin{array}{cc:c}1&0&T_x\\0&1&T_y\\\hdashline0&0&1\end{array}\right], Rotation = \left[\begin{array}{cc:c}cos(\alpha)&-sin(\alpha)&0\\sin(\alpha)&cos(\alpha)&0\\\hdashline0&0&1\end{array}\right]

.. math::
    Transformation = Translation * Rotation

.. math::
    p^* = Transformation * p => \begin{bmatrix}p^*_x\\p^*_y\\1\end{bmatrix} = \left[\begin{array}{cc:c}cos(\alpha)&-sin(\alpha)&T_x\\sin(\alpha)&cos(\alpha)&T_y\\\hdashline0&0&1\end{array}\right] * \begin{bmatrix}p_x\\p_y\\1\end{bmatrix}


Coordinate system and order of operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since matrix multiplications are not commutative, the order by which each transformation is applied matters.
Furthermore the concept of transformation also depends on the frame of reference of the user.
We can define two frames of reference:

**Extrinsic (space‑fixed) coordinates**: the axes stay put in the “world” frame, and each transformation is performed about one of those fixed axes.

**Intrinsic (body‑fixed) coordinates**: the axes ride along with the object, i.e. the next transformation is with respect to the new axes.

Let's look at some examples:

* **First rotation then translation (extrinsic coordinates)**:

.. math::
    Transformation = Translation * Rotation

.. note::
   Using ``scikit-image`` library, the equivalent is:

   .. math::
      Transformation = EuclideanTransform(rotation=R) + EuclideanTransform(translation=T)

   .. math::
      Transformation = EuclideanTransform(rotation=R, translation=T)

.. plot::
    :include-source: False

    import matplotlib.pyplot as plt
    import numpy as np
    from skimage import transform

    p0 = np.array([[0, 1],[0, 0],[1,0]])
    a = 30
    tf1 = transform.EuclideanTransform(rotation=np.deg2rad(a))
    p1 = transform.matrix_transform(p0, tf1.params)
    t = np.linspace(np.arctan2(p0[2,1],p0[2,0]), np.arctan2(p1[2,1],p1[2,0]), 100)
    r = np.linalg.norm(p0[2,:])
    ax1 = plt.subplot(1, 1, 1)
    plt.plot(p0[:,0], p0[:,1], 'ko-')
    plt.plot(p1[:,0], p1[:,1], 'ro-.')
    plt.plot(r*np.cos(t), r*np.sin(t), 'k--')
    plt.plot([0, p0[1,0]], [0, p0[1,1]], 'k:')
    plt.plot([0, p1[1,0]], [0, p1[1,1]], 'k:')
    ax1.annotate(r'$\alpha$', xy=(float(0.5 * p0[2,0] + 0.5 * p1[2,0]) - 0.15, float(0.5 * p0[2,1] + 0.5 * p1[2,1]) - 0.15))
    ax1.annotate('p', xy=(float(p0[1,0] - 0.0), float(p0[1,1] - 0.2)))

    tf2 = transform.EuclideanTransform(translation=(2,0))
    p2 = transform.matrix_transform(p1, tf2.params)
    plt.plot(p2[:,0], p2[:,1], 'ro-')
    ax1.annotate('', xy=(float(p2[1, 0]), float(p2[1, 1]-0.2)), xytext=(float(p1[1, 0]), float(p1[1, 1]-0.2)),
                 arrowprops=dict(facecolor='black', shrink=0.1, width=0.1, headlength=5, headwidth=5), )
    ax1.annotate('T', xy=(float(0.5 * p2[1, 0] + 0.5 * p1[1, 0] - 0.0), float(0.5 * p2[1, 1] + 0.5 * p1[1, 1] - 0.4)))
    ax1.annotate('p*', xy=(float(p2[1, 0] - 0.0), float(p2[1, 1] - 0.2)))

    plt.axis((-1, 4, -1, 2))
    ax1.set_aspect('equal')

    plt.show()


* **First translation then rotation (extrinsic coordinates)**:

.. math::
    Transformation = Rotation * Translation

.. note::
   Using ``scikit-image`` library, the equivalent is:

   .. math::
       Transformation = EuclideanTransform(translation=T) + EuclideanTransform(rotation=R)

.. plot::
    :include-source: False

    import matplotlib.pyplot as plt
    import numpy as np
    from skimage import transform
    p0 = np.array([[0, 1],[0, 0],[1,0]])
    a = 30
    tf1 = transform.EuclideanTransform(translation=(2,0))
    p1 = transform.matrix_transform(p0, tf1.params)
    ax1 = plt.subplot(1, 1, 1)
    plt.plot(p0[:,0], p0[:,1], 'ko-')
    plt.plot(p1[:,0], p1[:,1], 'ro-.')
    ax1.annotate('', xy=(float(p1[1, 0]), float(p1[1, 1]-0.2)), xytext=(float(p0[1, 0]), float(p0[1, 1]-0.2)),
                 arrowprops=dict(facecolor='black', shrink=0.1, width=0.1, headlength=5, headwidth=5), )
    ax1.annotate('T', xy=(float(0.5 * p0[1, 0] + 0.5 * p1[1, 0]) + 0.0, float(0.5 * p0[1, 1] + 0.5 * p1[1, 1] - 0.4)))
    ax1.annotate('p', xy=(float(p0[1, 0] - 0.15), float(p0[1, 1] - 0.15)))

    tf2 = transform.EuclideanTransform(rotation=np.deg2rad(a))
    p2 = transform.matrix_transform(p1, tf2.params)
    t = np.linspace(np.arctan2(p1[1,1],p1[1,0]), np.arctan2(p2[1,1],p2[1,0]), 100)
    r = np.linalg.norm(p1[1,:])
    plt.plot(p2[:,0], p2[:,1], 'ro-')
    plt.plot(r*np.cos(t), r*np.sin(t), 'k--')
    plt.plot([0, p1[1,0]], [0, p1[1,1]], 'k:')
    plt.plot([0, p2[1,0]], [0, p2[1,1]], 'k:')
    ax1.annotate(r'$\alpha$', xy=(float(0.5 * p2[1,0] + 0.5 * p1[1,0]) - 0.15, float(0.5 * p2[1,1] + 0.5 * p1[1,1]) - 0.15))
    ax1.annotate('p*', xy=(float(p2[1, 0] + 0.0), float(p2[1, 1] + 0.15)))

    plt.axis((-1, 4, -1, 3))
    ax1.set_aspect('equal')

    plt.show()

* **First translation then rotation (intrinsic coordinates)**:

.. math::
    Transformation = Rotation' * Translation

where ``Rotation'`` represents the intrinsic frame of reference

.. note::
   .. math::
      Transformation = R_{intrinsic} * T_{intrinsic} = T_{extrinsic} * R_{extrinsic}

.. plot::
    :include-source: False

    import matplotlib.pyplot as plt
    import numpy as np
    from skimage import transform
    p0 = np.array([[0, 1],[0, 0],[1,0]])
    a = 30
    tf1 = transform.EuclideanTransform(translation=(2,0))
    p1 = transform.matrix_transform(p0, tf1.params)
    ax1 = plt.subplot(1, 1, 1)
    plt.plot(p0[:,0], p0[:,1], 'ko-')
    plt.plot(p1[:,0], p1[:,1], 'ro-.')
    ax1.annotate('', xy=(float(p1[1, 0]), float(p1[1, 1]-0.2)), xytext=(float(p0[1, 0]), float(p0[1, 1]-0.2)),
                 arrowprops=dict(facecolor='black', shrink=0.1, width=0.1, headlength=5, headwidth=5), )
    ax1.annotate('T', xy=(float(0.5 * p0[1, 0] + 0.5 * p1[1, 0]) + 0.0, float(0.5 * p0[1, 1] + 0.5 * p1[1, 1] - 0.4)))
    ax1.annotate('p', xy=(float(p0[1, 0] + 0.0), float(p0[1, 1] - 0.2)))

    tf = transform.EuclideanTransform(rotation=np.deg2rad(a)) + tf1
    p2 = transform.matrix_transform(p0, tf.params)
    t = np.linspace(np.arctan2(p0[2,1],p0[2,0]), np.arctan2(p2[0,1],p2[0,0]), 100)
    r = np.linalg.norm(p0[0,:])
    plt.plot(p2[:,0], p2[:,1], 'ro-')
    plt.plot(p1[1,0]+r*np.cos(t), p1[1,1]+r*np.sin(t), 'k--')
    #plt.plot([0, p1[1,0]], [0, p1[1,1]], 'k:')
    #plt.plot([0, p2[1,0]], [0, p2[1,1]], 'k:')
    ax1.annotate(r'$\alpha$', xy=(float(0.5 * p2[2,0] + 0.5 * p1[2,0]) - 0.15, float(0.5 * p2[2,1] + 0.5 * p1[2,1]) - 0.15))
    ax1.annotate('p*', xy=(float(p2[1, 0] + 0.0), float(p2[1, 1] - 0.2)))

    plt.axis((-1, 4, -1, 2))
    ax1.set_aspect('equal')

    plt.show()


* **First rotation then translation (intrinsic coordinates)**:

.. math::
    Transformation = Translation' * Rotation = Rotation * Translation

where ``Translation'`` represents the intrinsic frame of reference

.. note::
   .. math::
      Transformation = T_{intrinsic} * R_{intrinsic} = R_{extrinsic} * T_{extrinsic}

.. plot::
    :include-source: False

    import matplotlib.pyplot as plt
    import numpy as np
    from skimage import transform
    p0 = np.array([[0, 1],[0, 0],[1,0]])
    a = 30
    tf1 = transform.EuclideanTransform(rotation=np.deg2rad(a))
    p1 = transform.matrix_transform(p0, tf1.params)
    t = np.linspace(np.arctan2(p0[2,1],p0[2,0]), np.arctan2(p1[2,1],p1[2,0]), 100)
    r = np.linalg.norm(p0[2,:])
    ax1 = plt.subplot(1, 1, 1)
    plt.plot(p0[:,0], p0[:,1], 'ko-')
    plt.plot(p1[:,0], p1[:,1], 'ro-.')
    plt.plot(r*np.cos(t), r*np.sin(t), 'k--')
    plt.plot([0, p0[1,0]], [0, p0[1,1]], 'k:')
    plt.plot([0, p1[1,0]], [0, p1[1,1]], 'k:')
    ax1.annotate(r'$\alpha$', xy=(float(0.5 * p0[2,0] + 0.5 * p1[2,0]) - 0.15, float(0.5 * p0[2,1] + 0.5 * p1[2,1]) - 0.15))
    ax1.annotate('p', xy=(float(p0[1,0] - 0.0), float(p0[1,1] - 0.2)))

    tf = transform.EuclideanTransform(translation=(2,0)) + tf1
    p2 = transform.matrix_transform(p0, tf.params)
    plt.plot(p2[:,0], p2[:,1], 'ro-')
    ax1.annotate('', xy=(float(p2[1, 0]), float(p2[1, 1]+0.1)), xytext=(float(p1[1, 0]), float(p1[1, 1]+0.1)),
                 arrowprops=dict(facecolor='black', shrink=0.1, width=0.1, headlength=5, headwidth=5), )
    ax1.annotate('T', xy=(float(0.5 * p2[1, 0] + 0.5 * p1[1, 0] + 0.0), float(0.5 * p2[1, 1] + 0.5 * p1[1, 1] + 0.2)))
    ax1.annotate('p*', xy=(float(p2[1, 0] - 0.0), float(p2[1, 1] - 0.2)))

    plt.axis((-1, 4, -1, 3))
    ax1.set_aspect('equal')

    plt.show()
