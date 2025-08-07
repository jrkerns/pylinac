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
    p1 = tf(p0)
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

Since matrix multiplications are not commutative, the order by which each individual transformation is applied matters.
Furthermore the concept of transformation also depends on the frame of reference of the user.
Two frames of reference can be defined:

**Extrinsic (space‑fixed) coordinates**: the axes stay put in the “world” frame, and each transformation is performed about one of those fixed axes.

**Intrinsic (body‑fixed) coordinates**: the axes ride along with the object, i.e. the next transformation is with respect to the new axes.

Let's look at some examples:

* **First rotation then translation (extrinsic coordinates)**:

.. math::
    Transformation = Translation * Rotation

.. plot::
    :include-source: False

    import matplotlib.pyplot as plt
    import numpy as np
    from skimage import transform

    p0 = np.array([[0, 1],[0, 0],[1,0]])
    a = 30
    tf1 = transform.EuclideanTransform(rotation=np.deg2rad(a))
    p1 = tf1(p0)
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
    p2 = tf2(p1)
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

.. plot::
    :include-source: False

    import matplotlib.pyplot as plt
    import numpy as np
    from skimage import transform
    p0 = np.array([[0, 1],[0, 0],[1,0]])
    a = 30
    tf1 = transform.EuclideanTransform(translation=(2,0))
    p1 = tf1(p0)
    ax1 = plt.subplot(1, 1, 1)
    plt.plot(p0[:,0], p0[:,1], 'ko-')
    plt.plot(p1[:,0], p1[:,1], 'ro-.')
    ax1.annotate('', xy=(float(p1[1, 0]), float(p1[1, 1]-0.2)), xytext=(float(p0[1, 0]), float(p0[1, 1]-0.2)),
                 arrowprops=dict(facecolor='black', shrink=0.1, width=0.1, headlength=5, headwidth=5), )
    ax1.annotate('T', xy=(float(0.5 * p0[1, 0] + 0.5 * p1[1, 0]) + 0.0, float(0.5 * p0[1, 1] + 0.5 * p1[1, 1] - 0.4)))
    ax1.annotate('p', xy=(float(p0[1, 0] - 0.15), float(p0[1, 1] - 0.15)))

    tf2 = transform.EuclideanTransform(rotation=np.deg2rad(a))
    p2 = tf2(p1)
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

where ``Rotation'`` represents the rotation in the intrinsic frame of reference

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
    p1 = tf1(p0)
    ax1 = plt.subplot(1, 1, 1)
    plt.plot(p0[:,0], p0[:,1], 'ko-')
    plt.plot(p1[:,0], p1[:,1], 'ro-.')
    ax1.annotate('', xy=(float(p1[1, 0]), float(p1[1, 1]-0.2)), xytext=(float(p0[1, 0]), float(p0[1, 1]-0.2)),
                 arrowprops=dict(facecolor='black', shrink=0.1, width=0.1, headlength=5, headwidth=5), )
    ax1.annotate('T', xy=(float(0.5 * p0[1, 0] + 0.5 * p1[1, 0]) + 0.0, float(0.5 * p0[1, 1] + 0.5 * p1[1, 1] - 0.4)))
    ax1.annotate('p', xy=(float(p0[1, 0] + 0.0), float(p0[1, 1] - 0.2)))

    tf = transform.EuclideanTransform(rotation=np.deg2rad(a)) + tf1
    p2 = tf(p0)
    t = np.linspace(np.arctan2(p0[2,1],p0[2,0]), np.arctan2(p2[0,1],p2[0,0]), 100)
    r = np.linalg.norm(p0[0,:])
    plt.plot(p2[:,0], p2[:,1], 'ro-')
    plt.plot(p1[1,0]+r*np.cos(t), p1[1,1]+r*np.sin(t), 'k--')
    ax1.annotate(r'$\alpha$', xy=(float(0.5 * p2[2,0] + 0.5 * p1[2,0]) - 0.15, float(0.5 * p2[2,1] + 0.5 * p1[2,1]) - 0.15))
    ax1.annotate('p*', xy=(float(p2[1, 0] + 0.0), float(p2[1, 1] - 0.2)))

    plt.axis((-1, 4, -1, 2))
    ax1.set_aspect('equal')

    plt.show()


* **First rotation then translation (intrinsic coordinates)**:

.. math::
    Transformation = Translation' * Rotation

where ``Translation'`` represents the translation in the intrinsic frame of reference

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
    p1 = tf1(p0)
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
    p2 = tf(p0)
    plt.plot(p2[:,0], p2[:,1], 'ro-')
    ax1.annotate('', xy=(float(p2[1, 0]), float(p2[1, 1]+0.1)), xytext=(float(p1[1, 0]), float(p1[1, 1]+0.1)),
                 arrowprops=dict(facecolor='black', shrink=0.1, width=0.1, headlength=5, headwidth=5), )
    ax1.annotate('T', xy=(float(0.5 * p2[1, 0] + 0.5 * p1[1, 0] + 0.0), float(0.5 * p2[1, 1] + 0.5 * p1[1, 1] + 0.2)))
    ax1.annotate('p*', xy=(float(p2[1, 0] - 0.0), float(p2[1, 1] - 0.2)))

    plt.axis((-1, 4, -1, 3))
    ax1.set_aspect('equal')

    plt.show()


In code
~~~~~~~

``scikit-image.transform`` can be used to implement these transformations using the following conventions:

* EuclideanTransform(r,t): performs the rotation first, then the translation, in extrinsic coordinates

.. code-block::

    tform = EuclideanTransform(rotation=r, translation=t)
    tform = EuclideanTransform(translation=t, rotation=r)  #order of parameters is irrelevant

* tform1 + tform2: the '+' operator is a magic method that performs
  EuclideanTransform1 first, then EuclideanTransform2, in extrinsic coordinates

.. code-block::

    tform = tform1 + tform2
    tform.matrix = tform2.matrix @ tform1.matrix

Example of ROI placement using rigid transformations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is an example for placing an ROI in the Catphan phantom:

1. ROI placement with respect to nominal phantom:
    1.1. Let's start using an ROI with width = 40 and height = 20

    1.2. Then rotate the ROI by 45 deg

    .. math::
      Tf_1 = R(45°)

    .. code-block::

        tform_1 = tform_r_45 = EuclideanTransform(rotation=np.deg2rad(45))

    1.3. Then translate in the radial direction by 60

    .. math::
      Tf_2 = T'(60) * Tf_1 = Tf_1 * T(60)

    .. code-block::

        tform_t_60 = EuclideanTransform(translation=[60,0])
        tform_2 = tform_t_60 + tform_1

    1.4. Then rotate in place by 90 to align the roi

    .. math::
      Tf_3 = R'(90°) * Tf_2 = Tf_2 * R(90°)

    .. code-block::

        tf_r_90 = EuclideanTransform(rotation=np.deg2rad(90))
        tform_3 = tf_r_90 + tform_2

    1.5. This is the ROI placement with respect to nominal phantom

    .. math::
      Tf_{roi}^{phantom} = Tf_3 = R(45°) * T(60) * R(90°)

    .. code-block::

        tform_roi_phantom = tform_3 = tf_r_90 + tform_t_60 + tform_r_45

.. plot::
    :include-source: False

    import matplotlib.pyplot as plt
    import numpy as np
    from pylinac.core.geometry import Rectangle
    from skimage import transform

    r_phantom = 100
    t = np.linspace(0, 2*np.pi, 100)
    p_phantom = r_phantom * np.vstack((np.sin(t), np.cos(t)))

    width = 50
    height = 20
    angle = 45
    rotation = 90
    radial_distance = 60
    lateral_distance = 0
    rect = Rectangle(width = width, height = height, center=(0,0))
    rect = np.array([v.as_array(("x","y")) for v in rect.vertices])
    rect = np.vstack((rect, rect[0,:]))
    tf1 = transform.EuclideanTransform(rotation=np.deg2rad(angle))
    tf2 = transform.EuclideanTransform(translation=(radial_distance, lateral_distance))
    tf3 = transform.EuclideanTransform(rotation=np.deg2rad(rotation))
    rect_rotated = tf1(rect)            # R
    rect_centered = (tf2 + tf1)(rect)     # T'*R = R*T = T+R
    rect_final = (tf3 + tf2 + tf1)(rect)  # R'(R*T) = R*T*R = R+T+R

    rect_rotated2 = tf3(rect)
    rect_translated = (tf3 + tf2)(rect)

    _, axs = plt.subplots(2, 4)
    axs[0,0].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[0,0].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[0,0].plot(p_phantom[0,:], p_phantom[1,:], 'k', linewidth=2)
    axs[0,0].plot(p_phantom[0,0], p_phantom[1,0], 'ro')
    axs[0,0].plot(rect[:,0], rect[:,1], 'b')
    axs[0,0].axis((-150, 150, -150, 150))
    axs[0,0].set_aspect('equal')

    axs[0,1].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[0,1].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[0,1].plot(p_phantom[0,:], p_phantom[1,:], 'k', linewidth=2)
    axs[0,1].plot(p_phantom[0,0], p_phantom[1,0], 'ro')
    axs[0,1].plot(rect_rotated[:,0], rect_rotated[:,1], 'b')
    axs[0,1].axis((-150, 150, -150, 150))
    axs[0,1].set_aspect('equal')

    axs[0,2].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[0,2].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[0,2].plot(p_phantom[0,:], p_phantom[1,:], 'k', linewidth=2)
    axs[0,2].plot(p_phantom[0,0], p_phantom[1,0], 'ro')
    axs[0,2].plot(rect_centered[:,0], rect_centered[:,1], 'b')
    axs[0,2].axis((-150, 150, -150, 150))
    axs[0,2].set_aspect('equal')

    axs[0,3].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[0,3].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[0,3].plot(p_phantom[0,:], p_phantom[1,:], 'k', linewidth=2)
    axs[0,3].plot(p_phantom[0,0], p_phantom[1,0], 'ro')
    axs[0,3].plot(rect_final[:,0], rect_final[:,1], 'b')
    axs[0,3].axis((-150, 150, -150, 150))
    axs[0,3].set_aspect('equal')

    axs[1,0].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[1,0].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[1,0].plot(p_phantom[0,:], p_phantom[1,:], 'k', linewidth=2)
    axs[1,0].plot(p_phantom[0,0], p_phantom[1,0], 'ro')
    axs[1,0].plot(rect[:,0], rect[:,1], 'b')
    axs[1,0].axis((-150, 150, -150, 150))
    axs[1,0].set_aspect('equal')

    axs[1,1].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[1,1].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[1,1].plot(p_phantom[0,:], p_phantom[1,:], 'k', linewidth=2)
    axs[1,1].plot(p_phantom[0,0], p_phantom[1,0], 'ro')
    axs[1,1].plot(rect_rotated2[:,0], rect_rotated2[:,1], 'b')
    axs[1,1].axis((-150, 150, -150, 150))
    axs[1,1].set_aspect('equal')

    axs[1,2].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[1,2].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[1,2].plot(p_phantom[0,:], p_phantom[1,:], 'k', linewidth=2)
    axs[1,2].plot(p_phantom[0,0], p_phantom[1,0], 'ro')
    axs[1,2].plot(rect_translated[:,0], rect_translated[:,1], 'b')
    axs[1,2].axis((-150, 150, -150, 150))
    axs[1,2].set_aspect('equal')

    axs[1,3].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[1,3].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[1,3].plot(p_phantom[0,:], p_phantom[1,:], 'k', linewidth=2)
    axs[1,3].plot(p_phantom[0,0], p_phantom[1,0], 'ro')
    axs[1,3].plot(rect_final[:,0], rect_final[:,1], 'b')
    axs[1,3].axis((-150, 150, -150, 150))
    axs[1,3].set_aspect('equal')

    axs[0,1].set_title('                          Intrinsic')
    axs[1,1].set_title('                          Extrinsic')

    plt.show()


2. Phantom placement with respect to image (global) coordinates:
    2.1. Let's start with a centered (nominal) phantom

    2.2. Then roll the phantom by 30 deg (exaggerated for visual purposes only)

    .. math::
      Tf_1 = R(30°)

    .. code-block::

        tform_1 = tform_r_30 = EuclideanTransform(rotation=np.deg2rad(30))

    2.3. Then translate the phantom to the image center (150, 150)

    .. math::
      Tf_2 = T(c) * Tf_1

    .. code-block::

        tform_t_c = EuclideanTransform(translation=c)
        tform_2 = tform_1 + tform_t_c

    2.4. This is phantom placement with respect to image coordinates

    .. math::
      Tf_{phantom}^{image} = Tf_2 = T(c) * R(30°)

    .. code-block::

        tform_phantom_image = tform_2 = tform_r_30 + tform_t_c
        # same as tform_phantom_image = EuclideanTransform(rotation=np.deg2rad(30), translation=c)

.. plot::
    :include-source: False

    import matplotlib.pyplot as plt
    import numpy as np
    from skimage import transform

    phantom_roll = 30
    phantom_center = (150,150)
    r_phantom = 100
    t = np.linspace(0, 2*np.pi, 100)
    p_phantom = r_phantom * np.vstack((np.sin(t), np.cos(t)))

    tf1 = transform.EuclideanTransform(rotation=np.deg2rad(phantom_roll))
    tf2 = transform.EuclideanTransform(translation=phantom_center)
    phantom_placement = tf1 + tf2
    phantom_rotated = tf1(p_phantom.T).T
    phantom_final = phantom_placement(p_phantom.T).T

    _, axs = plt.subplots(1, 3)
    axs[0].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[0].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[0].plot(p_phantom[0,:], p_phantom[1,:], 'k', linewidth=2)
    axs[0].plot(p_phantom[0,0], p_phantom[1,0], 'ro')
    axs[0].axis((-150, 300, -150, 300))
    axs[0].set_aspect('equal')

    axs[1].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[1].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[1].plot(phantom_rotated[0,:], phantom_rotated[1,:], 'k', linewidth=2)
    axs[1].plot(phantom_rotated[0,0], phantom_rotated[1,0], 'ro')
    axs[1].axis((-150, 300, -150, 300))
    axs[1].set_aspect('equal')

    axs[2].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[2].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[2].plot(phantom_final[0,:], phantom_final[1,:], 'k', linewidth=2)
    axs[2].plot(phantom_final[0,0], phantom_final[1,0], 'ro')
    axs[2].axis((-150, 300, -150, 300))
    axs[2].set_aspect('equal')

    plt.show()


3. ROI placement with respect to image (global) coordinates:
    3.1. The ROI transformation to global are the cascading transformations

    .. math::
      Tf_{roi}^{image} = Tf_{phantom}^{image} * Tf_{roi}^{phantom}

    .. code-block::

        tform_roi_image = tform_roi_phantom + tform_phantom_image

.. plot::
    :include-source: False

    import matplotlib.pyplot as plt
    import numpy as np
    from pylinac.core.geometry import Rectangle
    from skimage import transform

    phantom_roll = 30
    phantom_center = (150,150)
    r_phantom = 100
    t = np.linspace(0, 2*np.pi, 100)
    p_phantom = r_phantom * np.vstack((np.sin(t), np.cos(t)))

    width = 50
    height = 20
    angle = 45
    rotation = 90
    radial_distance = 60
    lateral_distance = 0
    rect = Rectangle(width = width, height = height, center=(0,0))
    rect = np.array([v.as_array(("x","y")) for v in rect.vertices])
    rect = np.vstack((rect, rect[0,:]))
    tf1 = transform.EuclideanTransform(rotation=np.deg2rad(angle))
    tf2 = transform.EuclideanTransform(translation=(radial_distance, lateral_distance))
    tf3 = transform.EuclideanTransform(rotation=np.deg2rad(rotation))
    roi_placement = tf3 + tf2 + tf1
    rect_phantom = roi_placement(rect)

    tf1 = transform.EuclideanTransform(rotation=np.deg2rad(phantom_roll))
    tf2 = transform.EuclideanTransform(translation=phantom_center)
    phantom_placement = tf1 + tf2
    phantom_final = phantom_placement(p_phantom.T).T

    roi_global = roi_placement + phantom_placement
    rect_final = roi_global(rect)

    _, axs = plt.subplots(1, 2)
    axs[0].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[0].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[0].plot(p_phantom[0,:], p_phantom[1,:], 'k', linewidth=2)
    axs[0].plot(p_phantom[0,0], p_phantom[1,0], 'ro')
    axs[0].plot(rect_phantom[:,0], rect_phantom[:,1], 'b')
    axs[0].axis((-150, 300, -150, 300))
    axs[0].set_aspect('equal')

    axs[1].annotate('', xy=(0, 125), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[1].annotate('', xy=(125, 0), xytext=(0, 0),
                 arrowprops=dict(facecolor='black', shrink=0.0, width=0.1, headlength=5, headwidth=5), )
    axs[1].plot(phantom_final[0,:], phantom_final[1,:], 'k', linewidth=2)
    axs[1].plot(phantom_final[0,0], phantom_final[1,0], 'ro')
    axs[1].plot(rect_final[:,0], rect_final[:,1], 'b')
    axs[1].axis((-150, 300, -150, 300))
    axs[1].set_aspect('equal')

    plt.show()
