.. _sec_internal:

**************
Internal Field
**************

Real case
=========

The internal field potential is given by

.. math::

   V_{int}(r,\theta,\phi) = a \sum_{nm} \left( \frac{a}{r} \right)^{n+1} \left( g_{nm} \cos{m \phi} + h_{nm} \sin{m \phi} \right)
   P_{nm}(\cos{\theta})

Defining :math:`\mathbf{B} = - \nabla V` gives

.. math::

   B_x = -B_{\theta} = \frac{1}{r} \frac{\partial V}{\partial \theta} &= \sum_{nm} \left( \frac{a}{r} \right)^{n+2}
   \left( g_{nm} \cos{m \phi} + h_{nm}\sin{m \phi} \right) \frac{\partial}{\partial \theta} P_{nm}(\cos{\theta}) \\
   B_y = B_{\phi} = -\frac{1}{r \sin{\theta}} \frac{\partial V}{\partial \phi} &= \frac{1}{\sin{\theta}} \sum_{nm}
   \left( \frac{a}{r} \right)^{n+2} m \left( g_{nm} \sin{m \phi} - h_{nm} \cos{m \phi} \right) P_{nm}(\cos{\theta}) \\
   B_z = -B_r = \frac{\partial V}{\partial r} &= -\sum_{nm} (n + 1) \left( \frac{a}{r} \right)^{n+2}
   \left( g_{nm} \cos{m \phi} + h_{nm}\sin{m \phi} \right) P_{nm}(\cos{\theta})

Complex case
============

The complex internal scalar potential is given in geocentric spherical coordinates by

.. math:: V_{int}(r,\theta,\phi) = a \sum_{nm} \left( \frac{a}{r} \right)^{n+1} g_{nm} Y_{nm}(\theta,\phi)

where the coefficients :math:`g_{nm}` are complex and :math:`Y_{nm} = P_{nm}(\cos{\theta}) \exp{im\phi}`.
Defining :math:`\mathbf{B} = - \nabla V` gives

.. math::

   \begin{pmatrix}
     B_x \\
     B_y \\
     B_z
   \end{pmatrix} =
   \sum_{nm} g_{nm} \left( \frac{a}{r} \right)^{n+2}
   \begin{pmatrix}
     \partial_{\theta} Y_{nm} \\
     -\frac{im}{\sin{\theta}} Y_{nm} \\
     -(n+1) Y_{nm}
   \end{pmatrix}
   
..   B_x = -B_{\theta} = \frac{1}{r} \frac{\partial V}{\partial \theta} &= \sum_{nm} \left( \frac{a}{r} \right)^{n+2}
..   g_{nm} \frac{\partial}{\partial \theta} Y_{nm} \\
..   B_y = B_{\phi} = -\frac{1}{r \sin{\theta}} \frac{\partial V}{\partial \phi} &= -\frac{1}{\sin{\theta}} \sum_{nm}
..   im \left( \frac{a}{r} \right)^{n+2} g_{nm} Y_{nm} \\
..   B_z = -B_r = \frac{\partial V}{\partial r} &= -\sum_{nm} (n + 1) \left( \frac{a}{r} \right)^{n+2} g_{nm} Y_{nm}
