************************
Spherical Harmonic Notes
************************

Here are some notes on spherical harmonic manipulations to avoid
needing to work them out each time they are needed.

Real and Complex Expansions
===========================

First, we will show the relationship between the complex and real
spherical harmonic representation of a real-valued function.
Consider the complex spherical harmonic decomposition of a
real-valued function:

.. math:: f(r,\theta,\phi) = \sum_{n=1}^N \sum_{m=-n}^n c_n(r) f_n^m Y_n^{m}(\theta,\phi)

where :math:`f_n^m` are complex coefficients, :math:`c_n(r)` is an arbitrary function of :math:`r` and

.. math:: Y_n^m(\theta,\phi) = P_n^{|m|}(\cos{\theta}) e^{im\phi}

Here, the associated Legendre function :math:`P_n^{|m|}(\cos{\theta})` can have any
desired normalization convention. For the magnetic potential :math:`V(r,\theta,\phi)`
of an internal source, the radial function would take the form

.. math:: c_n(r) = \left( \frac{a}{r} \right)^{n+1}

for some reference radius :math:`a`. For a magnetic field component, it would be of the
form

.. math:: c_n(r) = \left( \frac{a}{r} \right)^{n+2}

Now we wish to find the relationship between the representation of
:math:`f(r,\theta,\phi)` with complex coefficients :math:`f_n^m`, and its
representation in terms of real Gauss coefficients :math:`g_n^m` and :math:`h_n^m`.
First, because :math:`f(r,\theta,\phi)` is assumed real-valued (i.e., it is a
magnetic scalar potential function or a magnetic field component), this will
introduce a symmetry constraint on the :math:`f_n^m`. To see this, take the
complex conjugate of :math:`f(r,\theta,\phi)`:

.. math::

   f(r,\theta,\phi) &= f^*(r,\theta,\phi) \\
                    &= \sum_{n=1}^N \sum_{m=-n}^n c_n(r) (f_n^m)^* P_n^{|m|}(\cos{\theta}) e^{-im\phi} \\
                    &= \sum_{n=1}^N \sum_{m=n}^{-n} c_n(r) (f_n^{-m})^* P_n^{|-m|}(\cos{\theta}) e^{+im\phi} \\
   \sum_{n=1}^N \sum_{m=-n}^n c_n(r) f_n^m P_n^{|m|}(\cos{\theta}) e^{im\phi} &= \sum_{n=1}^N \sum_{m=-n}^n c_n(r) (f_n^{-m})^* P_n^{|m|}(\cos{\theta}) e^{im\phi}

In the third line we relabed the index :math:`m \rightarrow -m`. This implies

.. math:: f_n^{-m} = (f_n^m)^*

which is true for any real-valued Fourier transform. Now, we can relate the
complex coefficients to the real-valued Gauss coefficients :math:`g_n^m` and
:math:`h_n^m` as follows.

.. math::

   f(r,\theta,\phi) &= \sum_{n=1}^N \sum_{m=-n}^n c_n(r) f_n^m P_n^{|m|}(\cos{\theta}) e^{im\phi} \\
                    &= \sum_{n=1}^N c_n(r) \left( \sum_{m=-n}^{-1} f_n^m P_n^{|m|}(\cos{\theta}) e^{im\phi} + \sum_{m=0}^n f_n^m P_n^{|m|}(\cos{\theta}) e^{im\phi} \right) \\
                    &= \sum_{n=1}^N c_n(r) \left( \sum_{m=n}^{1} f_n^{-m} P_n^{|-m|}(\cos{\theta}) e^{-im\phi} + \sum_{m=0}^n f_n^m P_n^{|m|}(\cos{\theta}) e^{im\phi} \right) \\
                    &= \sum_{n=1}^N c_n(r) \left( \sum_{m=1}^n f_n^{-m} P_n^{|m|}(\cos{\theta}) e^{-im\phi} + \sum_{m=0}^n f_n^m P_n^{|m|}(\cos{\theta}) e^{im\phi} \right) \\
                    &= \sum_{n=1}^N c_n(r) \left( \sum_{m=1}^n (f_n^m)^* P_n^{|m|}(\cos{\theta}) e^{-im\phi} + \sum_{m=0}^n f_n^m P_n^{|m|}(\cos{\theta}) e^{im\phi} \right) \\
                    &= \sum_{n=1}^N c_n(r) \left( f_n^0 P_n^0 + \sum_{m=1}^n P_n^{|m|}(\cos{\theta}) \left\{ (f_n^m)^* e^{-im\phi} + f_n^m e^{im\phi} \right\} \right)

In the third line we again used the mapping :math:`m \rightarrow -m` in the first sum. In the fifth
line we used our symmetry property :math:`f_n^{-m} = (f_n^m)^*`. Now we will define our complex
coefficients :math:`f_n^m` in terms of real-valued Gauss coefficients as,

.. math::

   f_n^m = \left\{
             \begin{array}{cc}
               g_n^0 & m = 0 \\
               \frac{1}{2} \left( g_n^m - i h_n^m \right) & m > 0 \\
               \frac{1}{2} \left( g_n^{|m|} + i h_n^{|m|} \right) & m < 0
             \end{array}
           \right.

The reason for the factor of :math:`1/2` will become clear soon. Also note that
the Gauss coefficients :math:`g_n^m` are defined only for :math:`m \ge 0`,
with :math:`h_n^m` defined for :math:`m > 0`.
This definition of :math:`f_n^m` also satisfies our symmetry condition
:math:`f_n^{-m} = (f_n^m)^*`. Continuing our above derivation, we find,

.. math::

   f(r,\theta,\phi) &= \sum_{n=1}^N c_n(r) \left( f_n^0 P_n^0 + \sum_{m=1}^n P_n^{|m|}(\cos{\theta}) \left\{ (f_n^m)^* e^{-im\phi} + f_n^m e^{im\phi} \right\} \right) \\
                    &= \sum_{n=1}^N c_n(r) \left( f_n^0 P_n^0 + \sum_{m=1}^n P_n^{|m|}(\cos{\theta}) \left\{ \frac{1}{2} (g_n^m + i h_n^m) (\cos{m\phi} - i \sin{m \phi}) + \frac{1}{2} (g_n^m - i h_n^m) (\cos{m\phi} + i \sin{m\phi}) \right\} \right) \\
                    &= \sum_{n=1}^N c_n(r) \left( f_n^0 P_n^0 + \frac{1}{2} \sum_{m=1}^n P_n^{|m|}(\cos{\theta}) \left\{ 2 g_n^m \cos{m\phi} + 2 h_n^m \sin{m\phi} \right\} \right) \\
                    &= \sum_{n=1}^N \sum_{m=0}^nc_n(r) \left( g_n^m \cos{m\phi} + h_n^m \sin{m\phi} \right) P_n^m(\cos{\theta})

This last expression can be recognized as the standard internal magnetic scalar potential field
with the usual Gauss coefficients :math:`g_n^m,h_n^m`.

To summarize,

.. important::

   The complex and real spherical harmonic representations of a real-valued function are
   related by,

   .. math::

      \sum_{n=1}^N \sum_{m=-n}^n c_n(r) f_n^m P_n^{|m|}(\cos{\theta}) e^{im\phi} = \sum_{n=1}^N \sum_{m=0}^nc_n(r) \left( g_n^m \cos{m\phi} + h_n^m \sin{m\phi} \right) P_n^m(\cos{\theta})

   and the coefficients are related as

   .. math::

      f_n^m = \left\{
                \begin{array}{cc}
                  g_n^0 & m = 0 \\
                  \frac{1}{2} \left( g_n^m - i h_n^m \right) & m > 0 \\
                  \frac{1}{2} \left( g_n^{|m|} + i h_n^{|m|} \right) & m < 0
                \end{array}
              \right.

.. _sec_internal:

Internal Field
==============

Real case
---------

The internal field potential is given by

.. math::

   V_{int}(r,\theta,\phi) = a \sum_{n=1}^N \sum_{m=0}^n \left( \frac{a}{r} \right)^{n+1} \left( g_n^m \cos{m \phi} + h_n^m \sin{m \phi} \right)
   P_n^m(\cos{\theta})

If we define

.. math::

   \tilde{g}_n^m =
     \left\{
       \begin{array}{cc}
         g_n^m & m \ge 0 \\
         h_n^{|m|} & m < 0
       \end{array}
     \right.

and

.. math::

   S_n^m(\theta,\phi) =
     \left\{
       \begin{array}{cc}
         \cos{(m\phi)} P_n^m(\cos{\theta}) & m \ge 0 \\
         \sin{(|m|\phi)} P_n^{|m|}(\cos{\theta}) & m < 0
       \end{array}
     \right.

then we can write :math:`V_{int}(r,\theta,\phi)` more compactly as

.. math:: V_{int}(r,\theta,\phi) = \sum_{n=1}^N \sum_{m=-n}^n \tilde{g}_n^m V_n^m(r,\theta,\phi)

with

.. math:: V_n^m(r,\theta,\phi) = a \left( \frac{a}{r} \right)^{n+1} S_n^m(\theta,\phi)

Defining :math:`\mathbf{B} = - \nabla V` gives

.. math:: \mathbf{B}(r,\theta,\phi) = \sum_{n=1}^N \sum_{m=-n}^n \tilde{g}_n^m \mathbf{B}_n^m(r,\theta,\phi)

where :math:`\mathbf{B}_n^m = -\nabla V_n^m`. Writing out components, we have

.. math::

   \begin{pmatrix}
     B_r \\
     B_{\theta} \\
     B_{\phi}
   \end{pmatrix} =
   \sum_{n=1}^N \sum_{m=-n}^n \tilde{g}_n^m \left( \frac{a}{r} \right)^{n+2}
   \begin{pmatrix}
     (n+1) S_n^m \\
     -\partial_{\theta} S_n^m \\
     -\partial_{\phi} S_n^m
   \end{pmatrix}

or

.. math::

   \begin{pmatrix}
     B_x \\
     B_y \\
     B_z
   \end{pmatrix} =
   \sum_{n=1}^N \sum_{m=-n}^n \tilde{g}_n^m \left( \frac{a}{r} \right)^{n+2}
   \begin{pmatrix}
     \partial_{\theta} S_n^m \\
     -\partial_{\phi} S_n^m \\
     -(n+1) S_n^m
   \end{pmatrix}

Writing this out explicitly in terms of :math:`g_n^m,h_n^m` gives

.. math::

   B_x &= \sum_{n=1}^N \sum_{m=0}^n \left( \frac{a}{r} \right)^{n+2}
   \left( g_n^m \cos{m \phi} + h_n^m\sin{m \phi} \right) \frac{\partial}{\partial \theta} P_n^m(\cos{\theta}) \\
   B_y &= \frac{1}{\sin{\theta}} \sum_{n=1}^N \sum_{m=0}^n
   \left( \frac{a}{r} \right)^{n+2} m \left( g_n^m \sin{m \phi} - h_n^m \cos{m \phi} \right) P_n^m(\cos{\theta}) \\
   B_z &= -\sum_{n=1}^N \sum_{m=0}^n (n + 1) \left( \frac{a}{r} \right)^{n+2}
   \left( g_n^m \cos{m \phi} + h_n^m\sin{m \phi} \right) P_n^m(\cos{\theta})

Complex case
------------

The complex internal scalar potential is given in geocentric spherical coordinates by

.. math:: V_{int}(r,\theta,\phi) = a \sum_{n=1}^N \sum_{m=-n}^n \left( \frac{a}{r} \right)^{n+1} f_n^m Y_n^m(\theta,\phi)

where the coefficients :math:`f_n^m` are complex and

.. math:: Y_n^m = P_n^{|m|}(\cos{\theta}) e^{im\phi}

Defining :math:`\mathbf{B} = - \nabla V` gives

.. math::

   \begin{pmatrix}
     B_r \\
     B_{\theta} \\
     B_{\phi}
   \end{pmatrix} =
   \sum_{nm} f_n^m \left( \frac{a}{r} \right)^{n+2}
   \begin{pmatrix}
     (n+1) Y_n^m \\
     -\partial_{\theta} Y_n^m \\
     -\frac{im}{\sin{\theta}} Y_n^m
   \end{pmatrix}

or,

.. math::

   \begin{pmatrix}
     B_x \\
     B_y \\
     B_z
   \end{pmatrix} =
   \sum_{nm} f_n^m \left( \frac{a}{r} \right)^{n+2}
   \begin{pmatrix}
     \partial_{\theta} Y_n^m \\
     -\frac{im}{\sin{\theta}} Y_n^m \\
     -(n+1) Y_n^m
   \end{pmatrix}

Useful Identities
=================

The following identities are useful for developing regularization conditions in inverse problems.

.. math::
   
   \int_0^{2\pi} e^{i(m - m')\phi} d\phi &= 2 \pi \delta_{m m'} \\
   \int S_n^m(\theta,\phi) S_{n'}^{m'}(\theta,\phi) d\Omega &= \delta_{m m'} \delta_{n n'} \frac{4 \pi}{2n + 1}
