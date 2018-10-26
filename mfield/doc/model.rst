**********************
Model Parameterization
**********************

.. |epsiloni| replace:: :math:`\boldsymbol{\epsilon}_i`
.. |deltai| replace:: :math:`\boldsymbol{\delta}_i`
.. |fi| replace:: :math:`f_i`
.. |partialg| replace:: :math:`\frac{\partial}{\partial g_{nm}}`
.. |partialdg| replace:: :math:`\frac{\partial}{\partial \dot{g}_{nm}}`
.. |partialddg| replace:: :math:`\frac{\partial}{\partial \ddot{g}_{nm}}`
.. |partialgp| replace:: :math:`\frac{\partial}{\partial g_{n'm'}}`
.. |partialdgp| replace:: :math:`\frac{\partial}{\partial \dot{g}_{n'm'}}`
.. |partialddgp| replace:: :math:`\frac{\partial}{\partial \ddot{g}_{n'm'}}`
.. |partialeuler| replace:: :math:`\frac{\partial}{\partial \boldsymbol{\alpha}}`
.. |partialeulerp| replace:: :math:`\frac{\partial}{\partial \boldsymbol{\alpha}'}`
.. |partialc| replace:: :math:`\frac{\partial}{\partial \mathbf{c}}`
.. |partialcp| replace:: :math:`\frac{\partial}{\partial \mathbf{c}'}`
.. |partialk| replace:: :math:`\frac{\partial}{\partial k(t)}`
.. |partialkp| replace:: :math:`\frac{\partial}{\partial k(t')}`
.. |depsdg| replace:: :math:`-d\mathbf{B}^{int}_{nm}(\mathbf{r}_i)`
.. |depsdgv| replace:: :math:`-(t_i - t_0) d\mathbf{B}^{int}_{nm}(\mathbf{r}_i)`
.. |depsdga| replace:: :math:`-\frac{1}{2}(t_i - t_0)^2 d\mathbf{B}^{int}_{nm}(\mathbf{r}_i)`
.. |dfdg| replace:: :math:`\frac{1}{|| \mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k})||} \mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k}) \cdot d\mathbf{B}^{int}_{nm}(\mathbf{r}_i)`
.. |dfdgv| replace:: :math:`\frac{t_i - t_0}{|| \mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k})||} \mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k}) \cdot d\mathbf{B}^{int}_{nm}(\mathbf{r}_i)`
.. |dfdga| replace:: :math:`\frac{\frac{1}{2} (t_i - t_0)^2}{|| \mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k})||} \mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k}) \cdot d\mathbf{B}^{int}_{nm}(\mathbf{r}_i)`
.. |dfdc| replace:: :math:`-\frac{1}{F_i(\mathbf{c})} \mathbf{B}^{VFM}_i(\mathbf{c}) \cdot \frac{\partial}{\partial \mathbf{c}} \mathbf{B}^{VFM}_i(\mathbf{c})`
.. |depsdeuler| replace:: :math:`R_q \left[ \frac{\partial}{\partial \boldsymbol{\alpha}} R_3(\boldsymbol{\alpha}) \right] \mathbf{B}^{VFM}_i(\mathbf{c})`
.. |depsdc| replace:: :math:`R_q R_3(\boldsymbol{\alpha}) \frac{\partial}{\partial \mathbf{c}} \mathbf{B}^{VFM}_i(\mathbf{c})`
.. |depsdk| replace:: :math:`-d\mathbf{B}^{ext}(\mathbf{r}_i)`
.. |dfdk| replace:: :math:`\frac{1}{|| \mathbf{B}^{model}(\mathbf{r}_i; \mathbf{g},\mathbf{k})||} \mathbf{B}^{model}(\mathbf{r}_i; \mathbf{g},\mathbf{k}) \cdot d\mathbf{B}^{ext}(\mathbf{r}_i)`
.. |ddepsdeuler| replace:: :math:`R_q \left[ \frac{\partial^2}{\partial \boldsymbol{\alpha}^2} R_3(\boldsymbol{\alpha}) \right] \mathbf{B}^{VFM}_i(\mathbf{c})`
.. |ddfdg| replace:: :math:`\frac{1}{|| \mathbf{B}^{model}(\mathbf{r}_i; \mathbf{g},\mathbf{k})||} \left[ (\mathbf{b}^{model} \cdot d\mathbf{B}^{int}_{nm}(\mathbf{r}_i)) (\mathbf{b}^{model} \cdot d\mathbf{B}^{int}_{n'm'}(\mathbf{r}_i)) + d\mathbf{B}^{int}_{nm}(\mathbf{r}_i) \cdot d\mathbf{B}^{int}_{n'm'}(\mathbf{r}_i) \right]`
.. |ddfdgdgv| replace:: :math:`\frac{t_i - t_0}{|| \mathbf{B}^{model}(\mathbf{r}_i; \mathbf{g},\mathbf{k})||} \left[ (\mathbf{b}^{model} \cdot d\mathbf{B}^{int}_{nm}(\mathbf{r}_i)) (\mathbf{b}^{model} \cdot d\mathbf{B}^{int}_{n'm'}(\mathbf{r}_i)) + d\mathbf{B}^{int}_{nm}(\mathbf{r}_i) \cdot d\mathbf{B}^{int}_{n'm'}(\mathbf{r}_i) \right]`
.. |ddfdgdga| replace:: :math:`\frac{\frac{1}{2} (t_i - t_0)^2}{|| \mathbf{B}^{model}(\mathbf{r}_i; \mathbf{g},\mathbf{k})||} \left[ (\mathbf{b}^{model} \cdot d\mathbf{B}^{int}_{nm}(\mathbf{r}_i)) (\mathbf{b}^{model} \cdot d\mathbf{B}^{int}_{n'm'}(\mathbf{r}_i)) + d\mathbf{B}^{int}_{nm}(\mathbf{r}_i) \cdot d\mathbf{B}^{int}_{n'm'}(\mathbf{r}_i) \right]`
.. |xii| replace:: :math:`\xi_i`
.. |xiiv| replace:: :math:`(t_i - t_0) \xi_i`
.. |xiia| replace:: :math:`\frac{1}{2} (t_i - t_0)^2 \xi_i`
.. |xiivv| replace:: :math:`(t_i - t_0)^2 \xi_i`
.. |xiiva| replace:: :math:`\frac{1}{2} (t_i - t_0)^3 \xi_i`
.. |xiiaa| replace:: :math:`\frac{1}{4} (t_i - t_0)^4 \xi_i`

The penalty function which is minimized is

.. math:: \chi^2 = \sum_{i=1}^{N_{vec}} \boldsymbol{\epsilon}_i \cdot \boldsymbol{\epsilon}_i + \sum_{i=1}^{N_{scal}} f_i^2 + \sum_{i=1}^{N_{vec}^{SV}} \boldsymbol{\delta}_i \cdot \boldsymbol{\delta}_i

where the residuals are

.. math::

   \boldsymbol{\epsilon}_i & = R_q R_3(\boldsymbol{\alpha}) \mathbf{B}^{VFM}_i(\mathbf{c}) - \mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k}) \\
   f_i & = || \mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k}) || - F_i(\mathbf{c}) \\
   \boldsymbol{\delta}_i &= \dot{\mathbf{B}}_i - \dot{\mathbf{B}}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k})

Data measurements
=================

:math:`\mathbf{B}^{VFM}_i(\mathbf{c})` is the vector measurement in the VFM instrument frame which may be calibrated using
the parameters :math:`\mathbf{c} = (\mathbf{s},\mathbf{o},\mathbf{u})`, :math:`F_i(\mathbf{c})` is the calibrated scalar field
measurement, and :math:`\dot{\mathbf{B}}_i` is a vector secular variation measurement (i.e. derived from ground observatory
time series). The fluxgate calibration is as follows,

.. math:: \mathbf{B}^{VFM}_i(\mathbf{c}) = P^{-1}(\mathbf{u}) S(\mathbf{s}) (\mathbf{E}^{VFM}_i - \mathbf{o})

with

.. math::
   
   P(\mathbf{u}) &= \begin{pmatrix}
                      1 & 0 & 0 \\
                      -\sin{u_1} & \cos{u_1} & 0 \\
                      \sin{u_2} & \sin{u_3} & w(u_2,u_3)
                    \end{pmatrix} \\
   P^{-1}(\mathbf{u}) &= \begin{pmatrix}
                           1 & 0 & 0 \\
                           \frac{\sin{u_1}}{\cos{u_1}} & \frac{1}{\cos{u_1}} & 0 \\
                           -\frac{\sin{u_1} \sin{u_3} + \cos{u_1} \sin{u_2}}{w \cos{u_1}} & -\frac{\sin{u_3}}{w \cos{u_1}} & \frac{1}{w}
                         \end{pmatrix} \\
   w(u_2,u_3) &= \sqrt{1 - \sin^2{u_2} - \sin^2{u_3}} \\
   S(\mathbf{s}) &= \textrm{diag}(s_1, s_2, s_3) \\
   \mathbf{o} &= \left[ o_1, o_2, o_3 \right]^T

and :math:`\mathbf{E}^{VFM}_i` is the uncalibrated vector measurement in the VFM frame.
The matrix :math:`R_3(\boldsymbol{\alpha})` rotates a vector from the VFM frame to the CRF frame defined
by the star camera using the Euler angles :math:`\boldsymbol{\alpha} = (\alpha,\beta,\gamma)`. The matrix :math:`R_q` then rotates from CRF to NEC
(see Olsen et al, 2013).

The calibration parameters are given a time dependence represented by a B-spline,

.. math:: \mathbf{c}(t) = \sum_i \mathbf{c}_i N_{i,k}(t) =
                          \sum_i \begin{pmatrix}
                                   \mathbf{s}_i \\
                                   \mathbf{o}_i \\
                                   \mathbf{u}_i
                                 \end{pmatrix} N_{i,k}(t)

The basis splines :math:`N_{i,k}(t)` are defined by a knot vector which is constructed using
uniformly spaced knots over the data time period with a spacing of 30 days.

Model
=====

The vector model is given by

.. math:: \mathbf{B}^{model}(\mathbf{r}, t; \mathbf{g},\mathbf{k}) = \mathbf{B}^{int}(\mathbf{r}, t; \mathbf{g}) + \mathbf{B}^{crust}(\mathbf{r}) + \mathbf{B}^{ext}(\mathbf{r}, t) + \mathbf{B}^{ext,correction}(\mathbf{r}, t; \mathbf{k})

where :math:`\mathbf{B}^{int}(\mathbf{r}, t; \mathbf{g})` is an internal field model defined in :ref:`sec_internal`,
:math:`\mathbf{B}^{crust}(\mathbf{r})` is the MF7 crustal field model from degree 16 to 133,
:math:`\mathbf{B}^{ext}(\mathbf{r}, t)` is an external field model (POMME or CHAOS), and :math:`\mathbf{B}^{ext,correction}(\mathbf{r}, t; \mathbf{k})` is a daily
ring current correction, parameterized as

.. math:: \mathbf{B}^{ext,correction}(\mathbf{r}, t; \mathbf{k}) = k(t) \left( 0.7 \mathbf{B}^{ext,dipole}(\mathbf{r}) + 0.3 \mathbf{B}^{int,dipole}(\mathbf{r}) \right)

where :math:`\mathbf{B}^{ext,dipole}` and :math:`\mathbf{B}^{int,dipole}` are degree 1 external and internal dipole fields,
whose coefficients are aligned with the main field dipole (ie: :math:`g_{10},g_{11},h_{11}`). For each day, there is
one coefficient :math:`k(t)` for that day which is determined through the least
squares minimization. The crustal field term, :math:`\mathbf{B}^{crust}`, can
optionally be set to zero, in order to fit a high degree crustal field
in :math:`\mathbf{B}^{int}`.

The time derivative of the model vector is:

.. math:: \dot{\mathbf{B}}^{model}(\mathbf{r}, t; \mathbf{g},\mathbf{k}) = \dot{\mathbf{B}}^{int}(\mathbf{r}, t; \mathbf{g}) + \dot{\mathbf{B}}^{ext,correction}(\mathbf{r}, t; \mathbf{k})

We don't include :math:`\dot{\mathbf{B}}^{ext}(\mathbf{r}, t)` since an external field model is removed from the observatory
time series prior to calculating the secular variation measurements. The internal field model and its time derivative are
given below.

.. math::

   \mathbf{B}^{int}(\mathbf{r}, t; \mathbf{g}) &=
   \begin{pmatrix}
     B_x \\
     B_y \\
     B_z
   \end{pmatrix} =
   \sum_{n=1}^N \sum_{m=-n}^n \tilde{g}_n^m(t) \left( \frac{a}{r} \right)^{n+2}
   \begin{pmatrix}
     \partial_{\theta} S_n^m \\
     -\partial_{\phi} S_n^m \\
     -(n+1) S_n^m
   \end{pmatrix} \\
   \dot{\mathbf{B}}^{int}(\mathbf{r}, t; \mathbf{g}) &=
   \begin{pmatrix}
     \dot{B}_x \\
     \dot{B}_y \\
     \dot{B}_z
   \end{pmatrix} =
   \sum_{n=1}^N \sum_{m=-n}^n \dot{\tilde{g}}_n^m(t) \left( \frac{a}{r} \right)^{n+2}
   \begin{pmatrix}
     \partial_{\theta} S_n^m \\
     -\partial_{\phi} S_n^m \\
     -(n+1) S_n^m
   \end{pmatrix} = \mathbf{B}^{int}(\mathbf{r}, t; \dot{\mathbf{g}})

Jacobian
========

When minimizing :math:`\chi^2` with a nonlinear least squares algorithm, the Jacobian
is required.
For easy reference, we list the derivatives of the residuals with respect
to various model parameters, needed for the Jacobian calculation.

Internal field
--------------

The internal field model can be expressed as

.. math:: \mathbf{B}^{int}(\mathbf{r}, t; \mathbf{g}) = \sum_{nm} g_{nm}(t) d\mathbf{B}^{int}_{nm}(\mathbf{r})

where

.. math::

   d\mathbf{B}^{int}_{nm}(\mathbf{r}) =
   \left\{
   \begin{array}{cc}
   \left( \frac{a}{r} \right)^{n+2}
   \left(
   \begin{array}{c}
   \cos{(m\phi)} \partial_{\theta} P_{nm} \\
   \frac{m}{\sin{\theta}} \sin{(m\phi)} P_{nm} \\
   -(n+1) \cos{(m\phi)} P_{nm} \\
   \end{array}
   \right) & m \ge 0 \\
   \left( \frac{a}{r} \right)^{n+2}
   \left(
   \begin{array}{c}
   \sin{(m\phi)} \partial_{\theta} P_{nm} \\
   -\frac{m}{\sin{\theta}} \cos{(m\phi)} P_{nm} \\
   -(n+1) \sin{(m\phi)} P_{nm}
   \end{array}
   \right) & m < 0
   \end{array}
   \right.

Fluxgate calibration
--------------------

Let

.. math::
   
   s_j(t) &= \sum_k s_{jk} N_k(t) \\
   o_j(t) &= \sum_k o_{jk} N_k(t) \\
   u_j(t) &= \sum_k u_{jk} N_k(t)

where :math:`j = 1,2,3` and :math:`k` is summed from :math:`1` to the number of
control points in each spline. Then,

.. math::

   \frac{\partial}{\partial s_{jk}} \mathbf{B}^{VFM}_i(\mathbf{c}) &= N_k(t_i) \left( (\mathbf{E}^{VFM}_i)_j - o_j(t_i) \right) P^{-1}_j(\mathbf{u}(t_i)) \\
   \frac{\partial}{\partial o_{jk}} \mathbf{B}^{VFM}_i(\mathbf{c}) &= -N_k(t_i) s_j(t_i) P^{-1}_j(\mathbf{u}(t_i)) \\
   \frac{\partial}{\partial u_{jk}} \mathbf{B}^{VFM}_i(\mathbf{c}) &= N_k(t_i) \left[ \frac{\partial}{\partial u_j} P^{-1}(\mathbf{u}(t_i))\right] S(\mathbf{s}(t_i)) \left( \mathbf{E}^{VFM}_i - \mathbf{o}(t_i) \right)

where :math:`P^{-1}_j` is the :math:`j`-th column of :math:`P^{-1}`. Then,

.. math::
   
   \frac{\partial \boldsymbol{\epsilon}_i}{\partial ( \cdot )} &= R_q R_3(\boldsymbol{\alpha}) \frac{\partial}{\partial (\cdot)} \mathbf{B}^{VFM}_i(\mathbf{c}) \\
   \frac{\partial f_i}{\partial (\cdot)} &= -\frac{1}{F_i(\mathbf{c})} \mathbf{B}^{VFM}_i(\mathbf{c}) \cdot \frac{\partial}{\partial (\cdot)} \mathbf{B}^{VFM}_i(\mathbf{c}) \\
   \frac{\partial \boldsymbol{\delta}_i}{\partial ( \cdot )} &= 0

where :math:`(\cdot)` refers to one of :math:`s_{jk},o_{jk},u_{jk}`.

First derivatives
-----------------

The following table summarizes the first derivatives of the residuals needed for the Jacobian.

============== ========================== =========================== ========================
Derivative     Vector residual |epsiloni| Scalar residual :math:`f_i` Vector residual |deltai|
============== ========================== =========================== ========================
|partialg|     |depsdg|                   |dfdg|                      0
|partialdg|    |depsdgv|                  |dfdgv|                     |depsdg|
|partialddg|   |depsdga|                  |dfdga|                     |depsdgv|
|partialeuler| |depsdeuler|               0                           0
|partialc|     |depsdc|                   |dfdc|                      0
|partialk|     |depsdk|                   |dfdk|
============== ========================== =========================== ========================

Second derivatives
------------------

To use the geodesic acceleration method, we also need the second derivatives, given in
the tables below. For the vector residuals, we have

=============== ========== ============== ========== ==========
|epsiloni|      |partialg| |partialeuler| |partialc| |partialk|
=============== ========== ============== ========== ==========
|partialgp|     0          0              0          0
|partialeulerp| 0          |ddepsdeuler|  X          0
|partialcp|     0          X              X          0
|partialkp|     0          0              0          0
=============== ========== ============== ========== ==========

Therefore, the second directional derivative of the vector residual |epsiloni| is

.. math::

  D_v^2 \boldsymbol{\epsilon_i} = R_q
  \left[
    v_{\alpha}^2 \partial^2_{\alpha} + v_{\beta}^2 \partial^2_{\beta} + v_{\gamma}^2 \partial^2_{\gamma} +
    2 v_{\alpha} v_{\beta} \partial_{\alpha} \partial_{\beta} +
    2 v_{\alpha} v_{\gamma} \partial_{\alpha} \partial_{\gamma} +
    2 v_{\beta} v_{\gamma} \partial_{\beta} \partial_{\gamma}
  \right]
  R_3(\boldsymbol{\alpha}) \mathbf{B}^{VFM}_i

For the scalar residuals, we have

=============== ========== =========== ============ ============== ==========
|fi|            |partialg| |partialdg| |partialddg| |partialeuler| |partialk|
=============== ========== =========== ============ ============== ==========
|partialgp|     |xii|      |xiiv|      |xiia|       0              X
|partialdgp|    |xiiv|     |xiivv|     |xiiva|      0              X
|partialddgp|   |xiia|     |xiiva|     |xiiaa|      0              X
|partialeulerp| 0          0           0            0              0
|partialkp|     X          X           X            0              X
=============== ========== =========== ============ ============== ==========

In the above table,

.. math:: \xi_i = \frac{\partial^2 f_i}{\partial g_{nm} \partial g_{n'm'}} = \frac{1}{|| \mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k})||} \left[ (\mathbf{b}^{model} \cdot d\mathbf{B}^{int}_{nm}(\mathbf{r}_i)) (\mathbf{b}^{model} \cdot d\mathbf{B}^{int}_{n'm'}(\mathbf{r}_i)) + d\mathbf{B}^{int}_{nm}(\mathbf{r}_i) \cdot d\mathbf{B}^{int}_{n'm'}(\mathbf{r}_i) \right]

and

.. math:: \mathbf{b}^{model}(\mathbf{r}_i, t_i; \mathbf{g}, \mathbf{k}) = \frac{\mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g}, \mathbf{k})}{|| \mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g}, \mathbf{k}) || }

Therefore, the second directional derivative of the scalar residual |fi| is

.. math::

  D_v^2 f_i = \sum_{nm,n'm'} \xi_{i,nm,n'm'}
  & \left[
    v_{nm}^{MF} v_{n'm'}^{MF} + (t_i-t_0)^2 v_{nm}^{SV} v_{n'm'}^{SV} + \frac{1}{4} (t_i-t_0)^4 v_{nm}^{SA} v_{n'm'}^{SA} +
    \right. \\
  & \left.
    (t_i-t_0) (v_{nm}^{MF} v_{n'm'}^{SV} + v_{nm}^{SV} v_{n'm'}^{MF}) +
    \right. \\
  & \left.
    \left( \frac{1}{2} (t_i-t_0)^2 \right) (v_{nm}^{MF} v_{n'm'}^{SA} + v_{nm}^{SA} v_{n'm'}^{MF}) +
    \right. \\
  & \left.
    \left( \frac{1}{2} (t_i-t_0)^3 \right) (v_{nm}^{SV} v_{n'm'}^{SA} + v_{nm}^{SA} v_{n'm'}^{SV})
    \right]

Optimization
------------

Since the cost function :math:`\chi^2` depends on both vector and scalar
residuals, we can write the Jacobian as

.. math::

   \mathbf{J} =
   \left(
   \begin{array}{ccccc}
   \mathbf{J}_{MF}^{vec} & \mathbf{J}_{SV}^{vec} & \mathbf{J}_{SA}^{vec} & \mathbf{J}_{Euler}^{vec}(\mathbf{x}) & \mathbf{J}_{ext}^{vec}(\mathbf{x}) \\
   \mathbf{J}_{MF}^{scal}(\mathbf{x}) & \mathbf{J}^{scal}_{SV}(\mathbf{x}) & \mathbf{J}^{scal}_{SA}(\mathbf{x}) & 0 & \mathbf{J}^{scal}_{ext}(\mathbf{x}) \\
   0 & \dot{\mathbf{J}}_{SV}^{vec} & \dot{\mathbf{J}}_{SA}^{vec} & 0 & 0
   \end{array}
   \right)

where the top portion corresponds to vector residuals |epsiloni|, the middle portion
corresponds to scalar residuals :math:`f_i`, and the bottom portion to |deltai|.
Even if the vector and scalar residuals are "mixed", so that
the Jacobian does not separate vertically as shown above, we can consider the above matrix
without loss of generality, since we can always rearrange the rows of the matrix as needed.
For simplicity, we define

.. math::

   \mathbf{J}_{int} =
   \left(
   \begin{array}{ccc}
   \mathbf{J}_{MF} & \mathbf{J}_{SV} & \mathbf{J}_{SA}
   \end{array}
   \right)

and note :math:`\mathbf{J}_{SV} = t \mathbf{J}_{MF}` and :math:`\mathbf{J}_{SA} = \frac{1}{2} t^2 \mathbf{J}_{MF}`,
where :math:`t` is the timestamp of measurement :math:`i`. The Jacobian then becomes

.. math::

   \mathbf{J} =
   \left(
   \begin{array}{ccc}
   \mathbf{J}_{int}^{vec} & \mathbf{J}_{Euler}^{vec}(\mathbf{x}) & \mathbf{J}_{ext}^{vec}(\mathbf{x}) \\
   \mathbf{J}_{int}^{scal}(\mathbf{x}) & 0 & \mathbf{J}^{scal}_{ext}(\mathbf{x}) \\
   \dot{\mathbf{J}}_{int}^{vec} & 0 & 0
   \end{array}
   \right)

Note that for vector residuals, :math:`\mathbf{J}_{int}` does not depend on the model parameters
:math:`\mathbf{x}`. Also, the scalar residuals do not depend on the Euler angles, resulting in the
block of zeros in the above matrix. Additionally, while the matrices :math:`\mathbf{J}_{int}^{vec}`
and :math:`\mathbf{J}_{int}^{scal}(\mathbf{x})` are dense, the rest of the Jacobian corresponding
to the Euler angles and external field parameters has a lot of sparse structure.
During the nonlinear least squares iterations, we require the normal equations matrix
:math:`\mathbf{J}^T \mathbf{J}`. This matrix can be computed very efficiently by accounting
for the sparse structure in the above matrix. Writing it all out, we have:

.. math::

   \mathbf{J}^T \mathbf{J} =
   \left(
   \begin{array}{ccc}
   \mathbf{J}_{int}^T \mathbf{J}_{int}^{vec} + \mathbf{J}_{int}^T(\mathbf{x}) \mathbf{J}_{int}^{scal}(\mathbf{x}) + \dot{\mathbf{J}}_{int}^T \dot{\mathbf{J}}_{int}^{vec} & X & X \\
   \mathbf{J}_{Euler}^T(\mathbf{x}) \mathbf{J}_{int}^{vec} & \mathbf{J}_{Euler}^T(\mathbf{x}) \mathbf{J}_{Euler}(\mathbf{x}) & X \\
   \mathbf{J}_{ext}^T(\mathbf{x}) \mathbf{J}_{int}^{vec} + \mathbf{J}_{ext}^T(\mathbf{x}) \mathbf{J}_{int}^{scal}(\mathbf{x})  & \mathbf{J}_{ext}^T(\mathbf{x}) \mathbf{J}_{Euler}(\mathbf{x}) & \mathbf{J}_{ext}^T(\mathbf{x}) \mathbf{J}_{ext}^{vec}(\mathbf{x}) + \mathbf{J}_{ext}^T(\mathbf{x}) \mathbf{J}_{ext}^{scal}(\mathbf{x}) \\
   \end{array}
   \right)

The :math:`X` entries above indicate that the matrix is symmetric and so only the lower half needs to
be computed. The :math:`(1,1)` term :math:`\mathbf{J}_{int}^{T,vec} \mathbf{J}_{int}^{vec} + \dot{\mathbf{J}}_{int}^{T,vec} \dot{\mathbf{J}}_{int}^{vec}`
can be precomputed since it does not depend on :math:`\mathbf{x}`, which saves significant computations during the
iteration.
