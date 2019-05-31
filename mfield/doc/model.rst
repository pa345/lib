**********************
Model Parameterization
**********************

.. |epsiloni| replace:: :math:`\boldsymbol{\epsilon}_i`
.. |deltai| replace:: :math:`\boldsymbol{\delta}_i`
.. |fi| replace:: :math:`f_i`
.. |partialg| replace:: :math:`\frac{\partial}{\partial g_{n,k}^m}`
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
.. |depsdg| replace:: :math:`-N_k(t_i) \mathbf{B}_n^m(\mathbf{r}_i)`
.. |dfdg| replace:: :math:`-N_k(t_i) \mathbf{b}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k}) \cdot \mathbf{B}_n^m(\mathbf{r}_i)`
.. |dfdc| replace:: :math:`\frac{1}{F_i(\mathbf{c})} \mathbf{B}^{VFM}_i(\mathbf{c}) \cdot \frac{\partial}{\partial \mathbf{c}} \mathbf{B}^{VFM}_i(\mathbf{c})`
.. |ddeltadg| replace:: :math:`-\dot{N}_k(t_i) \mathbf{B}_n^m(\mathbf{r}_i)`
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
   f_i &= F_i(\mathbf{c}) - || \mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k}) || \\
   \boldsymbol{\delta}_i &= \dot{\mathbf{B}}_i - \dot{\mathbf{B}}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k})

Data measurements
=================

:math:`\mathbf{B}^{VFM}_i(\mathbf{c})` is the vector measurement in the VFM instrument frame which may be calibrated using
the parameters :math:`\mathbf{c} = (\mathbf{s},\mathbf{o},\mathbf{u})`, :math:`F_i(\mathbf{c}) = || \mathbf{B}^{VFM}(\mathbf{c}) ||`
is the calibrated scalar field measurement, and :math:`\dot{\mathbf{B}}_i` is a vector secular variation measurement (i.e. derived from ground observatory
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

The Euler angles are given a time dependence represented by a B-spline,

.. math:: \boldsymbol{\alpha}(t) = \sum_i \boldsymbol{\alpha}_i N_{i,k}(t) =
                                   \sum_i \begin{pmatrix}
                                            \alpha_i \\
                                            \beta_i \\
                                            \gamma_i
                                          \end{pmatrix} N_{i,k}(t)

The calibration parameters are also given a time dependence represented by a B-spline,

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
time series prior to calculating the secular variation measurements.

Internal Field Model
--------------------

The internal field model and its time derivative are
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

where

.. math::

   \tilde{g}_n^m(t) =
     \left\{
       \begin{array}{cc}
         g_n^m(t) & m \ge 0 \\
         h_n^{|m|}(t) & m < 0
       \end{array}
     \right.

and :math:`g_n^m(t), h_n^m(t)` are parameterized with B-splines:

.. math::

   g_n^m(t) &= \sum_i g_{n,i}^m N_{i,k}(t) \\
   h_n^m(t) &= \sum_i h_{n,i}^m N_{i,k}(t)

Jacobian
========

When minimizing :math:`\chi^2` with a nonlinear least squares algorithm, the Jacobian
is required.
For easy reference, we list the derivatives of the residuals with respect
to various model parameters, needed for the Jacobian calculation.

Internal field
--------------

The internal field model can be expressed as

.. math:: \mathbf{B}^{int}(\mathbf{r}, t; \mathbf{g}) = \sum_{nm} g_n^m(t) \mathbf{B}_n^m(\mathbf{r}) = \sum_{nmk} g_{n,k}^m N_k(t) \mathbf{B}_n^m(\mathbf{r})

where

.. math::

   \mathbf{B}_n^m(\mathbf{r}) =
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

Then,

.. math::

   \frac{\partial \boldsymbol{\epsilon}_i}{\partial g_{n,k}^m} &= -N_k(t_i) \mathbf{B}_n^m(\mathbf{r}_i) \\
   \frac{\partial f_i}{\partial g_{n,k}^m} &= -N_k(t_i) \frac{\mathbf{B}^{model}(\mathbf{r}_i,t_i;\mathbf{g},\mathbf{k})}{|| \mathbf{B}^{model}(\mathbf{r}_i,t_i;\mathbf{g},\mathbf{k}) ||} \cdot \mathbf{B}_n^m(\mathbf{r}_i) \\
   \frac{\partial \boldsymbol{\delta}_i}{\partial g_{n,k}^m} &= -\dot{N}_k(t_i) \mathbf{B}_n^m(\mathbf{r}_i)

Euler angles
------------

Let

.. math::
   
   \alpha_j(t) &= \sum_k \alpha_{jk} N_k(t) \\

where :math:`j = 1,2,3` and :math:`k` is summed from :math:`1` to the number of
control points in each spline. Then,

.. math::

   \frac{\partial \boldsymbol{\epsilon}_i}{\partial \alpha_{jk}} &= N_k(t_i) R_q \left[ \frac{\partial}{\partial \alpha_j} R_3(\alpha_1,\alpha_2,\alpha_3) \right] \mathbf{B}^{VFM}_i(\mathbf{c}) \\
   \frac{\partial f_i}{\partial \alpha_{jk}} &= 0 \\
   \frac{\partial \boldsymbol{\delta}_i}{\partial \alpha_{jk}} &= 0

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
   \frac{\partial f_i}{\partial (\cdot)} &= \frac{1}{F_i(\mathbf{c})} \mathbf{B}^{VFM}_i(\mathbf{c}) \cdot \frac{\partial}{\partial (\cdot)} \mathbf{B}^{VFM}_i(\mathbf{c}) \\
   \frac{\partial \boldsymbol{\delta}_i}{\partial ( \cdot )} &= 0

where :math:`(\cdot)` refers to one of :math:`s_{jk},o_{jk},u_{jk}`.

First derivatives
-----------------

The following table summarizes the first derivatives of the residuals needed for the Jacobian.

============== ========================== =========================== ========================
Derivative     Vector residual |epsiloni| Scalar residual :math:`f_i` Vector residual |deltai|
============== ========================== =========================== ========================
|partialg|     |depsdg|                   |dfdg|                      |ddeltadg|
|partialeuler| |depsdeuler|               0                           0
|partialc|     |depsdc|                   |dfdc|                      0
|partialk|     |depsdk|                   |dfdk|
============== ========================== =========================== ========================

where we define

.. math:: \mathbf{b}^{model}(\mathbf{r}, t; \mathbf{g}, \mathbf{k}) = \frac{\mathbf{B}^{model}(\mathbf{r}, t; \mathbf{g}, \mathbf{k})}{|| \mathbf{B}^{model}(\mathbf{r}, t; \mathbf{g}, \mathbf{k}) || }

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

.. math:: \xi_i = \frac{\partial^2 f_i}{\partial g_{nm} \partial g_{n'm'}} = \frac{1}{|| \mathbf{B}^{model}(\mathbf{r}_i, t_i; \mathbf{g},\mathbf{k})||} \left[ (\mathbf{b}^{model} \cdot \mathbf{B}_n^m(\mathbf{r}_i)) (\mathbf{b}^{model} \cdot \mathbf{B}_{n'}^{m'}(\mathbf{r}_i)) + \mathbf{B}_n^m(\mathbf{r}_i) \cdot \mathbf{B}_{n'}^{m'}(\mathbf{r}_i) \right]

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

The Jacobian matrix can be separated as

.. math::

   \mathbf{J}(\mathbf{x}) =
   \left(
   \begin{array}{cc}
     \mathbf{J}_1(\mathbf{x}) & \mathbf{J}_2(\mathbf{x})
   \end{array}
   \right)

where :math:`\mathbf{J}_1` is the Jacobian relating to the internal Gauss coefficients
:math:`g_{n,k}^m` and :math:`\mathbf{J}_2` relates to the Euler angles, fluxgate
calibration, and external field parameters. Specifically,

.. math::

   \mathbf{J}_1 &=
     \frac{\partial}{\partial g_{n,k}^m}
     \begin{pmatrix}
       \boldsymbol{\epsilon}_i \\
       f_i \\
       \boldsymbol{\delta}_i
     \end{pmatrix} =
     \begin{pmatrix}
       \mathbf{J}_{MF}^{vec} \\
       \mathbf{J}_{MF}^{scal}(\mathbf{x}) \\
       \dot{\mathbf{J}}_{MF}^{vec}
     \end{pmatrix} \\
   \mathbf{J}_2 &=
     \begin{pmatrix}
       \frac{\partial}{\partial \alpha_{jk}}
       \begin{bmatrix}
         \boldsymbol{\epsilon}_i \\
         f_i \\
         \boldsymbol{\delta}_i
       \end{bmatrix} &
       \frac{\partial}{\partial s_{jk}}
       \begin{bmatrix}
         \boldsymbol{\epsilon}_i \\
         f_i \\
         \boldsymbol{\delta}_i
       \end{bmatrix} &
       \frac{\partial}{\partial o_{jk}}
       \begin{bmatrix}
         \boldsymbol{\epsilon}_i \\
         f_i \\
         \boldsymbol{\delta}_i
       \end{bmatrix} &
       \frac{\partial}{\partial u_{jk}}
       \begin{bmatrix}
         \boldsymbol{\epsilon}_i \\
         f_i \\
         \boldsymbol{\delta}_i
       \end{bmatrix}
     \end{pmatrix} =
     \begin{pmatrix}
       \mathbf{J}_{Euler}^{vec}(\mathbf{x}) & \mathbf{J}_{Fluxgate}^{vec}(\mathbf{x}) \\
       0 & \mathbf{J}_{Fluxgate}^{scal}(\mathbf{x}) \\
       0 & 0
     \end{pmatrix}

The matrix :math:`\mathbf{J}_2` is quite sparse, and can be stored efficiently in CSR format
for example. During nonlinear least squares iterations, we must construct the normal
equations matrix

.. math::

   \mathbf{J}^T \mathbf{J} = \begin{pmatrix}
                               \mathbf{J}_1^T \mathbf{J}_1 & X \\
                               \mathbf{J}_2^T \mathbf{J}_1 & \mathbf{J}_2^T \mathbf{J}_2
                             \end{pmatrix}

The :math:`X` entry above indicates the matrix is symmetric and this portion does not need
to be computed. The upper left block is

.. math:: \mathbf{J}_1^T \mathbf{J}_1 = \mathbf{J}_{MF}^{T,vec} \mathbf{J}_{MF}^{vec} + \mathbf{J}_{MF}^{T,scal}(\mathbf{x}) \mathbf{J}_{MF}^{scal}(\mathbf{x}) + \dot{\mathbf{J}}_{MF}^{T,vec} \dot{\mathbf{J}}_{MF}^{vec}

Note that the term :math:`\mathbf{J}_{MF}^{T,vec} \mathbf{J}_{MF}^{vec} + \dot{\mathbf{J}}_{MF}^{T,vec} \dot{\mathbf{J}}_{MF}^{vec}`
does not depend on the model parameters :math:`\mathbf{x}` and can be precomputed.

Jacobian-vector products
________________________

For the normal equations method, we also
need to calculate matrix-vector products of the form :math:`\mathbf{J}^T u`. With the above structure, we have

.. math:: \mathbf{J}^T u = \begin{pmatrix}
                             \mathbf{J}_1^T \\
                             \mathbf{J}_2^T
                           \end{pmatrix} u =
                           \begin{pmatrix}
                             \mathbf{J}_1^T u \\
                             \mathbf{J}_2^T u
                           \end{pmatrix}

The matrix :math:`\mathbf{J}_2` is stored fully in CSR format, so :math:`\mathbf{J}_2^T u` can easily be
computed with a sparse BLAS operation. The matrix :math:`\mathbf{J}_1` is constructed row by row,
and so :math:`\mathbf{J}_1^T u` must be updated incrementally. Let row :math:`i` of :math:`\mathbf{J}_1`
be :math:`\mathbf{a}_i^T`. Then,

.. math:: \mathbf{J}_1^T u &= \begin{pmatrix}
                                \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n
                              \end{pmatrix} u \\
                           &= \sum_i u_i \mathbf{a}_i

With the B-spline representation of the Gauss coeffficients, the vectors :math:`\mathbf{a}_i` can
themselves be sparse, with :math:`k \times nnm` non-zero elements, where :math:`k` is the order
of the B-spline. This can be accounted for when updating the sum above.

Precomputation
______________

The vector residuals :math:`\boldsymbol{\epsilon}_i,\boldsymbol{\delta}_i` are linear
in the Gauss coefficients :math:`g_{n,k}^m` and so the portion of :math:`J^T J`
corresponding to these can be precomputed. The matrix :math:`\mathbf{J}_1`
corresponds to the core and crustal field parts of the model, which can be
further separated into vector and scalar parts:

.. math:: \mathbf{J}_1 =
            \begin{pmatrix}
              \mathbf{J}_{core} & \mathbf{J}_{crust}
            \end{pmatrix} =
            \begin{pmatrix}
              \mathbf{J}_{core}^{vec} & \mathbf{J}_{crust}^{vec} \\
              \mathbf{J}_{core}^{scal}(\mathbf{x}) & \mathbf{J}_{crust}^{scal}(\mathbf{x})
            \end{pmatrix}

The matrix :math:`\mathbf{J}_{core}` is somewhat sparse due to its B-spline parameterization,
while :math:`\mathbf{J}_{crust}` is dense. We have

.. math:: \mathbf{J}_1^T \mathbf{J}_1 =
            \begin{pmatrix}
              \mathbf{J}_{core}^{vec,T} \mathbf{J}_{core}^{vec} & X \\
              \mathbf{J}_{crust}^{vec,T} \mathbf{J}_{core}^{vec} & \mathbf{J}_{crust}^{vec,T} \mathbf{J}_{crust}^{vec}
            \end{pmatrix} +
            \begin{pmatrix}
              \mathbf{J}_{core}^{scal,T}(\mathbf{x}) \mathbf{J}_{core}^{scal}(\mathbf{x}) & X \\
              \mathbf{J}_{crust}^{scal,T}(\mathbf{x}) \mathbf{J}_{core}^{scal}(\mathbf{x}) & \mathbf{J}_{crust}^{scal,T}(\mathbf{x}) \mathbf{J}_{crust}^{scal}(\mathbf{x})
            \end{pmatrix}

The first term above can be precomputed, while the second must be computed during each
nonlinear least squares iteration.

Core field normal equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The term :math:`\mathbf{J}_{core}^{vec,T} \mathbf{J}_{core}^{vec}`
has a block representation as follows:

.. math:: \mathbf{J}_{core}^{vec,T} \mathbf{J}_{core}^{vec} =
            \begin{pmatrix}
              A_{11} & A_{12} & \cdots & A_{1n} \\
              \vdots & \vdots & \ddots & \vdots \\
              A_{n1} & A_{n2} & \cdots & A_{nn}
            \end{pmatrix}

where :math:`n` is the number of control points for each Gauss spline, and each
:math:`A_{ij}` is :code:`nnm_core`-by-:code:`nnm_core`. The matrix :math:`\mathbf{J}_{core}^{vec,T} \mathbf{J}_{core}^{vec}`
itself is :code:`nnm_core * ncontrol`-by-:code:`nnm_core * ncontrol`. The block entries are

.. math:: A_{ij} = \sum_{k=1}^{N_{vec}} N_i(t_k) N_j(t_k) B(\mathbf{r}_k) B^T(\mathbf{r}_k)

where :math:`B` is a :code:`nnm_core`-by-:code:`3` matrix given by

.. math:: B(\mathbf{r}) =
            \begin{pmatrix}
              X_n^m & Y_n^m & Z_n^m
            \end{pmatrix}

i.e. the Green's functions of the internal field model expansion. The sum for :math:`A_{ij}`
can be written as the sum of outer products of larger block matrices:

.. math:: A_{ij} = \sum_{k=1}^{\textrm{nblocks}} U_{ijk} U_{ijk}^T

where :math:`U_{ijk}` are :code:`nnm_core`-by-:code:`block_size` with the structure:

.. math:: U_{ijk} =
          \begin{pmatrix}
            \sqrt{N_i(t_1) N_j(t_1)} B(\mathbf{r}_1) & \cdots & \sqrt{N_i(t_b) N_j(t_b)} B(\mathbf{r}_b)
          \end{pmatrix}

Crustal field normal equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The crustal term :math:`\mathbf{J}_{crust}^{vec,T} \mathbf{J}_{crust}^{vec}` has the following
structure:

.. math:: \mathbf{J}_{crust}^{vec,T} \mathbf{J}_{crust}^{vec} = \sum_{k=1}^{N_{vec}} B_{crust}(\mathbf{r}_k) B_{crust}^T(\mathbf{r}_k)

which can be written in block form as

.. math:: \mathbf{J}_{crust}^{vec,T} \mathbf{J}_{crust}^{vec} = \sum_{k=1}^{\textrm{nblocks}} T_k T_k^T

where :math:`T_k` is :code:`nnm_crust`-by-:code:`block_size` and

.. math:: T_k = \begin{pmatrix} B_{crust}(\mathbf{r}_1) & \cdots & B_{crust}(\mathbf{r}_b) \end{pmatrix}

Mixed field normal equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cross term :math:`\mathbf{J}_{crust}^{vec,T} \mathbf{J}_{core}^{vec}` has the following block
structure:

.. math:: \mathbf{J}_{crust}^{vec,T} \mathbf{J}_{core}^{vec} =
            \begin{pmatrix}
              A_1 & A_2 & \cdots & A_{\textrm{ncontrol}}
            \end{pmatrix}

where,

.. math:: A_j = \sum_{k=1}^{N_{vec}} N_j(t_k) B_{crust}(\mathbf{r}_k) B_{core}^T(\mathbf{r}_k)

This can be written in block form as

.. math:: A_j = \sum_{k=1}^{\textrm{nblocks}} T_{jk} U_{jk}^T

where :math:`T_{jk}` is :code:`nnm_crust`-by-:code:`block_size`,
:math:`U_{jk}` is :code:`nnm_core`-by-:code:`block_size`, and

.. math::

   T_{jk} &= \begin{pmatrix} \sqrt{N_j(t_1)} B_{crust}(\mathbf{r}_1) & \cdots & \sqrt{N_j(t_b)} B_{crust}(\mathbf{r}_b) \end{pmatrix} \\
   U_{jk} &= \begin{pmatrix} \sqrt{N_j(t_1)} B_{core}(\mathbf{r}_1) & \cdots & \sqrt{N_j(t_b)} B_{core}(\mathbf{r}_b) \end{pmatrix}

Indexing
========

The parameter vector :math:`\mathbf{x}` is organized as follows:

.. math:: \mathbf{x} = \begin{pmatrix}
                         \mathbf{x}_{core} \\
                         \mathbf{x}_{crust} \\
                         \mathbf{x}_{Euler} \\
                         \mathbf{x}_{fluxcal}
                       \end{pmatrix}

The parameters :math:`\mathbf{x}_{core}` are the B-spline control points
for the Gauss coefficient splines up to degree :code:`nmax_core`. They are
organized as follows:

.. math:: \mathbf{x}_{core} = \begin{pmatrix}
                                \mathbf{g}_{n,0}^m \\
                                \mathbf{g}_{n,1}^m \\
                                \vdots \\
                                \mathbf{g}_{n,ncontrol}^m \\
                              \end{pmatrix}

where each :math:`\mathbf{g}_{n,k}^m` is a vector of Gauss coefficients
of length :code:`nnm_core`.

The parameters :math:`\mathbf{x}_{crust}` are the static Gauss coefficients
representing the crustal field, of length :code:`nnm_crust`.

The parameters :math:`\mathbf{x}_{Euler}` are the control points for the
Euler angle splines. For a given satellite, let the three Euler angle
splines be :math:`\alpha(t), \beta(t), \gamma(t)`, then

.. math::
   
     \alpha(t) &= \alpha^T N(t) \\
     \beta(t) &= \beta^T N(t) \\
     \gamma(t) &= \gamma^T N(t)

where :math:`\alpha, \beta, \gamma` are vectors of length :code:`ncontrol_euler`.
The coefficient vector is then organized as

.. math:: \mathbf{x}_{Euler} =
            \begin{pmatrix}
              \alpha_1 \\
              \beta_1 \\
              \gamma_1 \\
              \vdots \\
              \alpha_{nsat} \\
              \beta_{nsat} \\
              \gamma_{nsat}
            \end{pmatrix}

Each "satellite block" is of length :code:`3 * ncontrol_euler`.

Regularization
==============

Core Field Regularization
-------------------------

The core field is regularized by minimizing the third time derivative of
:math:`B_r` at the CMB:

.. math:: \left< \left| \frac{\partial^3 B_r}{\partial t^3} \right|^2 \right> = \frac{1}{\Delta t} \int dt \int d\Omega_c \left| \frac{\partial^3 B_r}{\partial t^3} \right|^2

The time integral is taken over the full time range of the Gauss coefficients with
:math:`\Delta t` equal to that time range. The surface integral is taken over the
core mantle boundary with radius :math:`c = 3485` km. Minimizing this quantity
is equivalent to minimizing the following:

.. math:: \mathbf{x}_{core}^T \Lambda_{core} \mathbf{x}_{core}

where :math:`\Lambda_{core} = G^{(3)} \otimes C` is a symmetric :code:`nnm * ncontrol`-by-:code:`nnm * ncontrol` matrix,
:math:`G^{(a)}` is a :code:`ncontrol`-by-:code:`ncontrol` B-spline Gram matrix, which is symmetric indefinite and banded.
It has entries

.. math:: G^{(a)}_{ij} = \frac{1}{\Delta t} \int dt \left( \frac{d^a N_i(t)}{dt^a} \right) \left( \frac{d^a N_j(t)}{dt^a} \right)

and :math:`C` is a diagonal :code:`nnm`-by-:code:`nnm` matrix with entries

.. math:: C_{nm,n'm'} = 4 \pi \left( \frac{a}{c} \right)^{2n+4} \frac{(n+1)^2}{2n + 1} \delta_{mm'} \delta_{nn'}

Euler Angle Regularization
--------------------------

The Euler angle splines are regularized by minimizing the second derivative (curvature) integrated
over the time interval. For each satellite, we minimize

.. math:: \textrm{curvature} = \frac{1}{\Delta t} \int dt \left( |\alpha''(t)|^2 + |\beta''(t)|^2 + |\gamma''(t)|^2 \right)

Since :math:`\alpha(t) = \alpha^T N(t), \beta(t) = \beta^T N(t), \gamma(t) = \gamma^T N(t)`, this is equivalent
to minimizing

.. math:: \mathbf{x}_{Euler}^T \Lambda_{Euler} \mathbf{x}_{Euler}

where :math:`\Lambda_{Euler}` is a block matrix with blocks for each satellite of the form :math:`G^{(2)} \otimes I_3`,
and :math:`G^{(2)}` is the Gram matrix of B-spline second derivatives using the knot vector and time period of
that particular satellite, size :code:`ncontrol_euler-by-ncontrol_euler`.

Implementation
--------------

The full minimization problem looks like

.. math:: \min_{\mathbf{x}} || \mathbf{f}(\mathbf{x}) ||^2 + \mathbf{x}^T \Lambda \mathbf{x}

where :math:`\mathbf{f}(\mathbf{x})` is a vector of all residuals and :math:`\Lambda` is the full
regularization matrix, which is block diagonal. This minimization problem is equivalent to

.. math:: \min_{\mathbf{x}} || \tilde{\mathbf{f}}(\mathbf{x}) ||^2

where the augmented residual vector is

.. math:: \tilde{\mathbf{f}}(\mathbf{x}) =
            \begin{pmatrix}
              \mathbf{f}(\mathbf{x}) \\
              L^T \mathbf{x}
            \end{pmatrix}

and :math:`L` is the Cholesky factor of :math:`\Lambda = L L^T`. The augmented Jacobian matrix is

.. math:: \tilde{J}(\mathbf{x}) =
            \begin{pmatrix}
              J(\mathbf{x}) \\
              L^T
            \end{pmatrix}
