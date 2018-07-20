***************
Electrodynamics
***************

Introduction
============

The basic equations are

.. math::

   \nabla \times \mathbf{E} &= 0 \\
   \mathbf{J} &= \sigma_0 \mathbf{E}_{\parallel} + \sigma_p \left( \mathbf{E}_{\perp} + \mathbf{u} \times \mathbf{B}_0 \right) + \sigma_h \mathbf{b} \times \left( \mathbf{E}_{\perp} + \mathbf{u} \times \mathbf{B}_0 \right) \\
   \nabla \cdot \mathbf{J} &= 0 \\
   \nabla \times \mathbf{B} &= \mu_0 \mathbf{J}

The Ohm's law equation can be simplified to

.. math:: \mathbf{J} = \sigma \left( \mathbf{E} + \mathbf{u} \times \mathbf{B}_0 \right)

where the conductivity tensor :math:`\sigma` is given by

.. math:: \sigma = \sigma_P I + \left( \sigma_0 - \sigma_p \right) \mathbf{b} \mathbf{b}^T + \sigma_h \left[ \mathbf{b} \right]_{\times}

where

.. math::

    \left[ \mathbf{b} \right]_{\times} =
    \begin{pmatrix}
      0 & -b_{\phi} & b_{\theta} \\
      b_{\phi} & 0 & -b_r \\
      -b_{\theta} & b_r & 0
    \end{pmatrix}

In a spherical coordinate basis :math:`(r,\theta,\phi)`, the conductivity tensor is given by

.. math::

   \sigma &=
   \begin{pmatrix}
     \sigma_{rr} & \sigma_{r \theta} & \sigma_{r \phi} \\
     \sigma_{\theta r} & \sigma_{\theta \theta} & \sigma_{\theta \phi} \\
     \sigma_{\phi r} & \sigma_{\phi \theta} & \sigma_{\phi \phi}
   \end{pmatrix} \\
   &=
   \begin{pmatrix}
     \sigma_p + (\sigma_0 - \sigma_p) b_r^2 & (\sigma_0 - \sigma_p) b_r b_{\theta} - \sigma_h b_{\phi} & (\sigma_0 - \sigma_p) b_r b_{\phi} + \sigma_h b_{\theta} \\
     (\sigma_0 - \sigma_p) b_{\theta} b_r + \sigma_h b_{\phi} & \sigma_p + (\sigma_0 - \sigma_p) b_{\theta}^2 & (\sigma_0 - \sigma_p) b_{\theta} b_{\phi} - \sigma_h b_r \\
     (\sigma_0 - \sigma_p) b_{\phi} b_r - \sigma_h b_{\theta} & (\sigma_0 - \sigma_p) b_{\phi} b_{\theta} + \sigma_h b_r & \sigma_p + (\sigma_0 - \sigma_p) b_{\phi}^2
   \end{pmatrix}

Sugiura and Poros 1969 use a simplified magnetic field model in which :math:`b_{\phi} = 0`. With this assumption,
:math:`b_r = \sin{I}, b_{\theta} = \cos{I}` where :math:`I` is the inclination angle. The conductivity components simplify to

.. math::

   \sigma_{rr} &= \sigma_0 \sin^2{I} + \sigma_p \cos^2{I} \\
   \sigma_{r \theta} &= \sigma_{\theta r} = (\sigma_0 - \sigma_p) \sin{I} \cos{I} \\
   \sigma_{r \phi} &= -\sigma_{\phi r} = \sigma_h \cos{I} XXX \\
   \sigma_{\theta \theta} &= \sigma_0 \cos^2{I} + \sigma_p \sin^2{I} \\
   \sigma_{\theta \phi} &= -\sigma_{\phi \theta} = -\sigma_h \sin{I} XXX \\
   \sigma_{\phi \phi} &= \sigma_p

2D Simplification
=================

The salient features of the EEJ can be obtained by a simplication to two dimensions, where all longitude
gradients are assumed to vanish:

.. math:: \frac{\partial \mathbf{J}}{\partial \phi} = \frac{\partial \mathbf{E}}{\partial \phi} = 0

With this assumption, the equation :math:`\nabla \times \mathbf{E} = 0` immediately becomes

.. math::

   \partial_r(r E_{\theta}) - \partial_{\theta}(E_r) &= 0 \\
   \left.
   \begin{array}{rr}
   \partial_{\theta}(\sin{\theta} E_{\phi}) & = 0 \\
   \partial_r(r E_{\phi}) & = 0
   \end{array} \right\} &
   \Rightarrow
   E_{\phi} = \frac{R E_{{\phi}_0}}{r \sin{\theta}}

where :math:`R` is a constant of integration and can be taken as a reference radius. The equation
:math:`\nabla \cdot \mathbf{J} = 0` becomes

.. math:: \frac{1}{r^2} \partial_r \left( r^2 J_r \right) + \frac{1}{r \sin{\theta}} \partial_{\theta} \left( \sin{\theta} J_{\theta} \right) = 0

The two components :math:`J_r` and :math:`J_{\theta}` can be written in terms of a single *current stream function* :math:`\psi` as

.. math::

   J_r &= \frac{c}{r^2 \sin{\theta}} \partial_{\theta} \psi \\
   J_{\theta} &= -\frac{c}{r \sin{\theta}} \partial_r \psi

where :math:`c` is an arbitrary constant. Let :math:`\mathbf{W} = \sigma \left( \mathbf{u} \times \mathbf{B}_0 \right)`. Then
Ohm's law :math:`\mathbf{J} = \sigma \left( \mathbf{E} + \mathbf{u} \times \mathbf{B}_0 \right)` becomes

.. math::

   J_r &= \sigma_{rr} E_r + \sigma_{r\theta} E_{\theta} + \sigma_{r \phi} E_{\phi} + W_r \\
   J_{\theta} &= \sigma_{\theta r} E_r + \sigma_{\theta \theta} E_{\theta} + \sigma_{\theta \phi} E_{\phi} + W_{\theta} \\
   J_{\phi} &= \sigma_{\phi r} E_r + \sigma_{\phi \theta} E_{\theta} + \sigma_{\phi \phi} E_{\phi} + W_{\phi}

Taking the combinations :math:`\sigma_{\theta \theta} J_r - \sigma_{r \theta} J_{\theta}` and
:math:`\sigma_{\theta r} J_r - \sigma_{rr} J_{\theta}` yields

.. math::
   
   \sigma_{\theta \theta} J_r - \sigma_{r \theta} J_{\theta} &= \alpha E_r + \beta E_{\phi} + \sigma_{\theta \theta} W_r - \sigma_{r \theta} W_{\theta} \\
   \sigma_{\theta r} J_r - \sigma_{rr} J_{\theta} &= -\alpha E_{\theta} + \gamma E_{\phi} + \sigma_{\theta r} W_r - \sigma_{rr} W_{\theta}

with

.. math::

   \alpha = \sigma_{rr} \sigma_{\theta \theta} - \sigma_{r \theta} \sigma_{\theta r} \\
   \beta = \sigma_{\theta \theta} \sigma_{r \phi} - \sigma_{r \theta} \sigma_{\theta \phi} \\
   \gamma = \sigma_{r \phi} \sigma_{\theta r} - \sigma_{r r} \sigma_{\theta \phi}

Rearranging slightly gives

.. math::

   E_r &= \frac{1}{\alpha} \left( \sigma_{\theta \theta} J_r - \sigma_{r \theta} J_{\theta} - \beta E_{\phi} - \sigma_{\theta \theta} W_r + \sigma_{r \theta} W_{\theta} \right) \\
   E_{\theta} &= -\frac{1}{\alpha} \left( \sigma_{\theta r} J_r - \sigma_{rr} J_{\theta} - \gamma E_{\phi} - \sigma_{\theta r} W_r + \sigma_{rr} W_{\theta} \right)

Finally, the condition :math:`\partial_r(r E_{\theta}) - \partial_{\theta}(E_r) = 0` gives a 2D linear
partial differential equation for the stream function :math:`\psi`,

.. math::

   f_1 \frac{\partial^2 \psi}{\partial r^2} + 2 f_2 \frac{\partial^2 \psi}{\partial r \partial \theta} + f_3 \frac{\partial^2 \psi}{\partial \theta^2} + f_4 \frac{\partial \psi}{\partial r} + f_5 \frac{\partial \psi}{\partial \theta} = g

with

.. math::

   f_1 &= r^2 \sigma_{rr} \\
   f_2 &= \frac{1}{2} r \left( \sigma_{r \theta} + \sigma_{\theta r} \right) \\
   f_3 &= \sigma_{\theta \theta} \\
   f_4 &= r \alpha \left[ r \frac{\partial}{\partial r} \left( \frac{\sigma_{rr}}{\alpha} \right) + \sin{\theta} \frac{\partial}{\partial \theta} \left( \frac{1}{\sin{\theta}} \frac{\sigma_{r \theta}}{\alpha} \right) \right] \\
   f_5 &= \alpha \left[ r^2 \frac{\partial}{\partial r} \left( \frac{1}{r} \frac{\sigma_{\theta r}}{\alpha} \right) + \sin{\theta} \frac{\partial}{\partial \theta} \left( \frac{1}{\sin{\theta}} \frac{\sigma_{\theta \theta}}{\alpha} \right) \right] \\
   g &= \frac{R}{c} r \alpha \left[ r \frac{\partial}{\partial r} \left( \frac{\gamma}{\alpha} \right) + \sin{\theta} \frac{\partial}{\partial \theta} \left( \frac{1}{\sin{\theta}} \frac{\beta}{\alpha} \right) \right] E_{\phi_0} + \\
     & \qquad \frac{r^2}{c} \alpha \sin{\theta} \left\{ \frac{\partial}{\partial r} \left[ \frac{r}{\alpha} \left( \sigma_{\theta r} W_r - \sigma_{rr} W_{\theta} \right) \right] + \frac{\partial}{\partial \theta} \left[ \frac{1}{\alpha} \left( \sigma_{\theta \theta} W_r - \sigma_{r \theta} W_{\theta} \right) \right] \right\}

Based on the expression for :math:`g`, we choose the constant :math:`c` to be

.. math:: c = R^2

which gives

.. math::

   g &= \frac{r}{R} \alpha \left[ r \frac{\partial}{\partial r} \left( \frac{\gamma}{\alpha} \right) + \sin{\theta} \frac{\partial}{\partial \theta} \left( \frac{1}{\sin{\theta}} \frac{\beta}{\alpha} \right) \right] E_{\phi_0} + \\
     & \qquad \left( \frac{r}{R} \right)^2 \alpha \sin{\theta} \left\{ \frac{\partial}{\partial r} \left[ \frac{r}{\alpha} \left( \sigma_{\theta r} W_r - \sigma_{rr} W_{\theta} \right) \right] + \frac{\partial}{\partial \theta} \left[ \frac{1}{\alpha} \left( \sigma_{\theta \theta} W_r - \sigma_{r \theta} W_{\theta} \right) \right] \right\}

and also

.. math::

   J_r &= \left( \frac{R}{r} \right)^2 \frac{1}{\sin{\theta}} \partial_{\theta} \psi \\
   J_{\theta} &= - \left( \frac{R}{r} \right) \frac{R}{\sin{\theta}} \partial_r \psi

This choice of :math:`c` means that :math:`\psi` has units of current density :math:`A/m^2`.

Nondimensionalization
=====================

To solve the 2D PDE numerically, we first remove the dimensions. Define

.. math::

   \tilde{\sigma}_{xy} &= \frac{\sigma_{xy}}{\sigma_s} \\
   \tilde{r} &= \frac{r}{r_s} \\
   \tilde{E}_{\phi_0} &= \frac{E_{\phi_0}}{E_s} \\
   \tilde{\psi} &= \frac{\psi}{\psi_s}

where :math:`\sigma_s, r_s, E_s, \psi_s` are constants carrying physical units and the
:math:`\tilde{x}` quantities are dimensionless. With these definitions, we find

.. math::

   f_1 &= r_s^2 \sigma_s \tilde{f}_1 \\
   f_2 &= r_s \sigma_s \tilde{f}_2 \\
   f_3 &= \sigma_s \tilde{f}_3 \\
   f_4 &= r_s \sigma_s \tilde{f}_4 \\
   f_5 &= \sigma_s \tilde{f}_5 \\
   g &= \sigma_s^2 E_s \tilde{g}

The nondimensionalized equation

.. math::

   \tilde{f}_1 \frac{\partial^2 \tilde{\psi}}{\partial \tilde{r}^2} + 2 \tilde{f}_2 \frac{\partial^2 \tilde{\psi}}{\partial \tilde{r} \partial \theta} + \tilde{f}_3 \frac{\partial^2 \tilde{\psi}}{\partial \theta^2} + \tilde{f}_4 \frac{\partial \tilde{\psi}}{\partial \tilde{r}} + \tilde{f}_5 \frac{\partial \tilde{\psi}}{\partial \theta} = \tilde{g}

will be satisfied, provided the following relationship holds:

.. math:: \frac{\sigma_s E_s}{\psi_s} = 1

This means, we have the freedom to choose :math:`r_s` and any two of :math:`\sigma_s,E_s,\psi_s` as we want.
Then the fourth scale is automatically determined by the equation above. Since
:math:`\mathbf{W} = \sigma \mathbf{u} \times \mathbf{B}_0` appears in the right hand side
:math:`g`, we can also define

.. math::

   \tilde{\mathbf{u}} &= \frac{\mathbf{u}}{u_s} \\
   \tilde{\mathbf{B}_0} &= \frac{\mathbf{B}_0}{B_s}

These factors are constrained by the relation

.. math:: u_s B_s = E_s

so we have the freedom to choose one of :math:`u_s,B_s`.
