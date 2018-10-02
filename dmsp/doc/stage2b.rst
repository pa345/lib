************************
Stage 2 Processing (NEW)
************************

Attitude correction
===================

The VFM X and Y components deviate significantly from the CHAOS model
rotated into the VFM frame using the quaternions computed in Stage 0.
I believe this is due to the spacecraft slowly rotating about the vertical
direction during its orbit. Before features like jumps and spikes can be
corrected in the VFM X and Y components, we need to reduce the residuals
with CHAOS in order to see these features more clearly.

.. _fig_attitude1:

.. figure:: /images/attitude1.png
   :scale: 60%

   DMSP F-16 VFM residuals before (left) and after (right) attitude correction applied.

:numref:`fig_attitude1` (left) shows the large VFM X and Y residuals for a few orbits
of DSMP F-16. The VFM Z residuals do not appear to suffer greatly from this issue.

To correct this, at each measurement point, we calculate a rotation angle about
the VFM Z axis to bring the horizontal measurement in agreement with the CHAOS
horizontal component. Therefore at each point along orbit, we solve the following
minimization problem:

.. math:: \min_{\gamma} || B_{CHAOS} - R_q R_3(\gamma) B_{VFM} ||^2

with

.. math::

   R_3(\gamma) =
     \left(
       \begin{array}{ccc}
         \cos{\gamma} & -\sin{\gamma} & 0 \\
         \sin{\gamma} & \cos{\gamma} & 0 \\
         0 & 0 & 1
       \end{array}
     \right)

and :math:`R_q` is the original rotation matrix based on the spacecraft-fixed
:math:`\hat{s}_1,\hat{s}_2,\hat{s}_3` vectors. :numref:`fig_attitude2` (blue)
shows the resulting angle :math:`\gamma` versus geocentric latitude for a few
orbits of DMSP F-16. Due to the smooth variation seen in :math:`\gamma` over these
orbits, we define a 3 parameter model,

.. math:: \gamma(\theta) = c_1 + \sin{(c_2 + c_3 \theta)}

The parameters :math:`c_1,c_2,c_3` are determined for each orbital track via
a least squares minimization. The resulting modeled angle is shown in
:numref:`fig_attitude2` (green).

.. _fig_attitude2:

.. figure:: /images/attitude2.png
   :scale: 60%

   Correction angle :math:`\gamma` versus geocentric latitude for several orbits in blue.
   In green is plotted the modeled correction angle.

Finally, for each measurement point, we compute a new quaternion attitude transformation,
defined by

.. math:: R_q' = R_q R_3(\gamma)

where :math:`R_q` is the original attitude based on the :math:`\hat{s}_1,\hat{s}_2,\hat{s}_3`
vectors. This new attitude matrix is used in subsequent calculations.

Jump Detection and Correction
=============================

There exist regular jumps in all three VFM components for each orbit. Some of these
are due to spacecraft fields (such as magnetotorquers).

.. _fig_jumps:

.. figure:: /images/jumps.png
   :scale: 60%

   VFM X, Y and Z components prior to jump correction (blue). Test statistics for
   jump detection shown in green. Jumps detected shown in light red. Dark red shows
   the corrected signal.

Spike Correction
================

The jump correction algorithm sometimes introduces single-sample spikes, or impulses,
in the data. Such impulses also exist in the original dataset. These are detected and
corrected with an impulse rejection filter.

.. _fig_spikes:

.. figure:: /images/spikes.png
   :scale: 60%

   Spike removal procedure for VFM X, Y and Z components. Blue squares indicate detected
   impulses, which are replaced by the window median.

Along-track RMS Filtering
=========================

After the jumps and spikes are corrected, the root-mean-square (rms) of the scalar residual
is calculated along each half orbit. Orbits with an rms larger than some threshold are
flagged and are not used in the scalar calibration step. The current rms threshold is
:code:`50 nT`. :numref:`fig_rms` shows an example scalar rms plot versus longitude for 1 year
of F-17 data (2015).

.. _fig_rms:

.. figure:: /images/rms_lon.png
   :scale: 60%

Scalar Calibration
==================
