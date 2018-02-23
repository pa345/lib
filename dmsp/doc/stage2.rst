******************
Stage 2 Processing
******************

Uncalibrated Data
=================

In :numref:`fig_res1` we plot the scalar residuals, Z residuals, and
H residuals of DMSP F-15 against CHAOS. The Z residuals are the difference between
the :math:`B_3` VFM axis measurement with the field model projected onto
the geodetic vertical direction. The H residuals are differences of the components
normal to the geodetic vertical direction.

.. _fig_res1:

.. figure:: /images/res1.png

   DMSP F-15 scalar, B_z, and H residuals against CHAOS for uncalibrated DMSP data

In :numref:`fig_res2` we plot the uncalibrated data against a
main field model as a function of QD latitude, separated by year. The left column
shows residuals against a field model radial component, and the right shows the residuals
against the geodetic normal component. In the radial residuals, we see larger residuals
at mid-latitudes compared to the geodetic residuals, indicating that the VFM :math:`B_3`
axis is indeed kept aligned with the geodetic vertical direction.

.. _fig_res2:

.. figure:: /images/res2.png
   :scale: 40%

   Left: Uncalibrated DMSP F-15 B_z residuals against CHAOS radial field. Right: Residuals against
   CHAOS projected into geodetic normal direction. Ascending tracks are shown in blue, descending
   tracks in green.

For 2012-2013 there is an obvious alignment problem in the :math:`B_3` component.
:numref:`fig_res3` shows the F-15 scalar residuals vs QD latitude and separated by year.

.. _fig_res3:

.. figure:: /images/res3.png
   :scale: 20%

   Uncalibrated DMSP F-15 scalar residuals against CHAOS, separated by year. Ascending tracks are shown
   in blue, descending tracks in green.

Finally, :numref:`fig_res4` shows H residuals plotted vs QD latitude and separated by year. Similar
to :numref:`fig_res2`, we plot residuals using the NEC H component on the left, and residuals against
the CHAOS vector projected onto the geodetic tangent direction on the right. We see flatter residuals
in the geodetic system, which further assures us that this is the correct system to use for the
spacecraft-fixed frame.

.. _fig_res4:

.. figure:: /images/res4.png
   :scale: 60%

   Left: Uncalibrated DMSP F-15 H residuals against CHAOS horizontal field. Right: Residuals against
   CHAOS projected onto geodetic tangent direction. Ascending tracks are shown in blue, descending
   tracks in green.

Spike removal
=============

.. _fig_spikes:

.. figure:: /images/spikes.png
   :scale: 60%

   Left: Residuals in VFM frame after removing field models (blue) with Hampel filtered signal (green).
   Right: Difference between VFM components before and after correcting for spikes. The figures
   show clearly visible single-sample spikes in all components.

In DMSP F-15 there are numerous single-sample "data spikes", especially in the VFM Z component,
but also in the other components. These spikes have typical amplitudes of a few tens of nT,
but can also be several thousand nT. They appear frequently in F-15 data, and much less often
in the other satellites (F-16, F-17, F-18), though they are found occasionally in their datasets
also.

In order to visualize the data spikes, we first remove a core field model from the VFM frame
components. The core field model (CHAOS) is first rotated into the VFM frame using the attitude
quaternions derived in the Stage 0 processing. Because the Z component is known more accurately
than X and Y, the Z residuals are much smaller than X and Y and it is easier to see the data spikes
on these curves. :numref:`fig_spikes` (left) shows the result of removing CHAOS from the VFM
components for F-15 during one day in 2008 in blue. Single sample spikes are clearly visible
in the Z plot (3rd row left). The spikes in X and Y are less visible due to the larger amplitude
of these signals (due to the less accurate attitude knowledge for these components).

One standard technique of removing single sample outliers is to use a 3-point median filter, which
is highly resistant to single-sample outliers. However, replacing the VFM data with its median
filtered version, would change a large number of data points. In the case of X and Y, since the
amplitude of the residuals is so large due to inaccurate attitude knowledge, median values could
potentially differ from the original value by tens or even hundreds of nT. I prefer a more cautious
approach, to keep as much of the original signal intact as possible, only modifying those samples
which are truly outliers.

For this purpose the Hampel filter seems to be a good choice. It works by computing the MAD
of a window containing the sample in question, and comparing the sample's deviation to the median
to the window's MAD. Samples with a large deviation are considered outliers. Currently I use
an 11 sample window (5 samples on each side plus the sample in question), and require a sample's
deviation to the median to be more than 10 times the MAD scale estimate before calling it an outlier.
Furthermore, the deviation (spike) must have a minimum amplitude of 2 nT (to avoid flagging sections of noisy oscillating
data as outliers). This procedure allows keeping nearly all of the original VFM signal intact,
and flagging only those points which are true outliers. The result of this filter is shown in green
on the left part of :numref:`fig_spikes`. This filtered curve is then used to modify the original
VFM dataset (i.e. the core field model is simply added back to the Hampel-filtered time series
to recover the corrected VFM time series).

The difference between the original and corrected VFM time series is shown in red in
:numref:`fig_spikes` (right). Most of the data points are zero, meaning the original VFM
data was left unchanged. Only a few data points are non-zero, indicating these were flagged
as spikes and corrected with the window median value by the Hampel filter. The Z component
shows the most spikes detected, most of which have similar amplitudes of about 20 nT. The
F component is also shown (bottom row), though this component is not directly filtered - it
is constructed from the VFM vector components.

Jump Correction
===============

There are regular jumps in all three VFM components for each orbit. Some of these are
due to other instruments, such as the magnetotorquers. Other jumps appear to occur when
the field changes by a certain amount, possibly due to an inadequate number of bits
in the A2D converter etc.

.. _fig_jumps:

.. figure:: /images/jumps.png
   :scale: 60%

   Top: VFM Z component prior to jump correction. Bottom: VFM Z component
   after jump correction.

:numref:`fig_jumps` (top) shows the VFM Z residual data prior after spike removal as a
function of QD latitude. We see regular jumps along the orbit, due to the range issue
discussed above. To correct for this, we again median-filter the data, and search for
differences between adjacent samples larger than 4 nT. If such a difference is found,
the difference is added to a running counter which is then added to all subsequent samples
to correct for the detected jump. This is performed on a track-by-track basis, so after
a full orbit, the counter is reset to zero to process the next track. The results
are shown in :numref:`fig_jumps` (bottom).

This procedure is done only for the VFM Z component, which is less affected by high-latitude
perturbations. The jumps are reliably detected and corrected for a large number of orbits.
For the X and Y components, we have found difficulty in identifying jumps at high-latitudes
due to the FAC/PEJ signals. Therefore for these components we set the threshold higher (1000 nT)
so very few jumps are actually corrected in these components for the time being. This is an
area which could use improvement in future work.

Initial Calibration
===================

.. _fig_res5:

.. figure:: /images/res5.png
   :scale: 60%

   DMSP F-15 scalar and B_z residuals against CHAOS after initial scalar calibration.

After the spike and jump detection/removal, we perform a scalar calibration
of the data. Currently, the set of 9 calibration parameters is computed yearly.
The resulting time series of residuals (scalar, VFM Z, and horizontal H) are shown
in :numref:`fig_res5`.

Quaternion correction
=====================

The final step in the processing is to correct the quaternions used to rotate
from the spacecraft-fixed frame into a local NEC frame. We have found that
the DMSP F-15 satellite tends to rotate slowly about the geodetic vertical
throughout its orbit, leading to misalignments in the X and Y components during
main field modeling which could not be resolved with fixed Euler angle rotations.

To correct this, at each measurement point, we calculate a rotation angle about
the VFM Z axis to bring the horizontal measurement in agreement with the CHAOS
horizontal component. Therefore at each point along orbit, we solve the following
minimization problem:

.. math:: \min || B_{CHAOS} - R_q R_3(\gamma) B_{VFM} ||^2

with

.. math::

   R_3(\gamma) =
     \left(
       \begin{array}{ccc}
         \cos{\gamma} & -sin{\gamma} & 0 \\
         \sin{\gamma} & cos{\gamma} & 0 \\
         0 & 0 & 1
       \end{array}
     \right)

and :math:`R_q` is the original rotation matrix based on the spacecraft-fixed
:math:`\hat{s}_1,\hat{s}_2,\hat{s}_3` vectors.

:numref:`fig_quat` (blue) shows the computed angle :math:`\gamma` to align the VFM
horizontal component with CHAOS over several orbits. We see that each orbit follows
the same general pattern, with a maximum deviation of around 4 degrees at the equator,
and small deviations at the poles.

.. _fig_quat:

.. figure:: /images/quat.png
   :scale: 60%

   Correction angle :math:`\gamma` versus geocentric latitude for several orbits in blue.
   In green is plotted the modeled correction angle.

Based on this figure, we attempt to model the latitude dependence of
the correction angle :math:`\gamma` with the following three parameter model,

.. math:: \gamma_M(\theta) = c_1 + c_2 \theta + c_3 \sin{\theta}

The figure shows a general :math:`\sin{\theta}` dependence, but we also allow for a constant
offset and a linear trend in the model. This model is fitted to the computed :math:`\gamma`
curve for each orbit to determine the coefficients :math:`c_1,c_2,c_3`. Then, at each point
the correction matrix :math:`R_3(\gamma_M)` is folded into the original quaternions to give
the final quaternions.
