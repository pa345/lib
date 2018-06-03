*************
Preprocessing
*************

Discard tracks with plasma bubbles
----------------------------------

If we plot N/S gradient data (Z) and additionally remove MF7, we
find a figure similar to that below:

.. figure:: /images/dZ_ns.png

   N/S gradient Z component after removing core, crust and external
   field models for 2.5 years of Swarm A data

In the equatorial region (between :math:`\pm 20` degrees QD latitude)
we see some tracks cause enhanced N/S gradients of up to about 8 nT.
The point corresponding to 8 nT near the equator in the figure above
corresponds to the following track:

.. figure:: /images/20160303_track.png

   Track recorded by Swarm A on 3 March 2016. Blue shows Z component
   after removing core, crust, and external field model. Green shows
   along-track 30s differences.

The high-frequency signal near the equator results in an even larger
N/S gradient signal. This high-frequency feature is likely due to a plasma
bubble since the data is selected for night-time and the feature is
located near the equator.
