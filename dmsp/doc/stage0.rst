******************
Stage 0 Processing
******************

In this step we combine DMSP ephemeris files with the MFR magnetometer
data files to produce files in the Swarm CDF format containing
positions and measurements. The scalar field measurement is set to the
modulus of the vector measurement.

We also calculate quaternions at each measurement point to rotate
from a spacecraft fixed frame into a local NEC frame. Since DMSP
does not carry a star camera, this is done by defining a set of
spacecraft axes as follows,

.. math::

   \begin{array}{ll}
     \textrm{geodetic downward} & \hat{s}_1 = -\hat{e}_{\mu} \\
     \textrm{velocity direction} & \hat{s}_2 = \hat{s}_3 \times \hat{s}_1 \\
     \textrm{orbit normal, toward sun} & \hat{s}_3 = \left( \hat{s}_1 \times \mathbf{v} \right) / \left| \hat{s}_1 \times \mathbf{v} \right|
   \end{array}

Here, :math:`\hat{e}_{\mu}` is the outward normal to the local ellipsoid as defined
by WGS84. So :math:`\hat{s}_3` is chosen as the local geodetic downward direction, since
the DMSP attitude control system is designed to keep the spacecraft fixed with
respect to the geodetic normal to within 0.01 degrees. The vector :math:`\hat{s}_2`
is defined to be normal to the velocity direction and geodetic downward. This is
approximately eastward on ascending orbits and westward on descending orbits.
The vector :math:`\hat{s}_1` then completes the right handed set, and is in
the velocity direction.
