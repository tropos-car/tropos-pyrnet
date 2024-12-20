.. PyrNet documentation master file, created by

Welcome to PyrNet's documentation!
==================================

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: Documentation:

   pyrnet.md
   campaigns.md
   howto.md


Calibration absolute
==================================
PyrNet stations pyranometer are absolut calibrated versus secondary standard pyranometers in field.
To verify comparability of the individual irradiance observations of each station, these stations are colectively calibrated within a time frame close to a major measurement campaign.

.. toctree::
   :maxdepth: 1
   :numbered:
   :caption: Calibration absolute:

   calibration/calibration_s2vsr.md
   calibration/calibration_metpvnet.md
   calibration/calibration_melcol.md
   calibration/calibration_ioprao.md
   calibration/calibration_hope-melpitz.md

Calibration cosine error
==================================
PyrNet stations pyranometer are compared to clear sky reference measurements. A cosine correction function is determined.

.. toctree::
   :maxdepth: 1
   :numbered:
   :caption: Calibration cosine error:

   calibration/mean_cosine_correction.md
   calibration/cos_corr_s2vsr.ipynb
   calibration/calibration_melcol_cc.md
   calibration/calibration_ioprao_cc.md


Code Documentation
==================================

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: Code:

   code_docs/cli.md
   nbs/pyrnet.ipynb
   nbs/utils.ipynb
   nbs/logger.ipynb
   nbs/reports.ipynb
   nbs/data.ipynb
   nbs/qcrad.ipynb
   code_docs/autodocs.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
