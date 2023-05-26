# How to ... 

## ... keep full resolution of geo coordinates?
From the processing step l1a -> l1b the decision whether to keep lat and lon coordinates time resolution (and interpolate to measurement sample), or average them over the maintenance interval has to be made.

In the configuration file provided for ```pyrnet process l1b -c``` the keyword "*average_latlon*" switches the behaviour.
The default is *True*, so the coordinates are averaged over the maintenance period.