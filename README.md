# K2_PSFphotometry
Perform psf photometry on K2 superstamp images

Using WCS calibrated K2 superstamp images, you can perform PSF photometry using a pre-generated effective psf map on a series of input stars. 

Example data is provided here for the Lagoon Nebula region (K2 Campaign 9). A full set of images can be downloaded from, e.g., https://archive.stsci.edu/prepds/k2superstamp/. The input list of stars must contain an id, a right ascension, and a declination (both in decimal degrees), separated by spaces (example provided here in members.txt). Image files (.fits format) are also input in a text file (example provided here in lagoonSuperstampsFile.txt).

The Python script outputs a .dat file with light curve data as well as a .png plot of theh resulting normalized flux versus time.
