import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import wcs
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats, gaussian_sigma_to_fwhm
from photutils.background import MADStdBackgroundRMS, MMMBackground, Background2D
from photutils.psf import extract_stars, EPSFBuilder, ImagePSF, SourceGrouper
from astropy.modeling.fitting import LevMarLSQFitter, SLSQPLSQFitter, SimplexLSQFitter, FittingWithOutlierRemoval
from astropy.visualization import simple_norm
from photutils.detection import IRAFStarFinder, DAOStarFinder
from photutils.psf import IterativePSFPhotometry
from astroquery.gaia import Gaia
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import SegmentationImage, detect_sources
from photutils.aperture import CircularAperture, aperture_photometry
from astropy.nddata import Cutout2D
from astropy import units as u
import time
import gc

def create_lightcurve(filelist, r, epsf, test_ra, test_dec, ID):
    counter = 0
    times = np.array([])
    fluxes= np.array([])
    fluxes2 = np.array([])
    cadences = np.array([])
    quality = np.array([])
    pos_corr = np.array([])
    xfits = np.array([])
    xfits2 = np.array([])
    yfits = np.array([])
    yfits2 = np.array([])
    residual_fluxes = np.array([])
    aperture_fluxes = np.array([])
    ymedian = 0.2710156237265835
    index = np.array([])
    uncertainty = np.array([])

    predicted_x = np.array([])
    predicted_y = np.array([])
    oversample = 4

    ra = r['ra']
    dec = r['dec']
    flux = r['phot_rp_mean_flux']
    mag = r['phot_rp_mean_mag']
    mag2 = r['phot_rp_mean_mag']
    ra2 = r['ra']
    dec2 = r['dec']
    neighbors = np.where(np.sqrt((dec - test_dec)**2 + ((ra - test_ra)*(np.cos(test_dec*np.pi/180.0)))**2) < 9.5/1125)
    the_star = np.where(np.sqrt((dec - test_dec)**2 + ((ra - test_ra)*(np.cos(test_dec*np.pi/180.0)))**2) < 0.5/1125)
    ra = ra[neighbors]
    mag = mag[neighbors]
    dec = dec[neighbors]
    mag2 = mag2[the_star]
    ra2 = ra2[the_star]
    dec2 = dec2[the_star]
    ra_fixed = []
    dec_fixed = []
    ra_flex = []
    dec_flex = []
    fix_center_pos = False
    for i in np.arange(len(mag)):
        if(len(mag2.data) ==0 or str(mag2.data[0]) == "--"):
            print("Bad Data")
            return
        elif(np.sqrt((dec.data[i] - test_dec)**2 + ((ra.data[i] - test_ra)*(np.cos(test_dec*np.pi/180.0)))**2) < 0.5/1125):
            continue
        elif(str(mag.data[i])== '--'):
            continue
        elif(np.sqrt((dec.data[i] - test_dec)**2 + ((ra.data[i] - test_ra)*(np.cos(test_dec*np.pi/180.0)))**2) < ((50/mag2.data[0])/1125) and mag.data[i]-2 > mag2.data[0]):
            continue
        elif(np.sqrt((dec.data[i] - test_dec)**2 + ((ra.data[i] - test_ra)*(np.cos(test_dec*np.pi/180.0)))**2) < ((50/mag2.data[0])/1125) and mag.data[i]-1 > mag2.data[0]):
            ra_fixed.append(ra.data[i])
            dec_fixed.append(dec.data[i])
        elif(mag.data[i] > 15):
            ra_fixed.append(ra.data[i])
            dec_fixed.append(dec.data[i])
        else:
            ra_flex.append(ra.data[i])
            dec_flex.append(dec.data[i])

    while counter<len(filelist):
     try:
        start = time.time()
        print(counter)
        file = filelist[counter]
        data = fits.getdata(file)
        header = fits.getheader(file)
        wcs_transform = wcs.WCS(header)
 
        xstar, ystar = wcs_transform.wcs_world2pix([[test_ra, test_dec]], 0)[0]
        if not((xstar - 8 > 0 and xstar + 8 < 110 and ystar - 8 > 10 and ystar + 8 < 134) or (xstar - 8 > 110 and xstar + 8 < 219 and ystar - 8 > 0 and ystar + 8 < 124) or (xstar > 102 and xstar < 110 and ystar <124 and ystar > 18) or (xstar > 110 and xstar < 118 and ystar > 8 and ystar < 116)):
            print("Out of Bounds")
            return
        
        #taking a 9 by 9 cutout from the image centered on the star we want to focus on
        box_d = 16
        position = (np.round(xstar), np.round(ystar))
        size = (box_d, box_d)  
        cutout = Cutout2D(data, position, size) 
        
        predicted_x = np.append(predicted_x, xstar)
        predicted_y = np.append(predicted_y, ystar)
        
        bkgrms = MADStdBackgroundRMS()
        std = bkgrms(data)

        min_sep = 1.24 * 6.4
        grouper = SourceGrouper(min_separation=min_sep)        
       
        xcoords = []
        ycoords = []
        xclose = []
        yclose = []
        
        fitter = LevMarLSQFitter()
        iraffind = IRAFStarFinder(threshold=3.0*std, fwhm=1.24, min_separation=1.24, roundness_range=(-3.0,5.0), sharpness_range=(-2.0,8.0))
        photometry = IterativePSFPhotometry(finder = iraffind, grouper = grouper, psf_model = epsf, fitter = fitter, fit_shape = (7, 7), aperture_radius = 7)
        mean, median, std = sigma_clipped_stats(data, sigma=2.0)
       
        center_x, center_y = box_d/2 - 0.5, box_d/2 - 0.5
        for i in range(0, len(ra_fixed)): 
            xstar1, ystar1 = wcs_transform.wcs_world2pix([[ra_fixed[i], dec_fixed[i]]], 0)[0]
            if(np.absolute(xstar1-xstar) > 0.1 and np.absolute(ystar1-ystar) > 0.1):
                if(np.sqrt((xstar1-xstar)**2 + (ystar1 - ystar)**2) > 0.25): 
                    xclose.append(xstar1-cutout.center_original[0]+center_x)
                    yclose.append(ystar1-cutout.center_original[1]+center_y)

        for i in range(0, len(ra_flex)):
            xstar1, ystar1 = wcs_transform.wcs_world2pix([[ra_flex[i], dec_flex[i]]], 0)[0]
            if(np.absolute(xstar1-xstar) > 0.1 and np.absolute(ystar1-ystar) > 0.1):
                if(np.sqrt((xstar1-xstar)**2 + (ystar1 - ystar)**2) > 0.25): 
                    xcoords.append(xstar1-cutout.center_original[0]+center_x)
                    ycoords.append(ystar1-cutout.center_original[1]+center_y)
            
        epsf.x_0.fixed = False
        epsf.y_0.fixed = False

        threshold = 3.0*std
        segm = detect_sources(cutout.data, threshold, 2)
        mask = segm.make_source_mask()

        use_med = True
        for i in range(0, len(cutout.data)):
            for j in range(0, len(cutout.data[0])):
                if mask[i][j] == False:
                    use_med = False

        if(use_med):
            mean, median, std = sigma_clipped_stats(data, sigma=2.0)
            cutout.data = cutout.data - median
        
        else: 
            bkg_estimator = MedianBackground()
            try: 
                for i in range(0, 5):
                    bkg_est = Background2D(data = cutout.data, box_size = 1, mask = mask, exclude_percentile = 99, bkg_estimator = bkg_estimator)
                    cutout.data = cutout.data - bkg_est.background
            except(ValueError):
                mean, median, std = sigma_clipped_stats(data, sigma=2.0)
                cutout.data = cutout.data - median

        residual_image = cutout.data

        x_pos = xstar - cutout.center_original[0] + center_x
        y_pos = ystar - cutout.center_original[1] + center_y
#        print(x_pos, y_pos)

        if(len(xcoords) > 0):
            posStar = Table(names = ['x_0', 'y_0'], data = [xcoords, ycoords])
#            print(cutout.data)
            result_tab = photometry(cutout.data, init_params=posStar)

# First remove the star from this results table:

            distances = np.sqrt((result_tab['x_fit'] - x_pos)**2 + (result_tab['y_fit'] - y_pos)**2)
            bad_star_idx = np.where(distances < 1)
#            print(bad_star_idx[0])
            if (len(bad_star_idx[0]) > 0):
                photometry.results.remove_row(bad_star_idx[0][0])

            residual_image = photometry.make_residual_image(cutout.data)

        epsf.x_0.fixed = True
        epsf.y_0.fixed = True

        
        if(len(xclose) > 0):
            posStar = Table(names = ['x_0', 'y_0'], data = [xclose, yclose])
            gc.collect()
            result_tab = photometry(residual_image, init_params=posStar)

# Again remove the star of interest if it is there:
            distances = np.sqrt((result_tab['x_fit'] - x_pos)**2 + (result_tab['y_fit'] - y_pos)**2)
            bad_star_idx = np.where(distances < 1)
#            print(bad_star_idx[0])
            if (len(bad_star_idx[0]) > 0):
                photometry.results.remove_row(bad_star_idx[0][0])

                residual_image = photometry.make_residual_image(residual_image)

        epsf.x_0.fixed = False
        epsf.y_0.fixed = False
        
        #focusing on the star of interest from the image
        xcoords = []
        ycoords = []
        #coordinates of the star of interest are just the center of the image
        xcoords.append(x_pos)
        ycoords.append(y_pos)
        posStar1 = Table(names = ['x_0', 'y_0'], data = [xcoords, ycoords])
        #doing photometry again on the residual with neighboring stars subtracted out
        gc.collect()
        result_tab_1 = photometry(residual_image, init_params = posStar1)
        residual_image_star1 = photometry.make_residual_image(residual_image)

# Isolate the star of interest again:
        distances = np.sqrt((result_tab_1['x_fit'] - x_pos)**2 + (result_tab_1['y_fit'] - y_pos)**2)
        good_star_idx = np.where(distances < 1)

#        if(len(result_tab_1[0]) == 12):
        if((len(result_tab_1) > 0) & (len(good_star_idx[0]) > 0)):
            xfits = np.append(xfits, result_tab_1['x_fit'][good_star_idx[0][0]])
        else:
            xfits = np.append(xfits, x_pos)
#        if(len(result_tab_1[0]) == 12):
        if((len(result_tab_1) > 0) & (len(good_star_idx[0]) > 0)):
            yfits = np.append(yfits, result_tab_1['y_fit'][good_star_idx[0][0]])
        else:
            yfits = np.append(yfits, y_pos)
#        if(len(result_tab_1[0]) == 12):
        if((len(result_tab_1) > 0) & (len(good_star_idx[0]) > 0)):
            uncertainty = np.append(uncertainty, result_tab_1['flux_err'][good_star_idx[0][0]])
        elif(len(result_tab_1) < 1 and counter != 0):
            uncertainty = np.append(uncertainty, uncertainty[counter-1]*100)
        else:
            uncertainty = np.append(uncertainty, 400000/mag2.data[0])

        position = [xfits[counter], yfits[counter]]
##        aperture = CircularAperture(position, r = 3.)
##        gc.collect()
##        phot_table = aperture_photometry(cutout.data, aperture)

##        aperture_fluxes = np.append(aperture_fluxes, phot_table['aperture_sum'])
##        residual_phot_table = aperture_photometry(residual_image_star1, aperture)

##        residual_fluxes = np.append(residual_fluxes, residual_phot_table['aperture_sum'])

#        print(len(result_tab_1))
        if((len(result_tab_1) > 0) & (len(good_star_idx[0]) > 0)):
          fluxes = np.append(fluxes, result_tab_1['flux_fit'][good_star_idx[0][0]])
        else:
          fluxes = np.append(fluxes, np.nan)

        times = np.append(times, header['MIDTIME'])
        cadences = np.append(cadences, header['CADENCEN'])
        quality = np.append(quality, header['QUALITY'])
        pos_corr = np.append(pos_corr, header['POSCORR1'])
        index = np.append(index, counter)

#        print(fluxes[counter])
        
        end = time.time()
        print(end - start)
        counter = counter+1

     except:
       print("Possibly bad input file!")
       counter = counter+1 

    outputFile = open('EPIC_' + str(ID.split("_")[0]) + '_psf.dat', 'w')
        
    line = "Dates" + " " + "Cadences" + " " + "PSFflux" + " " + "Uncert" + "Xpos" + " " + "Ypos" +  " " + "Quality" + "\n"
    outputFile.write(line)
    for d, c, f, u, x, y, q in zip(times, cadences, fluxes, uncertainty, xfits, yfits, quality): 
        line = str(d-2454833) + ' ' + str(int(c)) 
        line += ' ' + str(f) + ' ' + str(u) 
        line += ' ' + str(f) + ' ' + str(u) 
        line += ' ' + str(f) + ' ' + str(u) 
        line += ' ' + str(f) + ' ' + str(u) 
        line += ' ' + str(x) + ' ' + str(y) + ' ' + str(int(q)) + "\n"
        outputFile.write(line)
    
    outputFile.close()

#    print('lens: ', len(times), len(fluxes))

    plt.clf()
    plt.plot(times, fluxes/np.nanmedian(fluxes),'bo',markersize=1.5)
    plt.xlabel('Time [BJD-2454833]')
    plt.ylabel('Norm. flux')
    plt.savefig('lightcurve_' + str(ID.split("_")[0]) + '.png')

    plt.clf()
    plt.plot(times, fluxes/np.nanmedian(fluxes),'bo',markersize=1.5)
    plt.ylim(1-2.5*np.nanstd(fluxes/np.nanmedian(fluxes)),1+2.5*np.nanstd(fluxes/np.nanmedian(fluxes)))
    plt.xlabel('Time [BJD-2454833]')
    plt.ylabel('Norm. flux')
    plt.savefig('lightcurve_' + str(ID.split("_")[0]) + '(with_limits).png')

def main():
  files = 'lagoonSuperstampsFile.txt'
  filelist = np.loadtxt(files,dtype=str)


# Pull in the pre-computed effective PSF profile:

  epsf_data = fits.getdata('lagoon_PSF_8thmodel')
  epsf = ImagePSF(epsf_data, oversampling=4, origin=None)

# Get the locations and magnitudes of stars in the Lagoon superstamp region:
  
  coord = SkyCoord(ra=271.0473081, dec=-24.3894347, unit=(u.degree, u.degree), frame='icrs')
  radius = u.Quantity(0.15, u.deg)
  Gaia.ROW_LIMIT = 10000
  r = Gaia.query_object_async(coordinate=coord, radius = radius)

# Pull in the stars you want photometry for:
  lagoonmemberlist = open('members.txt', "r")
  lines = lagoonmemberlist.readlines()
  ra_positions = []
  dec_positions = []
  ID_numbers = []
  for x in lines:
     ra_positions.append(x.split(" ")[1])
     dec_positions.append(x.split(" ")[2].split("\n")[0])
     ID_numbers.append(x.split(" ")[0])

 
  counter = 0
  for i in range(counter, len(ra_positions)):
      test_ra, test_dec = float(ra_positions[i]), float(dec_positions[i])
      print(test_ra)
      print(test_dec)
      ID_number = str(ID_numbers[i])
      create_lightcurve(filelist, r, epsf, test_ra,test_dec, ID_number)
      counter += 1



# Standard boilerplate to call the main() function.
if __name__ == '__main__':
  main()
