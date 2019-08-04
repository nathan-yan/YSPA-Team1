import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from photutils import aperture_photometry, CircularAperture, CircularAnnulus, DAOStarFinder

class Image(fits.hdu.image.PrimaryHDU):
    """
    Image [fits.hdu.image.PrimaryHDU]
    An image primitive that combines multiple useful features, including wcs world/pix conversions into a single class.
    """

    def __init__(self, f, do_not_scale_data = False,
                          ignore_blank = False,
                          uint = True,
                          scale_back = None,
                          y_flip = True):
        with fits.open(f) as HDU:
            hdu = HDU[0]
            data = hdu.data.astype(np.float64)

            if y_flip:
                data = np.flipud(data)

            super(Image, self).__init__(data, hdu.header, do_not_scale_data, ignore_blank, uint, scale_back)

        # Create a wcs object
        self.wcs_obj = WCS(self.header)

        self.mean, self.median, self.stddev = None, None, None
        self.stars = None
        self.height, self.width = self.data.shape

    def get_world(self, *args):
        if len(args) == 1:
            return self.wcs_obj.all_pix2world(args[0], 1)

        return self.wcs_obj.all_pix2world(*args)

    def get_pix(self, *args):
        print(args)
        if len(args) == 1:
            pix = self.wcs_obj.all_world2pix(args[0], 1)
        else:
            pix = self.wcs_obj.all_world2pix(*args)

        pix[:, 0] %= self.height
        pix[:, 1] %= self.width

        return pix

    def set_statistics(self, statistics):
        self._stats_init = True
        self.mean, self.median, self.stddev = statistics

    def find_stars(self, stddev, fwhm):
        # Create a star finder
        if not self._stats_init:
            print("Need to initialize statistics")

        finder = DAOStarFinder(stddev * self.stddev, fwhm)

        self.stars = finder(self.data - self.median)

        return len(self.stars)

    def match_stars(self, coordinates, dist, verbose = False):
        xy = self.get_pix(coordinates)
        matched_stars = Table(names = self.stars.colnames + ['idx'])

        for j, pos in enumerate(xy):
            # Check all of the stars and evaluate their distances
            for i, star in enumerate(self.stars):
                distance = (pos[0] - star['xcentroid']) ** 2 +\
                           ((self.height - pos[1]) - star['ycentroid']) ** 2

                if distance <= dist ** 2:
                    if verbose:
                        print("Matched star #%d" % j)
                    matched_stars.add_row(list(self.stars[i]) + [j])

        return matched_stars

    def get_fluxes(self, coordinates, dist):
        """
        get_fluxes(self, coordinates: [list [ra, dec]])
        Each coordinate is an RA/Dec pair. If the image is plate solved, we can resolve the x and y pixel values for each coordinate, then find a star that matches it within a certain distance d.
        """

        xy = self.get_pix(coordinates)
        print(xy)

        #for pos in xy:
            #plt.plot([pos[0]], [pos[1]], 'o', ms = 30)

        #for star in self.stars:
    #        plt.plot([star['xcentroid']], [star['ycentroid']], 'o')


        for pos in xy:
            # Check all of the stars and evaluate their distances
            for star in self.stars:
                distance = (pos[0] - star['xcentroid']) ** 2 +\
                           ((self.height - pos[1]) - star['ycentroid']) ** 2

                if distance <= dist ** 2:
                    print("matched")
                    plt.plot([pos[0]], [1024 - pos[1]], 'o', ms = 10)


if __name__ == "__main__":
    import csv
    from astropy.coordinates import Angle

    with open("m36.csv") as cal_stars:
        reader = csv.reader(cal_stars)
        cal_star_data = []
        for row in reader:
            cal_star_data.append(row)

    cal_star_data = cal_star_data[1:]

    # Get the ra and dec of each star
    coordinates = []
    for s in range (len(cal_star_data)):
        cal_star_data[s][1] = Angle(cal_star_data[s][1] + " hours")
        cal_star_data[s][2] = Angle(cal_star_data[s][2] + " degrees")
        coordinates.append([cal_star_data[s][1].degree, cal_star_data[s][2].degree])

    a = Image("../SampleCMD/m36-R.new", y_flip = True)

    #coordinates = [[Angle("05 36 23.822 hours").degree, Angle("34 05 56.908 degrees").degree], []]

    stats = sigma_clipped_stats(a.data, sigma = 3.0)
    a.set_statistics(stats)
    a.find_stars(10., 3.)

    print(a.stars)

    import matplotlib.pyplot as plt
    print(coordinates)
    a.get_fluxes(coordinates, dist = 4)

    plt.imshow(a.data, vmin = 1000, vmax = 2600)
    plt.ylim(reversed(plt.ylim()))
    plt.show()
