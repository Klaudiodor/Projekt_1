import sys, argparse

from math import sin, cos, sqrt, atan, atan2, degrees, radians, pi, tan, floor, modf

from numpy import abs, array 

o = object()

class Transformacje:
    def __init__(self, file, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0  # semimajor_axis
            self.b = 6356752.31424518  # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2)  # eccentricity  WGS84:0.0818191910428
        self.ecc2 = (2 * self.flat - self.flat ** 2)  # eccentricity**2

        self.xyz = self.readfile(file)
        self.plh = self.xyz2plh(False)

    def xyz2plh(self, data_write, output='dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny.
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim,

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """

        """
        Processes a list of coordinate tuples (X, Y, Z) and converts them to geodetic coordinates.
        """
        results = []
        for X, Y, Z in self.xyz:
            r = sqrt(X ** 2 + Y ** 2)
            lat_prev = atan(Z / (r * (1 - self.ecc2)))
            lat = 0
            while abs(lat_prev - lat) > 0.000001 / 206265:
                lat_prev = lat
                N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev) ** 2)
                h = r / cos(lat_prev) - N
                lat = atan((Z / r) * (((1 - self.ecc2 * N / (N + h)) ** (-1))))
            lon = atan(Y / X)
            N = self.a / sqrt(1 - self.ecc2 * sin(lat) ** 2)
            h = r / cos(lat) - N

            if output == 'dec_degree':
                result = (degrees(lat), degrees(lon), h)
            elif output == 'dms':
                result = (self.deg2dms(degrees(lat)), self.deg2dms(degrees(lon)), h)
            else:
                raise NotImplementedError(f"{output} - output format not defined")

            results.append(result)

        if data_write:
            self.write2file("xyz2plh", results)

        return results
