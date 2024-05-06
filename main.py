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

    def plh2xyz(self):
        """
        Converts a list of geodetic coordinates (phi, lam, h) to Cartesian coordinates (X, Y, Z).

        Parameters
        ----------
        coordinates : list of tuples
            Each tuple contains (phi, lam, h) in degrees and meters.

        Returns
        -------
        list of tuples
            Each tuple contains (X, Y, Z) as floats, rounded to 3 decimal places.
        """
        results = []
        for phi, lam, h in self.plh:
            phi = radians(phi)
            lam = radians(lam)

            Rn = self.a / sqrt(1 - self.ecc2 * sin(phi) ** 2)
            q = Rn * self.ecc2 * sin(phi)

            X = (Rn + h) * cos(phi) * cos(lam)
            Y = (Rn + h) * cos(phi) * sin(lam)
            Z = (Rn + h) * sin(phi) - q

            results.append((round(X, 3), round(Y, 3), round(Z, 3)))

            self.write2file("plh2xyz", results)
        return results

    def deg2dms(self, deg):
        """
        Zamiana stopni stopni dziesiętnych na
        format stopnie minuty sekundy [st mm ss]

        Parameters
        ----------
        deg : float
            Decimal degree to be converted.

        Returns
        -------
        tuple
            A tuple (degrees, minutes, seconds) where:
            - degrees is an integer,
            - minutes is an integer,
            - seconds is a float.
        """
        negative = deg < 0
        deg = abs(deg)
        d = floor(deg)
        rem = (deg - d) * 60
        m = floor(rem)
        s = (rem - m) * 60

        if negative:
            d = -d

        return (d, m, s)

    def xyz2neu(self):
        """
        Trnasformacja współrzędnych do układu współrzędnych horyzontalnych
        north, east, up (n, e, u) z współrzędnych ortokartezjańskich (x, y, z).
        Transformacja odbywa się przy pomocy macierzy obrotu gdzie po przemnożeniu
        współrzędnych x, y, z, dostajemy współrzędne n, e u.

        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie ortokartezjańskim,

        Returns
        -------
        N, E, U : FLOAT
            współrzędne w ukladzie north, east, up
        """

        results = []
        for x, y, z in self.xyz:
            p = sqrt(x ** 2 + y ** 2)
            f = atan(z / (p * (1 - self.ecc2)))

            while True:
                Rn = self.a / (sqrt(1 - self.ecc2 * sin(f) ** 2))
                h = (p / cos(f)) - Rn
                fs = f
                f = atan(z / (p * (1 - (self.ecc2 * (Rn / (Rn + h))))))
                if abs(fs - f) < (0.000001 / 206265):  # zamienione sekundy na radiany
                    break
            l = atan2(y, x)

            R = array([[-sin(f) * cos(l), -sin(l), cos(f) * cos(l)],
                       [-sin(f) * sin(l), cos(l), sin(l) * cos(f)],
                       [cos(f), 0, sin(f)]])

            dX = array([x, y, z])
            dx = R.T @ dX
            n = dx[0]
            e = dx[1]
            u = dx[2]

            results.append((n, e, u))

        self.write2file("xyz2neu", results)

        return results