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

    def fl22000(self):
        """
        Algorytm transformacji szerokoci i wysokoci geodezyjnej (phi, lam) do układu współrzędnych
        płaskiech prostokątnych PL-2000 (x-2000, y-2000).
        Jest to proces wykorzystujący współrzędne w odwzorowaniu G-K.
        Układ PL-2000 wykorzystuje 4 pasy południkowe o szerokoci 3 stopni (15, 18, 21, 24 [st])

        Parameters
        ----------
        phi : FLOAT
            [stopnie dziesiętne]
        lam : FLOAT
            [stopnie dziesiętne]

        Returns
        -------
        X, Y : FLOAT
            Współrzedne w układzie 2000

        """

        m0_2000 = 0.999923

        b2 = (self.a ** 2) * (1 - self.ecc2)
        ep2 = (self.a ** 2 - b2) / b2

        results = []

        for phi, lam, h in self.plh:

            phi = radians(phi)
            lam = radians(lam)

            if lam < radians(16.5):
                nr = 5
                l_20 = 15 * pi / 180
            elif lam < radians(19.5) and lam > radians(16.5):
                nr = 6
                l_20 = 18 * pi / 180
            elif lam < radians(22.5) and lam > radians(19.5):
                nr = 7
                l_20 = 21 * pi / 180
            elif lam > radians(22.5):
                nr = 8
                l_20 = 24 * pi / 180

            dl = lam - l_20
            t = tan(phi)

            eta2 = ep2 * cos(phi) ** 2

            Rn = self.a / (sqrt(1 - self.ecc2 * sin(phi) ** 2))

            A0 = 1 - (self.ecc2 / 4) - (3 * self.ecc2 ** 2 / 64) - 5 * self.ecc2 ** 3 / 256
            A2 = (3 / 8) * (self.ecc2 + (self.ecc2 ** 2) / 4 + ((15 * self.ecc2 ** 3) / 128))
            A4 = (15 / 256) * ((self.ecc2 ** 2) + ((3 * self.ecc2 ** 3) / 4))
            A6 = (35 * (self.ecc2 ** 3)) / 3072

            sigma = self.a * (A0 * phi - A2 * sin(2 * phi) + A4 * sin(4 * phi) - A6 * sin(6 * phi))

            x_GK = sigma + (dl ** 2 / 2) * Rn * sin(phi) * cos(phi) * (
                        1 + (dl ** 2 / 12) * cos(phi) ** 2 * (5 - t ** 2 + 9 * eta2 + 4 * eta2 ** 2) + (
                            dl ** 4 / 360) * cos(phi) ** 4 * (61 - 58 * t ** 2 + 270 * eta2 - 330 * eta2 * t ** 2))
            y_GK = dl * Rn * cos(phi) * (
                        1 + (dl ** 2 / 6) * cos(phi) ** 2 * (1 - t ** 2 + eta2) + (dl ** 4 / 120) * cos(phi) ** 4 * (
                            5 - 18 * t ** 2 + t ** 4 + 14 * eta2 - 58 * eta2 * t ** 2))

            x_2000 = x_GK * m0_2000
            y_2000 = y_GK * m0_2000 + nr * 1000000 + 500000

            results.append((x_2000, y_2000))

        self.write2file("fl22000", results)

        return results

    def fl21992(self):
        """
        Algorytm transformacji szerokoci i wysokoci geodezyjnej (phi, lam) do układu współrzędnych
        płaskiech prostokątnych PL-1992 (x-1992, y-1992).
        Jest to proces wykorzystujący współrzędne w odwzorowaniu G-K.
        Układ PL-1992 wykorzystuje 1 pas południkowokwy 19 stopni.


        Parameters
        ----------
        phi : FLOAT
            [stopnie dziesiętne]
        lam : FLOAT
            [stopnie dziesiętne]

        Returns
        -------
        X, Y  : FLOAT
            Współrzedne w układzie 1992

        """

        l_92 = 19 * pi / 180
        m0_1992 = 0.9993

        b2 = (self.a ** 2) * (1 - self.ecc2)
        ep2 = (self.a ** 2 - b2) / b2

        results = []

        for phi, lam, h in self.plh:
            phi = radians(phi)
            lam = radians(lam)

            dl = lam - l_92
            t = tan(phi)

            eta2 = ep2 * cos(phi) ** 2

            Rn = self.a / (sqrt(1 - self.ecc2 * sin(phi) ** 2))

            A0 = 1 - (self.ecc2 / 4) - (3 * self.ecc2 ** 2 / 64) - 5 * self.ecc2 ** 3 / 256
            A2 = (3 / 8) * (self.ecc2 + (self.ecc2 ** 2) / 4 + ((15 * self.ecc2 ** 3) / 128))
            A4 = (15 / 256) * ((self.ecc2 ** 2) + ((3 * self.ecc2 ** 3) / 4))
            A6 = (35 * (self.ecc2 ** 3)) / 3072

            sigma = self.a * (A0 * phi - A2 * sin(2 * phi) + A4 * sin(4 * phi) - A6 * sin(6 * phi))

            x_GK = sigma + (dl ** 2 / 2) * Rn * sin(phi) * cos(phi) * (
                        1 + (dl ** 2 / 12) * cos(phi) ** 2 * (5 - t ** 2 + 9 * eta2 + 4 * eta2 ** 2) + (
                            dl ** 4 / 360) * cos(phi) ** 4 * (61 - 58 * t ** 2 + 270 * eta2 - 330 * eta2 * t ** 2))
            y_GK = dl * Rn * cos(phi) * (
                        1 + (dl ** 2 / 6) * cos(phi) ** 2 * (1 - t ** 2 + eta2) + (dl ** 4 / 120) * cos(phi) ** 4 * (
                            5 - 18 * t ** 2 + t ** 4 + 14 * eta2 - 58 * eta2 * t ** 2))

            x_1992 = x_GK * m0_1992 - 5300000
            y_1992 = y_GK * m0_1992 + 500000

            results.append((x_1992, y_1992))

        self.write2file("fl21992", results)

        return results

    def write2file(self, transformation, data):
        """
        Writes coordinate data to a file.

        Parameters
        ----------
        filename : str
            The name of the file where the data will be written.
        data : list of tuples
            A list of tuples containing coordinate data (latitude, longitude, height).
        mode : str, optional
            The file opening mode, 'w' for write (default) or 'a' for append.
        """
        filename = ''
        if transformation == "xyz2plh":
            filename = 'result_xyz2plh.txt'
        elif transformation == "plh2xyz":
            filename = 'result_plh2xyz.txt'
        elif transformation == "xyz2neu":
            filename = 'result_xyz2neu.txt'
        elif transformation == "fl22000":
            filename = 'result_fl22000.txt'
        elif transformation == "fl21992":
            filename = 'result_fl21992.txt'
        else:
            raise ValueError("Unsupported transformation type")

        with open(filename, "w") as file:
            if isinstance(data, tuple):  # Handle single tuple data
                data = [data]  # Convert to list for uniform processing

            if transformation == "xyz2plh":
                for lat, lon, h in data:
                    if isinstance(lat, tuple) and isinstance(lon, tuple):  # Handling 'dms' format
                        lat_str = f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}"
                        lon_str = f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}"
                        line = f"{lat_str} {lon_str} {h:.3f}\n"
                    else:  # Handling 'dec_degree' format
                        line = f"{lat:.8f} {lon:.8f} {h:.3f}\n"
                    file.write(line)
            elif transformation == "plh2xyz":
                for X, Y, Z in data:
                    line = f"{X:.3f} {Y:.3f} {Z:.3f}\n"
                    file.write(line)
            elif transformation == "xyz2neu":
                for n, e, u in data:
                    line = f"{n:.3f} {e:.3f} {u:.3f}\n"
                    file.write(line)
            elif transformation == "fl22000":
                for x_2000, y_2000 in data:
                    line = f"{x_2000:.3f} {y_2000:.3f}\n"
                    file.write(line)
            elif transformation == "fl21992":
                for x_1992, y_1992 in data:
                    line = f"{x_1992:.3f} {y_1992:.3f}\n"
                    file.write(line)