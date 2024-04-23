import sys

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
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2
        
        self.xyz = self.readfile(file)
        self.plh = self.xyz2plh(False)


    
    def xyz2plh(self, data_write, output = 'dec_degree'):
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
            r = sqrt(X**2 + Y**2)
            lat_prev = atan(Z / (r * (1 - self.ecc2)))
            lat = 0
            while abs(lat_prev - lat) > 0.000001/206265:
                lat_prev = lat
                N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
                h = r / cos(lat_prev) - N
                lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
            lon = atan(Y/X)
            N = self.a / sqrt(1 - self.ecc2 * sin(lat)**2)
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
        
        # r   = sqrt(X**2 + Y**2)           # promień
        # lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        # lat = 0
        # while abs(lat_prev - lat) > 0.000001/206265:    
        #     lat_prev = lat
        #     N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
        #     h = r / cos(lat_prev) - N
        #     lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        # lon = atan(Y/X)
        # N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        # h = r / cos(lat) - N       
        
        # if output == 'dec_degree':
        #     result = (degrees(lat), degrees(lon), h)
        # elif output == 'dms':
        #     result = (self.deg2dms(degrees(lat)), self.deg2dms(degrees(lon)), h)
        # else:
        #     raise NotImplementedError(f"{output} - output format not defined")
        
        # if data_write is True:
        #     write2file("xyz2plh", result)
        # return result
        
        # if output == "dec_degree":
        #     return degrees(lat), degrees(lon), h 
        # elif output == "dms":w
        #     lat = self.deg2dms(degrees(lat))
        #     lon = self.deg2dms(degrees(lon))
        #     return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        # else:
        #     raise NotImplementedError(f"{output} - output format not defined")

    # def deg2dms(self, x):
    #     sig = ' '
    #     if x < 0:
    #         sig = '-'
    #         x = abs(x)
    #     x = x * 180 / pi
    #     d = int(x)
    #     m = int(60 * (x - d))
    #     s = (x - d - m/60) * 3600
    #     print(f'{sig}{d:3d}{chr(176)}{abs(m):2d}\'{abs(s):7.5f}\"')
    
    def deg2dms(self, deg):
        """
        Converts decimal degrees to degrees, minutes, and seconds.
    
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
        if transformation == "xyz2plh":
            with open('result_xyz2plh.txt', 'w') as file:
                if isinstance(data, tuple):  # Handle single tuple data
                    data = [data]  # Convert to list for uniform processing
                for lat, lon, h in data:
                    if isinstance(lat, tuple) and isinstance(lon, tuple):  # Handling 'dms' format
                        lat_str = f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}"
                        lon_str = f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}"
                        line = f"{lat_str} {lon_str} {h:.3f}\n"
                    else:  # Handling 'dec_degree' format
                        line = f"{lat:.8f} {lon:.8f} {h:.3f}\n"
                    file.write(line)

    def readfile(self, filename):
        """
        Reads a file with coordinates in the format X,Y,Z where each is separated by a comma.
    
        Parameters
        ----------
        filename : str
            The name of the file to read.
    
        Returns
        -------
        coordinates : list of tuples
            A list where each tuple contains (X, Y, Z) as floats.
        """
        coordinates = []
        with open(filename, 'r') as file:
            for line in file:
                if line.strip():  # Ensure the line is not empty
                    parts = line.strip().split(',')
                    if len(parts) == 3:  # Ensure each line has exactly three values
                        try:
                            X = float(parts[0])
                            Y = float(parts[1])
                            Z = float(parts[2])
                            coordinates.append((X, Y, Z))
                        except ValueError as e:
                            print(f"Error converting line to floats: {line}. Error: {e}")
                    else:
                        print(f"Incorrect number of coordinates in line: {line}")
        return coordinates


if __name__ == "__main__":
    
    print("Pick model (one of: wgs84, grs80, mars)")
    model = input("")
    print("Choose transformation (one of: xyz2plh)")
    transformation = input("")
    print("Provide source file")
    file = input("")
    print("")
    
    geo = Transformacje(file, model)
    
    if transformation == "xyz2plh":
        print("Pick data output type for XYZ -> BLH transformation (one of: dec_degree, dms)")
        output_type = input("")
        print("")
        geo.xyz2plh(True, output_type)
        print("Transformation of XYZ -> BLH saved to result_xyz2plh.txt")
    else:
        print("Wrong transformation type")
    



   

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    