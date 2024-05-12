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
        self.model = model

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

    def get_model(self):
        return self.model

    def xyz2plh(self, data_write, output='dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokość elipsoidalna (phi, lam, h). Jest to proces iteracyjny.
        W wyniku 3-4-krotnej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnością ok. 1 cm.
        Parametry
        ----------
        X, Y, Z :
        [float]
             współrzędne w układzie orto-kartezjańskim,

        Zwracane wartości:
        -------
        results:
            [lista wartości]:
                lat:
                    [stopnie dziesiętne] - szerokość geodezyjna
                lon:
                    [stopnie dziesiętne] - długość geodezyjna.
                h :
                    [metry] - wysokość elipsoidalna
        output [plik .txt] :
            dec_degree - stopnie dziesiętne
            dms - stopnie, minuty, sekundy
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
        Konwertuje współrzędne geodezyjne (phi, lam, h) na współrzędne kartezjańskie (X, Y, Z).

        Parametry
        ----------
        Phi, lam, h:
        [lista krotek]
            Każda krotka zawiera (phi, lam, h) w stopniach dziesiętnych i metrach.

        Zwracane wartości
        -------
        X, Y, Z:
        [lista krotek]
            Każda krotka zawiera (X, Y, Z) jako liczby zmiennoprzecinkowe, zaokrąglone do 3 miejsc po przecinku.
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
        Zamiana stopni dziesiętnych na
        format stopnie minuty sekundy [st mm ss]

        Parametery
        ----------
        deg :
        [float]
            Stopnie dziesiętne do przeliczenia.

        Zwracane wartości
        -------
        dms :
        [krotka]
            Krotka (stopnie, minuty, sekundy), w której:
                - stopnie to liczba całkowita,
                - minuty to liczba całkowita,
                - sekundy to liczba zmiennoprzecinkowa.
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

    def xyz2neu(self, x0, y0, z0):
        """
        Transformacja współrzędnych do układu współrzędnych horyzontalnych
        north, east, up (n, e, u) z współrzędnych ortokartezjańskich (x, y, z).
        Transformacja odbywa się przy pomocy macierzy obrotu gdzie po przemnożeniu
        współrzędnych x, y, z, dostajemy współrzędne n, e u.

        Parametery
        ----------
        X, Y, Z :
        [float]
             współrzędne w układzie ortokartezjańskim.

        Zwracane wartości
        -------
        N, E, U :
        [float]
            współrzędne w układzie north, east, up.
        """
        # x0 = 3664940.500
        # y0 = 1409153.590
        # z0 = 5009571.170

        x0 = float(x0)
        y0 = float(y0)
        z0 = float(z0)
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

            #dX = array([x, y, z])
            dX = array([x-x0, y-y0, z-z0])
            dx = R.T @ dX
            n = dx[0]
            e = dx[1]
            u = dx[2]

            results.append((n, e, u))

        # print(results)
        self.write2file("xyz2neu", results)

        return results

    def fl22000(self):
        """
        Algorytm transformacji szerokości i wysokości geodezyjnej (phi, lam) do układu współrzędnych
        płaskich prostokątnych PL-2000 (x-2000, y-2000).
        Jest to proces wykorzystujący współrzędne w odwzorowaniu G-K.
        Układ PL-2000 wykorzystuje 4 pasy południkowe o szerokoci 3 stopni (15, 18, 21, 24 [st])

        Parametry
        ----------
        phi :
        [float]
            [stopnie dziesiętne]
        lam :
        [float]
            [stopnie dziesiętne]

        Zwracane wartości
        -------
        X, Y :
        [float]
            Współrzędne w układzie 2000.

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
        Algorytm transformacji szerokości i wysokości geodezyjnej (phi, lam) do układu współrzędnych
        płaskich prostokątnych PL-1992 (x-1992, y-1992).
        Jest to proces wykorzystujący współrzędne w odwzorowaniu G-K.
        Układ PL-1992 wykorzystuje 1 pas południkowy 19 stopni.


        Parametery
        ----------
        phi :
        [float]
            [stopnie dziesiętne]
        lam :
        [float]
            [stopnie dziesiętne]

        Zwracane wartości
        -------
        X, Y  :
        [float]
            Współrzędne w układzie 1992.

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
        Zapisuje współrzędne do pliku.

        Parametery
        ----------
        filename :
        [string]
            Nazwa pliku, w którym dane zostaną zapisane.
        data :
        [lista krotek]
            Lista krotek zawierających współrzędne.
        Zwracane wartości
        -------
        Plik z danymi po transformacji
        """
        filename = ''
        if transformation == "xyz2plh":
            filename = f'result_xyz2plh_{self.model}.txt'
        elif transformation == "plh2xyz":
            filename = f'result_plh2xyz_{self.model}.txt'
        elif transformation == "xyz2neu":
            filename = f'result_xyz2neu_{self.model}.txt'
        elif transformation == "fl22000":
            filename = f'result_fl22000_{self.model}.txt'
        elif transformation == "fl21992":
            filename = f'result_fl21992_{self.model}.txt'
        else:
            raise ValueError("Unsupported transformation type")

        with open(filename, "w") as file:
            if isinstance(data, tuple):  # Przypadek danych będących jedną krotką
                data = [data]  # Przekonwertowanie do listy dla ułatwienia uniwersalnego zapisu danych do pliku

            if transformation == "xyz2plh":
                for lat, lon, h in data:
                    if isinstance(lat, tuple) and isinstance(lon, tuple):  # Przypadek danych w formacie 'dms'
                        lat_str = f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}"
                        lon_str = f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}"
                        line = f"{lat_str} {lon_str} {h:.3f}\n"
                    else:  # Przypadek danych w formacie 'dec_degree'
                        line = f"{lat:.8f} {lon:.8f} {h:.3f}\n"
                    file.write(line)
            elif transformation == "plh2xyz":
                for X, Y, Z in data:
                    line = f"{X:.3f} {Y:.3f} {Z:.3f}\n"
                    file.write(line)
            elif transformation == "xyz2neu":
                for index, (n, e, u) in enumerate(data):
                    if index == 0:
                        line = f"{n} {e} {u}\n"
                    else:
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

    def readfile(self, filename):
        """
        Odczytuje plik ze współrzędnymi w formacie X,Y,Z, gdzie każda z nich oddzielona jest przecinkiem.

        Parametery
        ----------
        filename :
        [string]
            Nazwa pliku do odczytu/ ścieżka do pliku.

        Zwracane wartości
        -------
        coordinates :
        [lista krotek]
            Lista, w której każda krotka zawiera (X, Y, Z) jako liczby zmiennoprzecinkowe.
        """
        coordinates = []
        with open(filename, 'r') as file:
            for line in file:
                if line.strip():  # Upewnienie się, że linia w pliku nie jest pusta
                    parts = line.strip().split(',')
                    if len(parts) == 3:  # Upewnienie się, że linia ma dokładnie 3 wartości
                        try:
                            X = float(parts[0])
                            Y = float(parts[1])
                            Z = float(parts[2])
                            coordinates.append((X, Y, Z))
                        except ValueError as e:
                            print(f"Wystąpił błąd w trakcie przemiany linii do float : {line}. Error: {e}")
                    else:
                        print(f"Niepoprawna liczba współrzędnych w linii: {line}")
        return coordinates


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Zarządzanie narzędziem transformacji (w ramach CLI).")

    # Flaga decydująca o tym, czy program będzie interaktywny, czy będzie go można wywołać w 1 linijce
    parser.add_argument('--use_cli', action='store_true', help='Ustaw ten znacznik, aby używać argumentów wiersza poleceń zamiast interaktywnych komunikatów.')

    # Sczytaj wartości z flag
    parser.add_argument('--model', type=str, choices=['wgs84', 'grs80', 'mars'], help='Wybierz model elipsoidy.')
    parser.add_argument('--transformation', type=str, choices=['xyz2plh', 'plh2xyz', 'xyz2neu', 'fl22000', 'fl21992'], help='Wybierz model transformacji.')
    parser.add_argument('--file', type=str, help='Podaj nazwę/ ścieżkę do pliku z danymi źródłowymi.')
    parser.add_argument('--output_type', type=str, choices=['dec_degree', 'dms'], default='dec_degree', help='Podaj typ zapisu danych w wypadku transformacji xyz2plh.')
    parser.add_argument('--neu', type=float, nargs=3, help='Podaj wartości x0, y0 i z0 dla transformacji xyz2neu, w przeciwnym wypadku transformacja nie dojdzie do skutku.')

    # Sparsuj sczytane wartości
    args = parser.parse_args()

    if args.use_cli:

        # Upewnij się, że wszystkie wartości flag zostały podane przez użytkownika
        if not all([args.model, args.transformation, args.file]):
            parser.error("Przy podaniu flagi --use_cli wymagane są flagi: --model, --transformation, --file.")

        # Stwórz obiekt transformacji z wymaganym przez obiekt plikiem źródłowym i elipsoidą (modelem)
        geo = Transformacje(file=args.file, model=args.model)

        # Wykonaj wybraną transformację
        if args.transformation == "xyz2plh":
            geo.xyz2plh(True, args.output_type)
            print(f"Transformacja XYZ -> BLH została poprawnie zapisana do pliku: result_xyz2plh_{args.model}.txt")
        elif args.transformation == "plh2xyz":
            geo.plh2xyz()
            print(f"Transformacja BLH -> XYZ została poprawnie zapisana do pliku: result_plh2xyz_{args.model}.txt")
        elif args.transformation == "xyz2neu":
            x0, y0, z0 = args.neu
            geo.xyz2neu(x0, y0, z0)
            print(f"Transformacja of XYZ -> NEU została poprawnie zapisana do pliku: result_xyz2neu_{args.model}.txt")
        elif args.transformation == "fl22000":
            geo.fl22000()
            print(f"Transformacja of BL -> 2000 została poprawnie zapisana do pliku: result_fl22000_{args.model}.txt")
        elif args.transformation == "fl21992":
            geo.fl21992()
            print(f"Transformacja of BL -> 1992 została poprawnie zapisana do pliku: result_fl21992_{args.model}.txt")
        else:
            print("Podany typ transformacji jest niepoprawny.")

    else:

        print("Wybierz model elipsoidy (jeden z: wgs84, grs80, mars)")
        model = input("")
        print("Wybierz transformacje (jedną z: xyz2plh, plh2xyz, xyz2neu, fl22000, fl21992)")
        transformation = input("")
        print("Podaj nazwę pliku/ścieżkę do pliku z danymi źródłowymi:")
        file = input("")
        print("")

        geo = Transformacje(file, model)

        if transformation == "xyz2plh":
            print("Podaj typ zapisu danych dla transformacji XYZ -> BLH (jeden z: dec_degree, dms)")
            output_type = input("")
            print("")
            geo.xyz2plh(True, output_type)
            print(f"Transformacja of XYZ -> BLH została poprawnie zapisana do pliku: result_xyz2plh_{model}.txt")
        elif transformation == "plh2xyz":
            geo.plh2xyz()
            print(f"Transformacja of BLH -> XYZ została poprawnie zapisana do pliku: result_plh2xyz_{model}.txt")
        elif transformation == "xyz2neu":
            print("Podaj wartości x0, y0, z0 dla transformacji XYZ -> NEUp")
            x0 = input("x0: ")
            y0 = input("y0: ")
            z0 = input("z0: ")
            print("")
            geo.xyz2neu(x0, y0, z0)
            print(f"Transformacja of XYZ -> NEUp została poprawnie zapisana do pliku: result_xyz2neu_{model}.txt")
        elif transformation == "fl22000":
            geo.fl22000()
            print(f"Transformacja of BL -> 2000 została poprawnie zapisana do pliku: result_fl22000_{model}.txt")
        elif transformation == "fl21992":
            geo.fl21992()
            print(f"Transformacja of BL -> 1992 została poprawnie zapisana do pliku: result_fl21992_{model}.txt")
        else:
            print("Podany typ transformacji jest niepoprawny.")
