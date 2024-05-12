# Projekt_1 Informatyka Geodezyjna

## Cel i funkcjonalność programu

Program służy do transformacji współrzędnych pomiędzy różnymi układami. Podając współrzędne w układzie X,Y,Z skrypt pozwala na następujące transformacje:

- Transformacja X,Y,Z do fi, lambda, h (xyz2plh)
- Transformacja fi, lambda, h do X,Y,Z (plh2xyz)
- Transformacja X,Y,Z do NEU (xyz2neu)
- Transformacja X,Y,Z do układu 2000 (fl22000)
- Transformacja X,Y,Z do układu 1992 (fl21992)

Obsługiwane elipsoidy:
- **WGS'84** 
- **GRS'80**
- **Mars**

## Wymagania techniczne

Aby program poprawnie działał na docelowej maszynie, muszą być spełnione następujące warunki:

- Zainstalowany Python w wersji **3.12** 
- Zainstalowana biblioteka **numpy** (_pip install numpy_)
- Posiadać system operacyjny wspierający **CLI** (command line interface)

### Systemy operacyjne

Program działa na następujących systemach operacyjnych:

- **macOS Sonoma** (wersja 14.4.1)
- **Windows 11** (wersja 22H2, 23H2)

## Użycie

Aby mieć możliwość skorzystania z programu należy sklonować to repozytorium na docelową maszynę, a następnie uruchomić CLI w sklonowanym folderze.

Program można użyć na dwa sposoby:

### Interaktywne CLI

Interaktywne CLI polega na wywołaniu programu przy pomocy komendy **python main.py**, która skompiluje kod interaktywnie pytający o wymagane informacje potrzebne do wykonywania obliczeń i transformacji.

1. Program poprosi nas o podanie jednej z dostępnych elipsoid (_**modelu**_)
>Wybierz model elipsoidy (jeden z: wgs84, grs80, mars)
2. Program poprosi nas o podanie transformacji, którą chcemy przeprowadzić
>Wybierz transformacje (jedną z: xyz2plh, plh2xyz, xyz2neu, fl22000, fl21992)
3. Program poprosi nas o podanie nazwy pliku źródłowego, jeżeli znajduje się on w tym samym folderze co program, lub pełną ścieżkę do pliku, jeżeli plik znajduje się w innym folderze
>Podaj nazwę pliku/ścieżkę do pliku z danymi źródłowymi:
  
_W przypadku wybrania transformacji **xyz2plh** wymagane jest podanie trybu zapisania przeliczonych danych (stopnie lub stopnie dziesiętne)_
>Podaj typ zapisu danych w wypadku transformacji XYZ -> BLH (jeden z: dec_degree, dms)

_W przypadku wybrania transformacji **xyz2neu** wymagane jest podanie 3 dodatkowych współrzędnych w postaci float (np. **3664940.500**): **x0**, **y0** i **z0**_
>Podaj wartości x0, y0, z0 dla transformacji XYZ -> NEUp

Dane po przeliczeniu zostaną zapisane do pliku **.txt** w tym samym folderze gdzie program został uruchomiony z nazwą zależną od przeprowadzonej transformacji:

- xyz2plh - **result_xyz2plh.txt**
- plh2xyz - **result_plh2xyz.txt**
- xyz2neu - **result_xyz2neu.txt**
- fl22000 - **result_fl22000.txt**
- fl21992 - **result_fl21992.txt**

#### Przykładowe użycie programu:

1. Przechodzimy w CLI do folderu z programem przy użyciu komendy cd
>cd Projekt_1
2. Wywołujemy program
>python main.py
3. Wybieramy elipsoidę wgs84
>wgs84
4. Wybieramy transformację xyz2plh
>xyz2plh
5. Wybieramy zapisanie transformacji w formie stopni dziesiętnych
>dec_degree

Plik result_xyz2plh.txt zostaje zapisany do folderu gdzie został wywołany program z następującymi wartościami
>52.09727222 21.03153333 141.399  
>52.09727216 21.03153314 141.400  
>52.09727212 21.03153296 141.403  
>52.09727209 21.03153277 141.408  
>52.09727209 21.03153323 141.410  
>52.09727212 21.03153318 141.402  
>52.09727207 21.03153300 141.406  
>52.09727206 21.03153281 141.411  
>52.09727212 21.03153325 141.407  
>52.09727214 21.03153318 141.405  
>52.09727210 21.03153332 141.408  
>52.09727215 21.03153318 141.406
>> Wartości są rozdzielone znakiem " " (whitespace)


### CLI z flagami (wykonanie w jednej linii)

CLI z flagami pozwala na wywołanie programu przy pomocy tak zwanych _flag_, które pozwalają na ominięcie całęgo interaktywnego procesu i podanie wszystkich wymaganych informacji jednocześnie (w jednej linii).

Aby uruchomić program w ten sposób należy do komendy wywołującej dodać flagę **--use_cli**
>python main.py --use_cli (reszta flag)

Wywołując program przy użyciu tej flagi, należy podać wszystkie pozostałe flagi, inaczej otrzymamy błąd:

- **--model** (wgs84, grs80, mars)
>--model **wgs84** | --model **grs80** | --model **mars**
- **--transformation** (xyz2plh, plh2xyz, xyz2neu, fl22000, fl21992)
>--transformation **xyz2plh** | --transformation **plh2xyz** | --transformation **xyz2neu** | --transformation **fl22000** | --transformation **fl21992**
- **--file** (nazwa pliku / ścieżka do pliku)
>--file **file.txt**

_W przypadku wybrania transformacji **xyz2plh** można podać tryb zapisania przeliczonych danych **--output_type** (stopnie - **dms** lub stopnie dziesiętne - **dec_degree**), domyślnie (jeżeli flagi nie podamy) ta wartość to stopnie dziesiętne - **dec_degree**_
>--output_type **dec_degree** | --output_type **dms**

_W przypadku wybrania transformacji **xyz2plh** wymagane jest podanie 3 dodatkowych współrzędnych w postaci float (np. **3664940.500**): **x0**, **y0** i **z0**_
>--neu **x0_value y0_value z0_value**

#### Przykładowe użycie programu:

1. Przechodzimy w CLI do folderu z programem przy użyciu komendy cd
>cd Projekt_1
2. Wywołujemy program w modelu **wgs84**, z transformacją **xyz2plh** w formacie **stopni dziesiętnych** (_dec_degree_), używając pliku źródłowego **input_xyz.txt**
>python main.py --use_cli --model wgs84 --transformation xyz2plh --file input_xyz.txt --output_type dec_degree

Plik **result_xyz2plh.txt** zostaje zapisany do folderu gdzie został wywołany program z następującymi wartościami
>52.09727222 21.03153333 141.399  
>52.09727216 21.03153314 141.400  
>52.09727212 21.03153296 141.403  
>52.09727209 21.03153277 141.408  
>52.09727209 21.03153323 141.410  
>52.09727212 21.03153318 141.402  
>52.09727207 21.03153300 141.406  
>52.09727206 21.03153281 141.411  
>52.09727212 21.03153325 141.407  
>52.09727214 21.03153318 141.405  
>52.09727210 21.03153332 141.408  
>52.09727215 21.03153318 141.406
>> Wartości są rozdzielone znakiem " " (whitespace)

### Plik z danymi źródłowymi

Plik z danymi źródłowymi powinien być zapisany w układzie współrzędnych XYZ w następujący sposób (oddzielenie wartości znakiem "**,**")
>wartosc_X,wartosc_Y,wartosc_Z

Przykład (**_input_xyz.txt_**):
>3664940.500,1409153.590,5009571.170
>3664940.510,1409153.580,5009571.167
>3664940.520,1409153.570,5009571.167
>3664940.530,1409153.560,5009571.168
>3664940.520,1409153.590,5009571.170
>3664940.514,1409153.584,5009571.166
>3664940.525,1409153.575,5009571.166
>3664940.533,1409153.564,5009571.169
>3664940.515,1409153.590,5009571.170
>3664940.514,1409153.584,5009571.169
>3664940.515,1409153.595,5009571.169
>3664940.513,1409153.584,5009571.171
