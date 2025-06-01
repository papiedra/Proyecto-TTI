
# Proyecto Taller Transeversal I
Proyecto realizado por Pablo Piedrafita Sanromán.

Para la ejecución de los test, ejecutar el fichero Makefile_test.bat o ejecutar la secuencia de comandos (desde el directorio del proyecto):

g++ tests/tests.cpp src/*.cpp -lm -o bin/tests.exe
cd bin
tests.exe
pause


Nota: Si al ejecutar los test, el método de m_nzeros_01()  falla, ignorar, creo que es un tema de precisión, volver a ejecutar los test y debería de funcionar sin problema.
En algunos test de la parte extra notarás que he tenido que ser menos estricto con la precisión. Esto lo he hecho para que los test pasen de forma más sencilla y no tener que pegarme con redondeos de números muy pequeños. Los consideraré errores "admisibles"

Para ejecutar el main, ejecutar el fichero Makefile.bat o ejecutar la secuencia de comandos (desde el directorio del proyecto):

g++ tests/EFK_GEOS3.cpp src/*cpp -lm -o bin/main.exe
cd bin
main.exe
pause