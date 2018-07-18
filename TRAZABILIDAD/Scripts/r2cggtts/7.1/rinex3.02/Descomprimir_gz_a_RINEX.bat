ECHO ON

REM -------------DESCRIPCIÓ--------------
REM Bat que descomprimeix els fitxers 
REM RINEX 3 amb compressió binaria(*.crx.gz)
REM a format RINEX 3 (*.rnx). Elimina
REM els fitxers en format RINEX-Hatanaka.
REM
REM Crear un directori nou i copiar el fitxer bat, 
REM el gzip.exe, el crx2rnx.exe i tot els fitxers 
REM a descomprimir en el nou directori.
REM Accedir al directori creat i executar el 
REM fitxer bat. 
REM
REM -------------------------------------

REM --------FITXER NECESSARIS------------

REM gzip.exe -> Fitxer per descomprimir de binary a Hatanaka (http://www.gzip.org/#exe)

REM crx2rnx.exe -> Fitxer per descomprimir de Hatanaka a RINEX (http://terras.gsi.go.jp/ja/crx2rnx.html)

REM ------------------------------------- 


REM --------- Descomprimim els fitxers amb el GZIP ------------

For %%x in (.\??????????????????????????????????.crx.gz) do (

	.\gzip.exe -d %%x	
	
)


REM --------- Descomprimim els fitxers a RINEX 3 ------------


For %%x in (.\??????????????????????????????????.crx) do (

	.\crx2rnx.exe -f %%x
	
)

REM --------- Esborrem els fitxer RINEX 3 amb compressió Hatanaka ------------


DEL .\??????????????????????????????????.crx

