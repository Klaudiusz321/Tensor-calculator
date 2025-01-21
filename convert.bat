@echo off
REM -----------------------------------------------------
REM Skrypt: przetwarzanie.bat
REM Opis  : Auto-crop + (ew. invert) w IrfanView (tryb advancedbatch)
REM -----------------------------------------------------

:: 1. Ścieżka do IRFANVIEW (plik EXE, nie skrót!)
SET "IVPATH=C:\Program Files\IrfanView\i_view64.exe"

:: 2. Folder z oryginałami (uwzględnij rozszerzenie, np. *.png)
SET "INPATH=C:\Users\sorak\Desktop\screeny praca\zwykle\*.png"

:: 3. Folder docelowy z formatem wyjściowym
SET "OUTPATH=C:\Users\sorak\Desktop\screeny praca\latex\*.png"

:: 4. Wywołanie IrfanView w trybie advancedbatch:
"%IVPATH%" "%INPATH%" /advancedbatch /convert="%OUTPATH%"

echo -----------------------------------
echo Przetwarzanie zakonczone!
echo Znajdziesz pliki w: C:\Users\sorak\Desktop\screeny praca\latex
echo -----------------------------------
pause
