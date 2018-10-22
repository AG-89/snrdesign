@echo off

set MATLAB=C:\Assorted\Matlab

cd .

if "%1"=="" ("C:\Assorted\Matlab\bin\win64\gmake"  -f findpeaks_fast_rtw.mk all) else ("C:\Assorted\Matlab\bin\win64\gmake"  -f findpeaks_fast_rtw.mk %1)
@if errorlevel 1 goto error_exit

exit /B 0

:error_exit
echo The make command returned an error of %errorlevel%
An_error_occurred_during_the_call_to_make
