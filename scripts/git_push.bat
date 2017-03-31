@echo off
REM Script by Sean P

echo Comitting updates...
echo git add --all
cd .. & "C:\Program Files (x86)\SmartGit\git\bin\git" add --all
echo git commit -a
set /p input="Enter your commit message: "
"C:\Program Files (x86)\SmartGit\git\bin\git" commit -a -m "%input%"
echo Pushing updates to repository...
echo git push --all
"C:\Program Files (x86)\SmartGit\git\bin\git" push --all

pause