@echo off

setlocal ENABLEEXTENSIONS

REM   Amount of memory to allocate to RTG.  Use G suffix for gigabytes
REM   If not set, will default to 90% of available RAM
REM   for example following sets maximum memory to 4 GB
REM   SET RTG_MEM=4G

REM   Proxy host if needed, blank otherwise
SET PROXY_HOST=

REM   Proxy port if needed, blank otherwise
SET PROXY_PORT=

REM   Additional JVM options (e.g.: "-Djava.io.tmpdir=XXYY -XX:+UseLargePages")
SET RTG_JAVA_OPTS=

REM   Maximum memory for rtg to use (e.g. 48g)
SET RTG_MEM=

REM   If RTG_MEM is not defined use this percentage of total RAM
SET RTG_MEM_PCT=90

REM   Attempt to send crash logs to realtime genomics, "true" to enable
SET RTG_TALKBACK=

REM   Attempt to send simple usage logs to realtime genomics, "true" to enable
SET RTG_USAGE=

REM   Server URL when usage logging to server
SET RTG_USAGE_HOST=

REM   Destination directory when performing single-user file-based usage logging
SET RTG_USAGE_DIR=

REM   Optional fields to add to usage logging. e.g. username,hostname,commandline
SET RTG_USAGE_OPTIONAL=

REM   Directory in which to look for pipeline reference datasets
SET RTG_REFERENCES_DIR=

REM   Directory in which to look for pipeline reference datasets
SET RTG_MODELS_DIR=

REM   Default number of threads to use
SET RTG_DEFAULT_THREADS=

REM --- END OF USER SETTINGS ---
 
REM   location of RTG jar file
SET RTG_JAR="%~dp0RTG.jar"

REM   minimum memory used for RTG talkback/usage and other calls
SET RTG_MIN_MEM=64M

REM   path to java.exe
SET RTG_JAVA="%~dp0%jre\bin\java.exe"

REM check if we were started from windows explorer and provide a sensible message
echo %cmdcmdline% | find "cmd /c" >nul
if %ERRORLEVEL% NEQ 0 GOTO :start

SET PAUSE_ON_CLOSE=1
IF /I "%PROCESSOR_ARCHITECTURE%" NEQ "AMD64" GOTO :no64bit

echo RTG is a command line application, please follow the steps below
echo.  
echo 1) Open windows command prompt  
echo    e.g. on Windows Vista / 7, press start, search "cmd.exe" and press enter
echo.  
echo 2) Type "%~dp0%rtg.bat" to execute RTG command
echo.  
echo For more information please refer to "%~dp0%RTGOperationsManual.pdf"
echo.  
GOTO :no 
 
:start
IF /I "%PROCESSOR_ARCHITECTURE%" NEQ "AMD64" GOTO :no64bit

REM   Following checks for supported OS
REM   Windows XP 64 bit and Windows Server 2003
ver | FIND "5.2" > nul
IF %ERRORLEVEL% EQU 0 GOTO :os_ok

REM   Windows Vista and Windows Sever 2008 
ver | FIND "6.0" > nul
IF %ERRORLEVEL% EQU 0 GOTO :os_ok

REM   Windows 7 and Windows Sever 2008 R2 
ver | FIND "6.1" > nul
IF %ERRORLEVEL% EQU 0 GOTO :os_ok

REM   Windows 8 and Windows Sever 2012
ver | FIND "6.2" > nul
IF %ERRORLEVEL% EQU 0 GOTO :os_ok

REM   Windows 8.1 and Windows Sever 2012 R2 
ver | FIND "6.3" > nul
IF %ERRORLEVEL% EQU 0 GOTO :os_ok

REM   We have checked for all tested versions of Windows   
echo This is not a supported version of windows, please contact support@realtimegenomics.com for more information
GOTO :no

:os_ok

SET ACCEPTED_NAME="%~dp0%.license_accepted"

SET ACCEPTED_TALKBACK_NAME="%~dp0%.talkback_accepted"

SET ACCEPTED_USAGE_NAME="%~dp0%.usage_accepted"

SET EULA="%~dp0%EULA.txt"

IF EXIST %ACCEPTED_NAME% GOTO :setvars

REM  Skip asking about the EULA if we havn't bundled it
IF NOT EXIST %EULA% GOTO :yes

more %EULA%

SET /P ANSWER=Do you agree to the terms and conditions (y/n)? 

IF /i {%ANSWER%}=={y} (GOTO :yes)
IF /i {%ANSWER%}=={yes} (GOTO :yes)

echo You must agree with the license terms before you can use the software.
GOTO :no

:yes

echo %date% - %time% >%ACCEPTED_NAME%

SET /P TALKBACK_ANSWER=Would you like to enable automatic crash reporting (y/n)? 
echo.
IF /i {%TALKBACK_ANSWER%}=={y} (GOTO :talkback_yes)
IF /i {%TALKBACK_ANSWER%}=={yes} (GOTO :talkback_yes)

GOTO :talkback_no

:talkback_yes
echo true >%ACCEPTED_TALKBACK_NAME%

REM   Try and send a post-install talkback
echo Testing talkback connectivity...
%RTG_JAVA% %PROXY_HOST% %PROXY_PORT% -Xmx%RTG_MIN_MEM% -cp %RTG_JAR% com.rtg.util.diagnostic.SimpleTalkback "Post-install talkback test"
IF %ERRORLEVEL% NEQ 0 (
     echo Initial talkback did not succeed, probably due to firewall issues.
     echo You will be asked to manually submit any error logs.
)
echo.
GOTO :usage

:talkback_no
echo Automatic crash reporting disabled.
echo.
<nul (SET /P foo=false) >%ACCEPTED_TALKBACK_NAME%

:usage
SET /P USAGE_ANSWER=Would you like to enable automatic simple usage reporting (y/n)? 
echo.
IF /i {%USAGE_ANSWER%}=={y} (GOTO :usage_yes)
IF /i {%USAGE_ANSWER%}=={yes} (GOTO :usage_yes)

GOTO :usage_no

:usage_yes
<nul (SET /P FOO=true) >%ACCEPTED_USAGE_NAME%
GOTO :setvars

:usage_no
echo Automatic usage reporting disabled.
echo.
<nul (SET /P FOO=false) >%ACCEPTED_USAGE_NAME%
GOTO :setvars

:no64bit
echo RTG requires 64 bit version of windows

:no
IF DEFINED PAUSE_ON_CLOSE pause
exit /b 1

:setvars

IF "%RTG_TALKBACK%" == "" (
    IF EXIST %ACCEPTED_TALKBACK_NAME% SET /P TALKBACK_ANSWER=<%ACCEPTED_TALKBACK_NAME%
) ELSE (
    SET TALKBACK_ANSWER=%RTG_TALKBACK%
)
SET RTG_TALKBACK=-Dtalkback=%TALKBACK_ANSWER%

IF "%RTG_USAGE%" == "" (
    IF EXIST %ACCEPTED_USAGE_NAME% SET /P USAGE_ANSWER=<%ACCEPTED_USAGE_NAME%
) ELSE (
    SET USAGE_ANSWER=%RTG_USAGE%
)
SET RTG_USAGE=-Dusage=%USAGE_ANSWER%
IF NOT "%RTG_DEFAULT_THREADS%" == "" SET RTG_DEFAULT_THREADS=-Druntime.defaultThreads=%RTG_DEFAULT_THREADS%
IF NOT "%RTG_USAGE_HOST%" == "" SET RTG_USAGE=%RTG_USAGE% -Dusage.host=%RTG_USAGE_HOST%
IF NOT "%RTG_USAGE_DIR%" == "" SET RTG_USAGE=%RTG_USAGE% "-Dusage.dir=%RTG_USAGE_DIR%"
IF "%RTG_USAGE_OPTIONAL%" == "" (GOTO :skipoptional)
IF NOT "x%RTG_USAGE_OPTIONAL:username=x%" == "x%RTG_USAGE_OPTIONAL%" SET RTG_USAGE=%RTG_USAGE% -Dusage.log.username=true
IF NOT "x%RTG_USAGE_OPTIONAL:hostname=x%" == "x%RTG_USAGE_OPTIONAL%" SET RTG_USAGE=%RTG_USAGE% -Dusage.log.hostname=true
IF NOT "x%RTG_USAGE_OPTIONAL:commandline=x%" == "x%RTG_USAGE_OPTIONAL%" SET RTG_USAGE=%RTG_USAGE% -Dusage.log.commandline=true
:skipoptional

IF "%RTG_REFERENCES_DIR%" == "" (
    SET RTG_REFERENCES_DIR="-Dreferences.dir=%~dp0%references"
) ELSE (
    SET RTG_REFERENCES_DIR="-Dreferences.dir=%RTG_REFERENCES_DIR%"
)

IF "%RTG_MODELS_DIR%" == "" (
    SET RTG_MODELS_DIR="-Dmodels.dir=%~dp0%models"
) ELSE (
    SET RTG_MODELS_DIR="-Dmodels.dir=%RTG_MODELS_DIR%"
)

IF NOT "%PROXY_HOST%" == "" SET PROXY_HOST=-Dproxy.host=%PROXY_HOST%
IF NOT "%PROXY_PORT%" == "" SET PROXY_PORT=-Dproxy.port=%PROXY_PORT%

REM set memory
IF "%RTG_MEM%" == "" (
    FOR /F "usebackq" %%A in (`CALL %RTG_JAVA% -cp %RTG_JAR% com.rtg.util.ChooseMemory %RTG_MEM_PCT%`) DO SET RTG_MEM=%%A
)

%RTG_JAVA% -Xmx%RTG_MEM% %RTG_JAVA_OPTS% %RTG_REFERENCES_DIR% %RTG_MODELS_DIR% %RTG_USAGE% %RTG_TALKBACK% %PROXY_HOST% %PROXY_PORT% %RTG_DEFAULT_THREADS% -jar %RTG_JAR% %*

:exit_zero
exit /b 0
