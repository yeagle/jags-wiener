!define APP_NAME "JAGS-WIENER-MODULE"
!define PUBLISHER "JAGS-WIENER-MODULE"

;Name used for JAGS Wiener module registry keys
!define WIENER_KEYNAME "${APP_NAME}-${VERSION}"
;Name visible to users, for shortcuts, installation directories, etc.
!define WIENER_VISIBLE_NAME  "${APP_NAME} ${VERSION}"

!define MULTIUSER_MUI
!define MULTIUSER_EXECUTIONLEVEL Highest
!define MULTIUSER_INSTALLMODE_COMMANDLINE

!define INSTDIR_REG_ROOT "SHELL_CONTEXT"
!define INSTDIR_REG_KEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\${WIENER_KEYNAME}"

!define MULTIUSER_INSTALLMODE_INSTDIR "${PUBLISHER}\${WIENER_KEYNAME}"
!define MULTIUSER_INSTALLMODE_INSTDIR_REGISTRY_KEY "SOFTWARE\${PUBLISHER}\${WIENER_KEYNAME}"
!define MULTIUSER_INSTDIR_REGISTRY_VALUENAME "InstallDir"

!addincludedir ${JAGSINC}
!include AdvUninstLog.nsh
!include MultiUser64.nsh
!include "MUI2.nsh"
!include "Sections.nsh"
!include "x64.nsh"
!include LogicLib.nsh
 
Name "${WIENER_VISIBLE_NAME}"
OutFile "${APP_NAME}-${VERSION}.exe"

!define APP_REG_KEY "Software\${PUBLISHER}\${WIENER_KEYNAME}"
!define PUB_REG_KEY "Software\${PUBLISHER}"


; Installer pages

!insertmacro UNATTENDED_UNINSTALL

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE ${LICENSE}
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

; Uninstaller pages

!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_UNPAGE_FINISH

!insertmacro MUI_LANGUAGE "English"

Section #Default section

   WriteRegStr ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "InstallDir" "$INSTDIR"
   WriteRegStr ${INSTDIR_REG_ROOT} "${APP_REG_KEY}"     "InstallDir" "$INSTDIR"

   # Information for uninstaller
   WriteRegStr ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "DisplayName" "${WIENER_VISIBLE_NAME}"
   WriteRegStr ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "UninstallString" "${UNINST_EXE}"
   WriteRegStr ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "Publisher" "${PUBLISHER}"
   WriteRegStr ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "DisplayVersion" "${VERSION}"
   WriteRegStr ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "URLInfoAbout" "http://jags-wiener.sourceforge.net"
   WriteRegDWORD ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "NoModify" 1
   WriteRegDWORD ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "NoRepair" 1
   WriteRegStr ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "InstallLocation" "$INSTDIR"

   Call GetInstalledSize
   pop $0
   WriteRegDWORD ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "EstimatedSize" "$0"

SectionEnd

Section "32-bit installation" Sec32

#   SetOutPath "$INSTDIR\i386\lib"
#   !insertmacro UNINSTALL.LOG_OPEN_INSTALL
#   File inst32\lib\*.dll.a
#   File inst32\lib\*.la
#   !insertmacro UNINSTALL.LOG_CLOSE_INSTALL

   SetOutPath "$INSTDIR\i386\modules"
   !insertmacro UNINSTALL.LOG_OPEN_INSTALL
   File /r inst32\lib\JAGS\modules-${MAJOR}\*
   !insertmacro UNINSTALL.LOG_CLOSE_INSTALL

   Push @JAGS_HOME@               #text to be replaced
   Push $INSTDIR\i386             #replace with
   Push all                       #replace all occurrences
   Push all                       #replace all occurrences
   Call AdvReplaceInFile

SectionEnd #32-bit installation

Section "64-bit installation" Sec64

#   SetOutPath "$INSTDIR\x64\lib"
#   !insertmacro UNINSTALL.LOG_OPEN_INSTALL
#   File inst64\lib\*.dll.a
#   File inst64\lib\*.la
#   !insertmacro UNINSTALL.LOG_CLOSE_INSTALL

   SetOutPath "$INSTDIR\x64\modules"
   !insertmacro UNINSTALL.LOG_OPEN_INSTALL
   File /r inst64\lib\JAGS\modules-${MAJOR}\*
   !insertmacro UNINSTALL.LOG_CLOSE_INSTALL

   Push @JAGS_HOME@               #text to be replaced
   Push $INSTDIR\x64              #replace with
   Push all                       #replace all occurrences
   Push all                       #replace all occurrences
   Call AdvReplaceInFile

SectionEnd #64-bit installation

#Section "Header files" SecHeader
#
#   SetOutPath "$INSTDIR\include"
#   !insertmacro UNINSTALL.LOG_OPEN_INSTALL
#   File inst32\include\JAGS\*.h
#   File /r inst32\include\JAGS\*
#   !insertmacro UNINSTALL.LOG_CLOSE_INSTALL
#
#SectionEnd

!insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
  !insertmacro MUI_DESCRIPTION_TEXT ${Sec32} "Files for 32-bit Windows"
  !insertmacro MUI_DESCRIPTION_TEXT ${Sec64} "Files for 64-bit Windows"
  !insertmacro MUI_DESCRIPTION_TEXT ${SecHeader} "For developers who need to compile programs linked to JAGS"
!insertmacro MUI_FUNCTION_DESCRIPTION_END

Function .onInit
   !insertmacro MULTIUSER_INIT
   !insertmacro UNINSTALL.LOG_PREPARE_INSTALL
   ${If} ${RunningX64}
      ;Nothing to do
   ${Else}
      ; Deselect and hide 64-bit section
      Push $0
      SectionGetFlags ${Sec64} $0
      IntOp $0 $0 & ${SECTION_OFF}
      SectionSetFlags ${Sec64} $0
      SectionSetText  ${Sec64} ""
      Pop $0
      ; Enforce 32-bit selection
      Push $1
      SectionGetFlags ${Sec32} $1
#      IntOp $1 $1 & ${SF_SELECTED}
      IntOp $1 $1 | ${SF_RO}
      SectionSetFlags ${Sec32} $1
      Pop $1
   ${EndIf}
FunctionEnd

Function .onInstSuccess
   ;create/update log always within .onInstSuccess function
   !insertmacro UNINSTALL.LOG_UPDATE_INSTALL
FunctionEnd

Section "Uninstall"

   ;uninstall from path, must be repeated for every install logged path individually
   !insertmacro UNINSTALL.LOG_UNINSTALL "$INSTDIR\i386\modules"
   !insertmacro UNINSTALL.LOG_UNINSTALL "$INSTDIR\x64\modules"
   !insertmacro UNINSTALL.LOG_END_UNINSTALL

   DeleteRegKey /ifempty ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}"
   DeleteRegKey /ifempty ${INSTDIR_REG_ROOT} "${APP_REG_KEY}"
   DeleteRegKey /ifempty ${INSTDIR_REG_ROOT} "${PUB_REG_KEY}"

SectionEnd # end of uninstall section

Function un.onInit
   !insertmacro MULTIUSER_UNINIT
   !insertmacro UNINSTALL.LOG_BEGIN_UNINSTALL
FunctionEnd

#This is a function taken from the NSIS Wiki to replace one text string
#with another.  It doesn't work properly if there is more than one instance
#of the replacement string on a line

Function AdvReplaceInFile
         Exch $0 ;file to replace in
         Exch
         Exch $1 ;number to replace after
         Exch
         Exch 2
         Exch $2 ;replace and onwards
         Exch 2
         Exch 3
         Exch $3 ;replace with
         Exch 3
         Exch 4
         Exch $4 ;to replace
         Exch 4
         Push $5 ;minus count
         Push $6 ;universal
         Push $7 ;end string
         Push $8 ;left string
         Push $9 ;right string
         Push $R0 ;file1
         Push $R1 ;file2
         Push $R2 ;read
         Push $R3 ;universal
         Push $R4 ;count (onwards)
         Push $R5 ;count (after)
         Push $R6 ;temp file name
         GetTempFileName $R6
         FileOpen $R1 $0 r ;file to search in
         FileOpen $R0 $R6 w ;temp file
                  StrLen $R3 $4
                  StrCpy $R4 -1
                  StrCpy $R5 -1
        loop_read:
         ClearErrors
         FileRead $R1 $R2 ;read line
         IfErrors exit
         StrCpy $5 0
         StrCpy $7 $R2
 
        loop_filter:
         IntOp $5 $5 - 1
         StrCpy $6 $7 $R3 $5 ;search
         StrCmp $6 "" file_write2
         StrCmp $6 $4 0 loop_filter
 
         StrCpy $8 $7 $5 ;left part
         IntOp $6 $5 + $R3
         StrCpy $9 $7 "" $6 ;right part
         StrCpy $7 $8$3$9 ;re-join
 
         IntOp $R4 $R4 + 1
         StrCmp $2 all file_write1
         StrCmp $R4 $2 0 file_write2
         IntOp $R4 $R4 - 1
 
         IntOp $R5 $R5 + 1
         StrCmp $1 all file_write1
         StrCmp $R5 $1 0 file_write1
         IntOp $R5 $R5 - 1
         Goto file_write2
 
        file_write1:
         FileWrite $R0 $7 ;write modified line
         Goto loop_read
 
        file_write2:
         FileWrite $R0 $R2 ;write unmodified line
         Goto loop_read
 
        exit:
         FileClose $R0
         FileClose $R1
 
         SetDetailsPrint none
         Delete $0
         Rename $R6 $0
         Delete $R6
         SetDetailsPrint both
 
         Pop $R6
         Pop $R5
         Pop $R4
         Pop $R3
         Pop $R2
         Pop $R1
         Pop $R0
         Pop $9
         Pop $8
         Pop $7
         Pop $6
         Pop $5
         Pop $4
         Pop $3
         Pop $2
         Pop $1
         Pop $0
FunctionEnd

; Return on top of stack the total size of the selected (installed) sections, formated as DWORD
; Assumes no more than 256 sections are defined
Var GetInstalledSize.total
Function GetInstalledSize
	Push $0
	Push $1
	StrCpy $GetInstalledSize.total 0
	${ForEach} $1 0 256 + 1
		${if} ${SectionIsSelected} $1
			SectionGetSize $1 $0
			IntOp $GetInstalledSize.total $GetInstalledSize.total + $0
		${Endif}
 
		; Error flag is set when an out-of-bound section is referenced
		${if} ${errors}
			${break}
		${Endif}
	${Next}
 
	ClearErrors
	Pop $1
	Pop $0
	IntFmt $GetInstalledSize.total "0x%08X" $GetInstalledSize.total
	Push $GetInstalledSize.total
FunctionEnd
