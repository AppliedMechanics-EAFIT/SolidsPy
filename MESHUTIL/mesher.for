C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C    PROGRAM mesher                                                    C
C                                                                      C
C    Reads a .mesh file created with GMESH and writes a new mesh to    C
C    be used within the PYTHON code SOLIDS_ISO                         C
C                                                                      C   
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      PROGRAM mesher
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*10 FILENAME, FGAR
C
      IOUT=4
      INOD=3
      IELE=2
      IIN=5
      WRITE(*,*) 'INPUT THE JOB NAME(max 10 characters):'
      READ(*,*) FILENAME
      WRITE(*,*) 'INPUT THE ELEMENT TYPE(1-Q; 2 2-tri; 3 1-tri):'
      READ(*,*)  IE
      LST=LEN_TRIM(FILENAME)
      OPEN(UNIT=IIN,FILE =FILENAME(1:LST)//".msh",FORM='FORMATTED')
      OPEN(UNIT=INOD,FILE="nodes.txt",FORM='FORMATTED')
      OPEN(UNIT=IELE,FILE="eles.txt",FORM='FORMATTED')
C
      CALL NODSECT(IIN,INOD)
C
      READ(IIN, *) FGAR
      READ(IIN, *) FGAR
      READ(IIN, *) NUMEL
C
      SELECT CASE(IE)
C
          CASE(1)
          CALL LINQUAD(IIN,IELE,NUMEL,IE)

          CASE(2)
          CALL QADTRIAN(IIN,IELE,NUMEL,IE)

          CASE(3)
          CALL LINTRIAN(IIN,IELE,NUMEL,IE)
C
      END SELECT
C
      END PROGRAM