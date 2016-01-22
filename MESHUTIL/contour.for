C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C    PROGRAM contour                                                   C
C                                                                      C
C    Reads a SOLIDS_ISO results file and wrirtes .msh files ready to  C
C    be visualized by Gmesh                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      PROGRAM contour
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*10 FILENAME, FGAR
C
C
      ALLOCATABLE U(:),V(:)
C
      WRITE(*,*) 'INPUT THE JOB NAME(max 10 characters):'
      READ(*,*) FILENAME
      LST=LEN_TRIM(FILENAME)
      WRITE(*,*) 'INPUT THE NUMBER OF NODAL POINTS:'
      READ(*,*) NMNP
C
      INOD =3
      IOUT1=4
      IOUT2=5
      IOUT3=6
C
      OPEN(UNIT=INOD,FILE="out.txt",FORM='FORMATTED')
      OPEN(UNIT=IOUT1,FILE =FILENAME(1:LST)//"H.msh",position='append')
      OPEN(UNIT=IOUT2,FILE =FILENAME(1:LST)//"V.msh",position='append')
      OPEN(UNIT=IOUT3,FILE =FILENAME(1:LST)//"F.msh",position='append')
C
      ALLOCATE(U(NMNP),V(NMNP))
C
      DO I=1,NMNP
         READ(INOD,*) U(I), V(I)
      END DO
C
      WRITE(IOUT1,*)
      WRITE(IOUT1,1000) '$NodeData'
      WRITE(IOUT1,*) 1
      WRITE(IOUT1,*) '"Horizontal Displacement field"'
      WRITE(IOUT1,*) 1
      WRITE(IOUT1,*) 0
      WRITE(IOUT1,*) 3
      WRITE(IOUT1,*) 0
      WRITE(IOUT1,*) 1
      WRITE(IOUT1,*) NMNP
      do i =1,NMNP
         WRITE(IOUT1,*) i, U(I)
      end do
      WRITE(IOUT1,*) '$EndNodeData'
C
      REWIND(INOD)
      WRITE(IOUT2,*)
      WRITE(IOUT2,1000) '$NodeData'
      WRITE(IOUT2,*) 1
      WRITE(IOUT2,*) '"Vertical Displacement field"'
      WRITE(IOUT2,*) 1
      WRITE(IOUT2,*) 0
      WRITE(IOUT2,*) 3
      WRITE(IOUT2,*) 0
      WRITE(IOUT2,*) 1
      WRITE(IOUT2,*) NMNP

      do i =1,NMNP
         WRITE(IOUT2,*) I, V(I)
      end do
      WRITE(IOUT2,*) '$EndNodeData'
C
      REWIND(INOD)
      WRITE(IOUT3,*)
      WRITE(IOUT3,1000) '$NodeData'
      WRITE(IOUT3,*) 1
      WRITE(IOUT3,*) '"Displacement amplitude"'
      WRITE(IOUT3,*) 1
      WRITE(IOUT3,*) 0
      WRITE(IOUT3,*) 3
      WRITE(IOUT3,*) 0
      WRITE(IOUT3,*) 3
      WRITE(IOUT3,*) NMNP

      do i =1,NMNP
         WRITE(IOUT3,*) I, U(I), V(I), 0.0
      end do
      WRITE(IOUT3,*) '$EndNodeData'
C
 1000 FORMAT (1A9)
C
      END PROGRAM