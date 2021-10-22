!     ######spl
      FUNCTION futr(w)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      IMPLICIT NONE

! input:
      REAL w

! function value:
      REAL futr

! local:
      REAL a1, a2, a3, a4, a5,                                          &
     &     x1, x2, x3, x4, x5,                                          &
     &     t1, t2, t3, t4, t5,                                          &
     &     b1, b2, b3, b4, b5,                                          &
     &     p
!_______________________________________________________________________

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FUTR',0,ZHOOK_HANDLE)
      a1 = -10.91
      a2 = - 0.86
      a3 = - 8.60
      a4 = - 9.36
      a5 = -13.15

      x1 = 270.
      x2 = 302.
      x3 = 334.
      x4 = 367.
      x5 = 400.

      t1 = (w-x2)*(w-x3)*(w-x4)*(w-x5)
      t2 = (w-x1)*(w-x3)*(w-x4)*(w-x5)
      t3 = (w-x1)*(w-x2)*(w-x4)*(w-x5)
      t4 = (w-x1)*(w-x2)*(w-x3)*(w-x5)
      t5 = (w-x1)*(w-x2)*(w-x3)*(w-x4)

      b1 = (x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)
      b2 = (x2-x1)*(x2-x3)*(x2-x4)*(x2-x5)
      b3 = (x3-x1)*(x3-x2)*(x3-x4)*(x3-x5)
      b4 = (x4-x1)*(x4-x2)*(x4-x3)*(x4-x5)
      b5 = (x5-x1)*(x5-x2)*(x5-x3)*(x5-x4)

      p = a1*t1/b1 + a2*t2/b2 + a3*t3/b3 + a4*t4/b4 + a5*t5/b5

      futr  = EXP(p)
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('FUTR',1,ZHOOK_HANDLE)
      RETURN
      ENDFUNCTION FUTR
