MODULE decimal
  !
  ! Modulo donde se declaran los tipos de las variables a usar
  ! dp = doble  precision
  ! sp = simple precision
  !
  ! Este modulo es utilizado por todos los programas
  !
  IMPLICIT NONE
  INTEGER, PARAMETER::dp=KIND(1.D0)
  INTEGER, PARAMETER::sp=KIND(1.0)
  !
END MODULE DECIMAL
