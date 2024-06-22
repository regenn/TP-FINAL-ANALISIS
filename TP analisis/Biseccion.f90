PROGRAM BISECCION

  IMPLICIT NONE

! DECLARACION DE VARIABLES

real(8) a,b,cota
integer(4) maxiter

! SECCION EJECUTABLE

call grafico
print *,'Con los datos del gráfico'
print *, 'Ingrese límite inferior y superior del intervalo donde cree que se encuentra la raíz'
read *, a,b
print *,'Ingrese cota máxima de error'
read *,cota
!print *,'Ingrese numero maximo de iteraciones'
!read *,maxiter
If (func(a)*func(b).lt.0) then                           !TEOREMA DE BOLZANO lt=menor que  
   call algoritmobisec(a,b,cota,maxiter)
  else
   print *,'Hay un nro. par de raices o no existe raiz'
end if    

CONTAINS

! SUBRUTINAS Y FUNCIONES INTERNAS
FUNCTION func(x)
 
 real(8) x,func
 func=2-cosh(x)-sin(2*x**2) 
 !func=log(x**2.+1.)-exp(x/2.)*cos(3.1415926536*x)
                                   !cambiarloooooooooooooooooooooooooooooooooooooooo
 
END FUNCTION func 

SUBROUTINE intervalo(a,b)

  real(8) a,b,xmedio,valora,valorx
  
  xmedio=(a+b)/2.0
  valora=func(a)
  valorx=func(xmedio)
  If ((valora*valorx).lt.0) then                          !TEOREMA DE BOLZANO
     b=xmedio
   else
     a=xmedio
  end if
  
END SUBROUTINE intervalo

SUBROUTINE algoritmobisec(a,b,cota,maxiter)

 real(8) a,b,cota,valor,m
 integer(4) iter,maxiter
   
 OPEN (UNIT=2,FILE='Biseccion.dat', STATUS='REPLACE')
 write(2,'(3F10.5)',advance='no') a,b,(a+b)/2
 write(2,*)
  iter=1
 m=(a+b)/2.                                          !PUNTO MEDIO DEL INTERVALO
 Do while (abs(func(a)-func(m)).gt.cota)    !.and.(iter.le.maxiter)) seria lo mismo hacer abs(b-m)   BORRAR SEGÚN CUAL ME ESTÉ PIDIENDO
   write(2,'(2F10.5)',advance='no') a,b
   write(2,*)
   call intervalo(a,b)
   m=(a+b)/2.
   iter=iter+1
 end do
 write(2,'(2F30.15)',advance='no') a,b
 print *,'La raiz mas aproximada es ',m               !ES LO MISMO QUE SEA CUALQUIERA DE LOS DOS LÍMITES, POR ESO DESIGNO M
 valor=abs(func(m))
 print *,'El error en y es'
 print *, valor                                       
 print*, 'Se halló en',iter,'iteraciones'
 CLOSE(2,STATUS='KEEP')

END SUBROUTINE algoritmobisec 

SUBROUTINE grafico
 
 real(8) xi,xf,h,i

OPEN(2,FILE='Grafbisec.dat',STATUS='replace')
print *,'Ingrese x inicial y final para graficar'
read *, xi,xf
print *,'Indique a que paso desea graficar'
read *,h
i=xi
Do while ((i.ge.xi).and.(i.le.xf))
   write(2,'(2F10.5)',advance='no') i,func(i)
   write(2,*)
   i=i+h
end do
CLOSE (2,STATUS='KEEP')
print *,'EL grafico de la funcion es '
CALL SYSTEM ("gnuplot -persist myScript.p")

END SUBROUTINE grafico
          
END PROGRAM

