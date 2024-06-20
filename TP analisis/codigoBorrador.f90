PROGRAM codigoBorrador
    implicit none

    integer, parameter:: N=3!grado coeficiente principal
    integer:: i
    real (8),dimension(0:N)::e,e_ant,p
    real(8),dimension(1:N)::q,q_ant
    real(8):: tol=0.001
    logical:: cond
    complex,dimension(1:N)::raices
    !asumo que solo nos interesan los ant y los actuales.

    !q(1) valores de la iteracion actual del coeficiente 1 (x^1)
    !q_ant(1) valores de la iteracion anterior del coeficiente 1 (x^1)

    !INICIALIZACION DE LOS COEFICIENTES DEL POLINOMIO.
    !p(0)=50
    !p(1)=10
    !p(2)=3
    !p(3)=4
    !p(4)=1

    p(3)=1
    p(2)=2
    p(1)=-5
    p(0)=6

    !ANALISIS SI EL POLINOMIO ESTA COMPLETO
    i=0
    do while (cond .AND. i<=N)
        cond=p(i).EQ.0
        i=i+1
    end do

    if(.NOT.cond) then
        !llama al metodo de traslacion
        !traslacion(p,N)
    end if

    open(unit=3,file='Raices.txt',status='replace')


    !CALL Bairstow(p,0,N)
    CALL QD(q,e,q_ant,e_ant,p,N,tol,raices)
    CALL imprimeRaices(raices,N)

    CONTAINS


    SUBROUTINE imprimeRaices(raices,N)
        integer::N,i
        complex,dimension(1:N)::raices
        write(*,'(A)')"PARTE REAL    PARTE IMAGINARIA"
        do i=1,N
            write(*,'(2F10.4)')REAL(raices(i)),IMAG(raices(i))
        end do
    END SUBROUTINE

    SUBROUTINE traslacion(p,N)

        real(8),dimension(0:N)::p
        integer::N,a,i

        !VALOR A PARA LA TRASLACION
        a=1
        !si a<0, se traslada hacia la derecha. Si a>0, se traslada hacia la izq

        !sustituir x por z=x-a
        do i=0,N
            !p(i)=p(i)*(())
        end do


    end SUBROUTINE

    SUBROUTINE QD(q,e,q_ant,e_ant,p,N,tol,raices)
        integer,INTENT(IN):: N
        real (8),dimension(0:N)::e,e_ant,p
        real(8),dimension(1:N)::q,q_ant
        integer:: iter,i=0,max_iter
        real (8)::error,tol,u,v
        complex:: raiz,raiz2
        complex,dimension(1:N)::raices

        !ITERACION INICIAL PARA Q
        do i=1,N
            if (i .EQ. 1) then
                 q(i)=-(p(N-1)/p(N))
             else
                 q(i)=0.0
             end if
             write (3,'(A,F10.5)',ADVANCE='NO')"          ",q(i)
             !dejo un espacio de 10 digitos para mantener la forma de la matriz
        end do
        write(3,'(A)')"          "!espacio en blanco y bajada de linea

        !ITERACION INICIAL PARA E
        do i=0,N
            !Si pertenece a las ultimas columnas
            if (i==N .OR. i==0) then
                e(i)=0.0
            else
                e(i)=p(n-i-1)/p(n-i)
            end if
            if (i .NE. N) then
                write (3,'(F10.5,A)',ADVANCE='NO')e(i),"          "
            else
                write(3,'(F10.5)')e(i)!valor de e_ant(n) y bajada de linea
            end if
        end do

        iter=1
        error=999!para que entre en la primera iteracion
        !SIGUIENTES ITERACIONES
        max_iter=1000
        do while (error>=tol .AND. iter<max_iter)
            !la condicion de corte tiene que ver con el error entre q y q_ant? que tolerancia deberia considerar?
            q_ant=q
            e_ant=e
            !ITERACION ACTUAL PARA Q
            do i=1,N
                q(i)=e_ant(i)-e_ant(i-1)+q_ant(i)
                write (3,'(A,F10.5)',ADVANCE='NO')"          ",q(i)
            end do
                write(3,'(A)')"          "
        
            !ITERACION ACTUAL PARA E
            do i=0,N
                if(i==0 .OR. i==N) then
                    e(i)=0.0
                else
                    e(i)=(q(i+1)/q(i))*e_ant(i)
                    !DUDA: q(0) no existe, lo inicializo igual?
                end if
                if (i .NE. N) then
                    write (3,'(F10.5,A)',ADVANCE='NO')e(i),"          "
                else
                    write(3,'(F10.5)')e(i)!valor de e(n) y bajada de linea
                end if
            end do
            error=maxval(abs(e))
            !error=calculaError(q,q_ant,n)!saca el promedio del error
            iter=iter+1
        end do

        do i=0,N-1
            if (abs(e(i))>tol) then
                write(*,*)"SE CUMPLIO RAIZ IMAGINARIA"
                call Bairstow(p,u,v,N)
                call Resolvente(-u,-v,raiz,raiz2)
                raices(i)=raiz
                raices(i+1)=raiz2
                write(*,'(4F10.4)')REAL(raiz),IMAG(raiz),REAL(raiz2),IMAG(raiz2)
            else
                write(*,*)"SE CUMPLIO RAIZ REAL"
                write(*,'(F10.4)')q(i+1)
                raices(i+1)=DCMPLX(q(i+1),0)!es una raiz real, tendra parte imaginaria nula.
                !la funcion DCMPLX crea un numero complejo a partir de 2 parametros: la parte real del complejo y la parte imaginaria del complejo.
            end if
        end do


    END SUBROUTINE QD


SUBROUTINE Resolvente(u,v,raiz,raiz2)
    real(8)::u,v
    complex::raiz,raiz2,aux
    !x^2-u*x-v
    !-b +- sqrt(b^2-4ac)/2a
    write(*,'(A,2F10.4)')"u y v dentro de la resolvente",u,v
    aux=u**2.0-4.0*v
    !write(*,'(A,F10.4)')"resolvente",aux
    raiz=(-u+sqrt(aux))/2.0
    raiz2=(-u-sqrt(aux))/2.0

END SUBROUTINE

    !la idea es en un ciclo de repetion, como un do, ir revisando que coeficientes oscilan en su valor de e(i) respecto de e_ant(i), y con ese indice i (para el coeficiente p(i)) llamar a la subrutina.
    !paso tambien p porque necesito los coeficientes p(n) y p(n-1).
    SUBROUTINE Bairstow(pol,u,v, N)
        integer,intent(IN)::N
        real(8),dimension(0:N),intent(IN):: pol!polinomio
        !real(8),dimension(-2:N)::q,p
        integer:: t
        real(8)::tol,error,h,k,q,p,q_ant1,q_ant2,p_ant1,p_ant2,p_ant3
        real(8),intent(INOUT)::u,v


        !TOLERANCIA
        tol=0.001
        error=10.0*tol

        !PASO INICIAL
        u=-pol(N-1)/pol(N)!u0
        v=-pol(N-1)/pol(N)!v0

        do while(error>=tol)

            !inicializacion
            q_ant1=0.0
            q_ant2=0.0
            p_ant1=0.0
            p_ant2=0.0

            do t=N,1,-1
                q = pol(t) + u * q_ant1 + v * q_ant2
                q_ant2 = q_ant1
                q_ant1 = q
                p = q + u * p_ant1 + v * p_ant2
                p_ant3 = p_ant2
                p_ant2 = p_ant1
                p_ant1 = p
            end do
            q = pol(0) + u * q_ant1 + v * q_ant2!ultima iteracion de q


            !CALCULO DE H Y K
            h = (q * p_ant3 - q_ant1 * p_ant2)/(p_ant2**2.0 - p_ant1 * p_ant3)
            k = (q_ant1 * p_ant1 - q * p_ant2) / (p_ant2**2.0 - p_ant1 * p_ant3)
             
            !ACTUALIZACION DE U Y V
            u=u+h
            v=v+k

            !CALCULO DEL NUEVO ERROR
            if(abs(q)>abs(q_ant1))then
                error=abs(q)
            else
                error=abs(q_ant1)
            end if
            write(*,'(A,F10.4)')"error",error
        end do
            
        !INICIALIZACION Q Y P
        !q(-2)=0
        !q(-1)=0
        !p(-2)=0
        !p(-1)=0

        !PASO INICIAL
        !u=-pol(N-1)/pol(N)!u0
        !v=-pol(N-1)/pol(N)!v0
        
        !do while(error>tol)

         !   do t=N,1,-1
         !       q(t)=pol(t)+(u*q(t-1))+(v*q(t-2))
         !       p(t)=q(t)+(u*p(t-1))+(v*p(t-2))
         !   end do
         !   q(0)=pol(0)+(u*q(-1))+(v*q(-2))!p deja de iterar en N-1

            !CALCULO DE H Y K
         !   h=( ( q(N) * p(N-3) ) - ( q(N-1) * p(N-2) ) )/( (p(N-2)**2) - ( p(N-1) * p(N-3) ) )!requiere verificacion de que N>0
        !    k=((q(N-1)*p(N-1))-(q(N)*p(N-2)))/((p(N-2)**2)-(p(N-1)*p(N-3)))

            !ACTUALIZACION DE U Y V
        !    u=u+h
        !    v=v+k
        !    write(*,'(2F10.4)')u,v
            !CALCULO DEL NUEVO ERROR
        !    error=max(abs(q(N)),abs(q(N-1)))
            !error=abs(q(N)-q(N-1))
       ! end do



    END SUBROUTINE Bairstow

    FUNCTION calculaError(q,q_ant,N)
        real(8)::calculaError
        integer::N,i
        real(8),dimension(1:N)::q,q_ant

        calculaError=0
        do i=1,N
            calculaError=calculaError+abs(q(i)-q_ant(i))
        end do

        calculaError=calculaError/N

    END FUNCTION

    
END PROGRAM