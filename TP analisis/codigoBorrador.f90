PROGRAM codigoBorrador
    implicit none

    integer,parameter:: N=3!grado coeficiente principal
    integer:: i,Naux
    real (8),dimension(0:N)::e,e_ant,p,pDerivada
    real(8),dimension(1:N)::q,q_ant
    real(8):: x,c,tol=0.000001
    logical:: cond
    complex,dimension(1:N)::raices
    !asumo que solo nos interesan los ant y los actuales.

    !q(1) valores de la iteracion actual del coeficiente 1 (x^1)
    !q_ant(1) valores de la iteracion anterior del coeficiente 1 (x^1)

    !INICIALIZACION DE LOS COEFICIENTES DEL POLINOMIO.
    !p(0)=1
    !p(1)=-32
    !p(2)=160
    !p(3)=-256
    !p(4)=128

    !p(3)=1
    !p(2)=2
    !p(1)=-5
    !p(0)=6
    
    p(3)=1
    p(2)=2
    p(1)=2
    p(0)=-2

    Naux=N

    !x=1.0
    !CALL Derivada(pDerivada,Naux)
    !write(*,*)PolinomioEnPuntoX(P,N,x)S
    !CALL ImprimePolinomio(pDerivada,Naux)

    c=1.0
    !CALL traslacionIndeterminada(p,Naux,c)
    !CALL homoteciaIndeterminada(p,Naux,c)
    CALL ImprimePolinomio(p,Naux)
    CALL QD(q,e,q_ant,e_ant,p,Naux,tol,raices)
    CALL ImprimeRaices(raices,Naux)

    !ANALISIS SI EL POLINOMIO ESTA COMPLETO
    !i=0
    !do while (cond .AND. i<=N)
        !cond=p(i).EQ.0
        !i=i+1
    !end do

    !if(.NOT.cond) then
        !llama al metodo de traslacion
        !traslacion(p,N)
    !end if


    !CALL Bairstow(p,0,N)

    !CALL QD(q,e,q_ant,e_ant,p,N,tol,raices)
    !CALL imprimeRaices(raices,N)

    CONTAINS

    SUBROUTINE ImprimePolinomio(P,N)
        integer::N,i
        real(8),dimension(0:N)::P

        do i=0,N
            write(*,'(F10.4,A,I10)')P(i),'x^',i
        end do

    END SUBROUTINE


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

    !los coeficientes del polinomio se multiplican por el factor c, las raices no se alteran. 
    subroutine homoteciaPolinomio(p,N,c)
        real(8),dimension(0:N)::p
        integer::N
        real(8)::c

        P=P*c

    end subroutine
    
    !Genera un cambio de escala por un factor abs(c). Si c<0, invierte el sentido del eje y.
    subroutine homoteciaIndeterminada(p,N,c)
        real(8),dimension(0:N)::p
        integer::N,i
        real(8)::c

        do i=1,N
            !no se considera termino independiente P(0) porque c^0=1
            P(i)=(c**i)*P(i)
        end do

    end subroutine

    !lleva las raices a su complemento y viceversa. Para una raiz x, la lleva a 1/x.
    subroutine traslacionReciproca(p,N)
        real(8),dimension(0:N)::p,paux
        integer::N,i

        Paux=P
        do i=0,N
            !a_(n-i)*x^i=a_i*x^(n-i)
            P(N-i)=Paux(i)
        end do

    end subroutine

    subroutine traslacionIndeterminada(p,N,c)
        real(8),dimension(0:N)::p,auxDerivada
        integer::N,i,gradoDerivada
        real(8)::c

        auxDerivada=p
        gradoDerivada=N

        do i=0,N
            P(i)=PolinomioEnPuntoX(auxDerivada,gradoDerivada,c)
            P(i)=P(i)/Factorial(i)
            call Derivada(auxDerivada,gradoDerivada)
        end do
    end subroutine

    function Factorial(valor)
        real(8)::Factorial
        integer::i,valor

        Factorial=1
        if (valor.NE.0 .AND. valor.NE.1) then
            do i=valor,1,-1
                Factorial=Factorial*i
            end do
        end if

    end function


    function PolinomioEnPuntoX(P,N,x)
        real(8),dimension(0:N)::p
        integer::N,i
        real(8)::x,PolinomioEnPuntoX

        PolinomioEnPuntoX=P(0)

        do i=1,N
            PolinomioEnPuntoX=PolinomioEnPuntoX + P(i)*(x**i)
        end do

    end function

    subroutine Derivada(P,N)
        integer::N,i,aux
        real(8),dimension(0:N)::P
        !el grado del polinomio producto de la derivada se va a reducir
        !ejemplo resultado del ciclo:
        !polinomio original: 3x^2 + 5x + 6
        !su derivada: 6x + 5
        do i=0, N-1
            P(i)=P(i+1)*(i+1)
        end do

        N=N-1
        
    end subroutine

    SUBROUTINE QD(q,e,q_ant,e_ant,p,N,tol,raices)
        integer,INTENT(IN):: N
        real (8),dimension(0:N)::e,e_ant,p
        real(8),dimension(1:N)::q,q_ant
        integer:: iter,i=0,max_iter
        real (8)::error,tol,u,v
        complex:: raiz,raiz2
        complex,dimension(1:N)::raices

        open(unit=3,file="Raices.txt",status='replace')

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
        max_iter=200
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
        real(8)::u,v


        !TOLERANCIA
        tol=0.00001
        error=5.0*tol

        !PASO INICIAL
        write(*,'(A,F10.4)')"pol(N-1)",pol(N-1)
        write(*,'(A,F10.4)')"pol(N)",pol(N)
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
                !write(*,'(A,F10.4)')"pol(t)",pol(t)
                write(*,'(A,F10.4)')"u",u
                write(*,'(A,F10.4)')"v",v
                write(*,'(A,F10.4)')"q",q
                write(*,'(A,F10.4)')"q_ant1",q_ant1
            end do
            q = pol(0) + u * q_ant1 + v * q_ant2!ultima iteracion de q


            !CALCULO DE H Y K
            write(*,'(A,F10.4)')"p_ant2",p_ant1
            write(*,'(A,F10.4)')"p_ant1",p_ant2
            write(*,'(A,F10.4)')"p_ant3",p_ant3
            h = (q * p_ant3 - q_ant1 * p_ant2)/((p_ant2**2.0) - (p_ant1 * p_ant3))
            k = (q_ant1 * p_ant1 - q * p_ant2) / ((p_ant2**2.0) - (p_ant1 * p_ant3))
            write(*,'(A,F10.4)')"h",h
            write(*,'(A,F10.4)')"k",k
             
            !ACTUALIZACION DE U Y V
            u=u+h
            v=v+k

            !CALCULO DEL NUEVO ERROR
            if(abs(q)>abs(q_ant1))then
                !q y q_ant1 deben ser menores a la cota para salir, tomo el mayor
                error=abs(q)
            else
                error=abs(q_ant1)
            end if
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