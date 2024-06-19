PROGRAM codigoBorrador
    implicit none

    integer, parameter:: N=4!grado coeficiente principal
    real (8),dimension(0:N)::e,e_ant,p
    real(8),dimension(1:N)::q,q_ant
    real(8):: tol=0.005
    !asumo que solo nos interesan los ant y los actuales.

    !q(1) valores de la iteracion actual del coeficiente 1 (x^1)
    !q_ant(1) valores de la iteracion anterior del coeficiente 1 (x^1)

    !INICIALIZACION DE LOS COEFICIENTES DEL POLINOMIO.
    p(0)=1
    p(1)=-32
    p(2)=160
    p(3)=-256
    p(4)=128

    open(unit=3,file='Raices.txt',status='replace')

    CALL QD(q,e,q_ant,e_ant,p,N,tol)

    CONTAINS

    SUBROUTINE QD(q,e,q_ant,e_ant,p,N,tol)
        integer,INTENT(IN):: N
        real (8),dimension(0:N)::e,e_ant,p
        real(8),dimension(1:N)::q,q_ant
        integer:: iter,i=0
        real (8)::error,tol

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
        do while (error>=tol)
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
            error=calculaError(q,q_ant,n)
        end do


    END SUBROUTINE QD

    !la idea es en un ciclo de repetion, como un do, ir revisando que coeficientes oscilan en su valor de e(i) respecto de e_ant(i), y con ese indice i (para el coeficiente p(i)) llamar a la subrutina.
    !paso tambien p porque necesito los coeficientes p(n) y p(n-1).
    SUBROUTINE Bairstow(pol, i, N)
        integer,intent(IN)::i,N
        real(8),dimension(0:N),intent(IN):: pol!polinomio
        real(8),dimension(-2:N)::q,p
        integer:: t
        real(8)::u,v,tol,error,h,k


        !TOLERANCIA
        tol=0.005

        !INICIALIZACION Q Y P
        q(-2)=0
        q(-1)=0
        p(-2)=0
        p(-1)=0

        !PASO INICIAL
        u=-(pol(N-1)/pol(N))!u0
        v=-(pol(N-1)/pol(N))!v0
        
        do while(error>tol)

            do t=0,N!nose si es el mismo N
                q(t)=pol(t)+(u*q(t-1))+(v*q(t-2))
                p(t)=q(t)+(u*p(t-1))+(v*p(t-2))
            end do

            !CALCULO DE H Y K
            h=((q(N)*p(N-3))-(q(N-1)*p(N-2)))/((p(N-2)**2)-(p(N-1)*p(N-3)))!requiere verificacion de que N>0
            k=((q(N-1)*p(N-1))-(q(N)*p(N-2)))/((p(N-2)**2)-(p(N-1)*p(N-3)))

            !ACTUALIZACION DE U Y V
            u=u+h
            v=v+k

            !CALCULO DEL NUEVO ERROR
            error=abs(q(N)-q(N-1))
        end do

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