PROGRAM codigoBorrador
    use utils
    use newtonModule
    implicit none

    integer,parameter:: N=4!grado coeficiente principal
    integer:: i,Naux
    real (8),dimension(0:N)::e,e_ant,p,pDerivada
    real(8),dimension(1:N)::q,q_ant
    real(8):: x,c,tol=0.00001
    logical:: cond
    character::opcion
    complex,dimension(1:N)::raices
    !asumo que solo nos interesan los ant y los actuales.

    !q(1) valores de la iteracion actual del coeficiente 1 (x^1)
    !q_ant(1) valores de la iteracion anterior del coeficiente 1 (x^1)

    !INICIALIZACION DE LOS COEFICIENTES DEL POLINOMIO.
    p(4)=-100.0
    p(3)=100.0
    p(2)=50.0
    p(1)=-30.0
    p(0)=-1.0

    Naux=N

    CALL ImprimePolinomio(p,Naux)!imprime el polinomio original en pantalla

     !ANALISIS SI EL POLINOMIO ESTA COMPLETO
    if (.NOT.completo(p,Naux)) then
        write(*,*)"El polinomio no esta completo, se debe aplicar una traslacion sobre la indeterminada."
        write(*,*)"'S'-Realizar la transformacion requerida."
        write(*,*)"'N'-Salir."
        read(*,*)opcion
        if(opcion.EQ.'S'.OR.opcion.EQ.'s') then
            write(*,*)"Ingrese el factor c para aplicar la transformacion"
            read(*,*)c
            CALL traslacionIndeterminada(P,Naux,c)
            CALL ImprimePolinomio(p,Naux)!imprime el polinomio transformado en pantalla
        end if
    end if

    if (completo(p,Naux)) then
        CALL menu(p,Naux)
        CALL QD(q,e,q_ant,e_ant,p,Naux,tol,raices)
        CALL ImprimeRaices(raices,Naux)
    end if

    c=0.5
    CONTAINS

    function completo(P,N)
        integer::N
        real(8),dimension(0:N)::P
        logical:: completo
        !ANALISIS SI EL POLINOMIO ESTA COMPLETO
        i=0
        completo = .TRUE.
        do while (completo .AND. i<=N)
            completo=p(i).NE.0
            i=i+1
        end do
        completo = i.GT.N
        !Imrpimo completo
    end function 

    SUBROUTINE menu(P,N)
        integer::N
        real(8),dimension(0:N)::P
        real(8)::c , auxNewton
        integer::opcion
        opcion = 1
        do while(opcion<=5 .AND. opcion>=1)
            write(*,*)"Elija la transformacion a realizar"
            write(*,*)"#####################################"
            write(*,*)"1: Traslacion sobre la indeterminada."
            write(*,*)"2: Traslacion reciproca."
            write(*,*)"3: Homotecia sobre un polinomio."
            write(*,*)"4: Homotecia sobre la indeterminada."
            write(*,*)"5: Graficar el polinomio."
            write(*,*)"6: No aplicar ninguna transformacion."
            write(*,*)"7: Aproximar raiz con Newton."
            write(*,*)"#####################################"
            read(*,*)opcion

            write(*,'(A,I3)')"La opcion elegida es la opcion ",opcion

            if(opcion.NE.6) then
                if (opcion.NE.2 .AND. opcion.NE.5 .AND. opcion.NE.7) then
                    write(*,*)"Ingrese el factor c para aplicar la transformacion"
                    read(*,*)c
                end if
                select case(opcion)
                case(1)
                    call TraslacionIndeterminada(P,Naux,c)
                case (2)
                    call TraslacionReciproca(P,Naux)
                case (3)
                    call homoteciaPolinomio(P,Naux,c)
                case (4)
                    call homoteciaIndeterminada(P,Naux,c)
                case (5)
                    call graficar_funcion(P,Naux)
                case (7)
                    write(*,*)"Ingrese el valor de x para aproximar la raiz"
                    read(*,*)x
                    auxNewton = newton(Naux,P,x,tol,100)
                    write(*,'(A,F10.4)')'La raiz aproximada es:',auxNewton
                case default
                    !No se eligio ninguna transformacion 
                end select

                if (opcion<=4 .AND. opcion>=1) then
                    write(*,*)"El nuevo polinomio obtenido tras la transformacion es:"
                    call imprimePolinomio(P,N)

                end if
        end if

        end do
        

    end subroutine

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

    SUBROUTINE QD(q,e,q_ant,e_ant,p,N,tol,raices)
        integer,INTENT(IN):: N
        real (8),dimension(0:N)::e,e_ant,p
        real(8),dimension(1:N)::q,q_ant
        integer:: iter,i=0,max_iter
        real (8)::error,tol,u,v
        complex(8):: raiz,raiz2
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
        write(*,'(A,I3)')'Q-D finalizado, iteraciones:',iter
        do i=0,N-1
            if (abs(e(i))>tol) then
                !u y v iniciales para Bairstow
                u = q(i) + q(i+1)
                !u=0.0
                !v=2.0
                !raices del factor cuadratico
                v = -1. * q_ant(i) * q(i+1)
                call Bairstow(p,u,v,N)
                call Resolvente(-u,-v,raiz,raiz2)
                raices(i)=raiz
                raices(i+1)=raiz2
                !write(*,'(4F10.4)')REAL(raiz),IMAG(raiz),REAL(raiz2),IMAG(raiz2)
            else
                raices(i+1)=DCMPLX(q(i+1),0)!es una raiz real, tendra parte imaginaria nula.
                !la funcion DCMPLX crea un numero complejo a partir de 2 parametros: la parte real del complejo y la parte imaginaria del complejo.
            end if
        end do

        close(3)

    END SUBROUTINE QD

    SUBROUTINE Resolvente(u, v, raiz1, raiz2)

        REAL(8) a,u,v
        COMPLEX(8) aux,raiz1,raiz2
            a = 1.0
            aux = u**2 - 4. * a * v
            raiz1 = (-u + sqrt(aux)) / (2.*a)
            raiz2 = (-u - sqrt(aux)) / (2.*a)
    END SUBROUTINE Resolvente
    

    !la idea es en un ciclo de repetion, como un do, ir revisando que coeficientes oscilan en su valor de e(i) respecto de e_ant(i), y con ese indice i (para el coeficiente p(i)) llamar a la subrutina.
    !paso tambien p porque necesito los coeficientes p(n) y p(n-1).
    SUBROUTINE Bairstow(pol,u,v, N)
        integer,intent(IN)::N
        real(8),dimension(0:N),intent(IN):: pol!polinomio
        !real(8),dimension(-2:N)::q,p
        integer:: t
        real(8)::tol,error,h,k,q,p,q_ant1,q_ant2,p_ant1,p_ant2,p_ant3
        real(8)::u,v

        write(*,'(A)')'Iniciando Teorema de Bairstow por tener raices comodulares'


        !TOLERANCIA
        tol=0.00001
        error=2.0*tol

        !PASO INICIAL
        !u=-pol(N-1)/pol(N)!u0
        !v=-pol(N-1)/pol(N)!v0

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
            end do
            q = pol(0) + u * q_ant1 + v * q_ant2!ultima iteracion de q

            !CALCULO DE H Y K
            h = (q * p_ant3 - q_ant1 * p_ant2)/((p_ant2**2.0) - (p_ant1 * p_ant3))
            k = (q_ant1 * p_ant1 - q * p_ant2) / ((p_ant2**2.0) - (p_ant1 * p_ant3))
             
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

    subroutine graficar_funcion(polinomio, N)
        REAL(8), DIMENSION(0:N) :: polinomio
        character(len=1000) :: polinomio_string
        character(len=80) :: script, N_char, term_str
        integer N, i, length
        ! Crear un archivo de script de Gnuplot
        script = 'grafico_polinomio.txt'
        
        ! Convertir el entero N a una cadena de caracteres
        write(N_char, '(I0)') N
    
        ! Inicializar el string del polinomio
        polinomio_string = ""
    
        ! Construir el string del polinomio
        do i = N, 1, -1
            write(term_str, '(F10.4,A,I2,A)') polinomio(i), ' * x**', i, ' + '
            length = len_trim(term_str)
            polinomio_string = trim(polinomio_string)//trim(term_str(1:length))
        end do
        write(term_str, '(F10.4)') polinomio(0)
        polinomio_string = trim(polinomio_string)//trim(term_str)
    
        ! Crear el archivo de script de Gnuplot
        open(unit=2, file=script, status='replace')
        write(2,*) 'set autoscale'
        write(2,*) 'set grid'
        write(2, *) 'set title "Polinomio de grado ' // trim(adjustl(N_char)) // '"'
        write(2, *) 'set xlabel "Eje X"'
        write(2, *) 'set ylabel "Eje Y"'
        write(2, *) 'set label 1 "' // trim(polinomio_string) // '" at graph 0.02, graph 0.9'

        ! Escribir la primera parte del comando plot
        write(2, '(A)', advance='NO') 'plot [-2:2] [-2:2] ' // trim(polinomio_string)

        ! Escribir la segunda parte del comando plot en una nueva línea para evitar truncamiento
        write(2, '(A)') ' with lines lc 1 title "Polinomio", 0 title "Eje X" with lines lc 3 lw 2'
    
        close(2)
        ! Llamar a Gnuplot para mostrar la gráfica
        call system('gnuplot -persist ' // script)
    
        ! Imprimir el polinomio en pantalla
        print *, "Polinomio: ", trim(polinomio_string)
    
    end subroutine graficar_funcion
    
END PROGRAM