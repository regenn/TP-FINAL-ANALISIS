module utils
    contains


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

end module utils