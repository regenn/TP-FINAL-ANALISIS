module NewtonModule
    use utils
    implicit none

    contains

    function newton_iter(grado, coefs, raiz0)
        implicit none
        integer, intent(in) :: grado
        integer  gradoAux
        real(8), dimension(0:grado) :: coefs, coefsAux
        real(8), intent(in) :: raiz0
        real(8) newton_iter, y0
        coefsAux = coefs 
        gradoAux = grado
        y0 = PolinomioEnPuntoX(coefs,grado,raiz0)
        CALL Derivada(coefsAux,gradoAux)
        newton_iter = raiz0 - y0 / PolinomioEnPuntoX(coefsAux,gradoAux,raiz0)
    end function newton_iter

    function newton(grado, coefs, raiz0, tol, maxIter)
        implicit none
        !Write hello world
        integer, intent(in) :: grado, maxIter
        real(8), dimension(0:grado) :: coefs
        real(8), intent(in) :: raiz0, tol
        real(8) newton, raiz, raizAnterior
        integer i
        i=0
        raizAnterior = raiz0 *0.5
        raiz = raiz0
        do while(i<maxIter .AND. abs(raiz - raizAnterior) > tol)
            raizAnterior = raiz
            raiz = newton_iter(grado, coefs, raiz)
            i=i+1 
        end do
        write(*,*) "Iteraciones: ", i
        newton = raiz
    end function newton


END module NewtonModule