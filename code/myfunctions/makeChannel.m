function [Hcomp,H]=makeChannel(m,n)

    Hcomp=(randn(m,n)+1i*randn(m,n))/sqrt(2);

    Hreal=real(Hcomp);
    Himag=imag(Hcomp);
    H=[Hreal -Himag;Himag Hreal];
end