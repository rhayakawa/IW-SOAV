function h=proxfw(z,gamma,w)
% w: w^{+}
    d=2*w-1;
    h=z+gamma;
    index1=(z>-1-gamma);
    h(index1)=-1;
    index2=(z>-1-gamma.*d);
    h(index2)=z(index2)+gamma(index2).*d(index2);
    index3=(z>1-gamma.*d);
    h(index3)=1;
    index4=(z>1+gamma);
    h(index4)=z(index4)-gamma(index4);

end