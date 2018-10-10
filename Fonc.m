function f=fun_v(q)
   deff([v]=vecf(v0),v=v0/norm(v0)^3);
   soleil=1:3;jupiter=4:6;saturne=7:9;
   f(soleil)     =-G*m0*m1*vecf(q(soleil) -q(jupiter))-G*m0*m2*vecf(q(soleil) -q(saturne));
   f(jupiter)    =-G*m1*m0*vecf(q(jupiter)-q(soleil))-G*m1*m2*vecf(q(jupiter)-q(saturne));
   f(saturne)    =-G*m2*m0*vecf(q(saturne)-q(soleil))-G*m2*m1*vecf(q(saturne)-q(jupiter));
   f(uranus)     =-G*m3*m0*vecf(q(uranus)-q(soleil))-G*m2*m1*vecf(q(uranus)-q(jupiter));
   f(neptune)    =-G*m4*m0*vecf(q(neptune)-q(soleil))-G*m2*m1*vecf(q(neptune)-q(jupiter));
   f(pluton)     =-G*m5*m0*vecf(q(pluton)-q(soleil))-G*m2*m1*vecf(q(pluton)-q(jupiter));
function

function f=fun_u(p)
   f=[p(1:3)/m0;p(4:6)/m1;p(7:9)/m2];
end
function [vp,vq]=euler_symplectique(n,h,p,q)
vp=p;vq=q;
for i=1:n
   q=q+h*fun_u(p);   
   p=p+h*fun_v(q);
   vq=[vq,q];           
   vp=[vp,q];
end