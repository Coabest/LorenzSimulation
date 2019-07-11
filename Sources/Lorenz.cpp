#include "../Headers/Lorenz.hpp"

Lorenz::Lorenz(){
	//Posiciones iniciales de los puntos 1, 2, 3 y 4 (son la misma)
	x1 = x2 = x3 = x4 = 0.1f;
	y1 = y2 = y3 = y4 = 0.0f;
	z1 = z2 = z3 = z4 = 0.0f;
	
	//Parametros iniciales de las constantes
  sigma = 10, rho = 28, beta = 8.f/3.f;
  
  //Escojo un salto de tiempo arbitrario
  dt = 0.005;
  
  point1xy.setRadius(0.17f); point1xy.setFillColor(Color(255,255,0,250));
  point2xy.setRadius(0.17f); point2xy.setFillColor(Color(255,0,255,200));
  point3xy.setRadius(0.17f); point3xy.setFillColor(Color(255,0,0,170));
  point4xy.setRadius(0.17f); point4xy.setFillColor(Color(0,255,255,150));
  
  point1xz.setRadius(0.17f); point1xz.setFillColor(Color(255,255,0,250));
  point2xz.setRadius(0.17f); point2xz.setFillColor(Color(255,0,255,200));
  point3xz.setRadius(0.17f); point3xz.setFillColor(Color(255,0,0,170));
  point4xz.setRadius(0.17f); point4xz.setFillColor(Color(0,255,255,150));
  
  point1yz.setRadius(0.17f); point1yz.setFillColor(Color(255,255,0,250));
  point2yz.setRadius(0.17f); point2yz.setFillColor(Color(255,0,255,200));
  point3yz.setRadius(0.17f); point3yz.setFillColor(Color(255,0,0,170));
  point4yz.setRadius(0.17f); point4yz.setFillColor(Color(0,255,255,150));
};

ld Lorenz::dx(const ld x, const ld y)            { return sigma*(y - x); }
ld Lorenz::dy(const ld x, const ld y, const ld z){ return x*(rho - z) - y; }
ld Lorenz::dz(const ld x, const ld y, const ld z){ return x*y - beta*z; }

void Lorenz::Euler(ld &x, ld &y, ld &z){
	ld x1 = dt*dx(x, y);
	ld y1 = dt*dy(x, y, z);
	ld z1 = dt*dz(x, y, z);
	x += x1;
	y += y1;
	z += z1;
}

void Lorenz::Runge_Kutta(ld &x, ld &y, ld &z){
	ld A1 = dt*dx(x, y);
	ld A2 = dt*dy(x, y, z);
	ld A3 = dt*dz(x, y, z);
	
	ld B1 = dt*dx(x+A1/2, y+A2/2);
	ld B2 = dt*dy(x+A1/2, y+A2/2, z+A3/2);
	ld B3 = dt*dz(x+A1/2, y+A2/2, z+A3/2);
	
	ld C1 = dt*dx(x+B1/2, y+B2/2);
	ld C2 = dt*dy(x+B1/2, y+B2/2, z+B3/2);
	ld C3 = dt*dz(x+B1/2, y+B2/2, z+B3/2);
	
	ld D1 = dt*dx(x+C1/2, y+C2/2);
	ld D2 = dt*dy(x+C1/2, y+C2/2, z+C3/2);
	ld D3 = dt*dz(x+C1/2, y+C2/2, z+C3/2);
	
	x += (A1 + 2*B1 + 2*C1 + D1)/6;
	y += (A2 + 2*B2 + 2*C2 + D2)/6;
	z += (A3 + 2*B3 + 2*C3 + D3)/6;
}

void Lorenz::Kutta_Merson(ld &x, ld &y, ld &z){
	ld A1 = dt*dx(x, y);
	ld A2 = dt*dy(x, y, z);
	ld A3 = dt*dz(x, y, z);
	
	ld B1 = dt*dx(x + A1/3, y + A2/3);
	ld B2 = dt*dy(x + A1/3, y + A2/3, z + A3/3);
 	ld B3 = dt*dy(x + A1/3, y + A2/3, z + A3/3);
 	
 	ld C1 = dt*dx(x + A1/6 + B1/6, y + A2/6 + B2/6);
 	ld C2 = dt*dy(x + A1/6 + B1/6, y + A2/6 + B2/6, z + A3/6 + B3/6);
 	ld C3 = dt*dz(x + A1/6 + B1/6, y + A2/6 + B2/6, z + A3/6 + B3/6);
 	
 	ld D1 = dt*dx(x + A1/8 + 3*C1/8, y + A2/8 + 3*C2/8);
 	ld D2 = dt*dy(x + A1/8 + 3*C1/8, y + A2/8 + 3*C2/8, z + A3/8 + 3*C3/8);
 	ld D3 = dt*dz(x + A1/8 + 3*C1/8, y + A2/8 + 3*C2/8, z + A3/8 + 3*C3/8);
 	
 	ld E1 = dt*dx(x + A1/2 - 3*C1/2 + 2*D1, y + A2/2 - 3*C2/2 + 2*D2);
 	ld E2 = dt*dy(x + A1/2 - 3*C1/2 + 2*D1, y + A2/2 - 3*C2/2 + 2*D2, z + A3/2 - 3*C3/2 + 2*D3);
 	ld E3 = dt*dz(x + A1/2 - 3*C1/2 + 2*D1, y + A2/2 - 3*C2/2 + 2*D2, z + A3/2 - 3*C3/2 + 2*D3);
 	
 	x += A1/6 + 2*D1/3 + E1/6;
 	y += A2/6 + 2*D2/3 + E2/6;
 	z += A3/6 + 2*D3/3 + E3/6;
 	
 	/*
 	//Errores del metodo
 	ld xe = x + A1/2 - 3*C1/2 + 2*D1;
 	ld ye = y + A2/2 - 3*C2/2 + 2*D2;
 	ld ze = z + A3/2 - 3*C3/2 + 2*D3;
 	*/
}

void Lorenz::Adams_Moulton(ld &x, ld &y, ld &z, int pasos){
	//Couldn't do it :C
}

//Sirve para 1, 2 y 3 pasos
void Lorenz::Adams_Bashforth(ld &x, ld &y, ld &z, int pasos){
	// 1 paso (Metodo de Euler)
	ld x1 = dt*dx(x, y);
	ld y1 = dt*dy(x, y, z);
	ld z1 = dt*dz(x, y, z);
	if (pasos == 1){ x += x1; y += y1; z += z1; return; }
	
	// 2 pasos
	ld x2 = x1 + dt*(3*dx(x1, y1)/2 - dx(x, y)/2);
	ld y2 = y1 + dt*(3*dy(x1, y1, z1)/2 - dy(x, y, z)/2);
	ld z2 = z1 + dt*(3*dz(x1, y1, z1)/2 - dz(x, y, z)/2);
	if (pasos == 2){ x += x2; y += y2; z += z2; return; }
	
	// 3 pasos
	ld x3 = x2 + dt*(23*dx(x2, y2)/12 - 4*dx(x1, y1)/3 + 5*dx(x2, y2)/12);
	ld y3 = y2 + dt*(23*dy(x2, y2, z2)/12 - 4*dy(x1, y1, z1)/3 + 5*dy(x, y, z)/12);
	ld z3 = z2 + dt*(23*dz(x2, y2, z2)/12 - 4*dz(x1, y1, z1)/3 + 5*dz(x, y, z)/12);
	if (pasos == 3){ x += x3; y += y3; z += z3; return; }
}

void Lorenz::run(){
	// Realiza 1 paso (2 en el caso de Runge-Kutta) con cada tick
	Euler(x1, y1, z1);
	Runge_Kutta(x2, y2, z2);
	Kutta_Merson(x3, y3, z3);
	Adams_Bashforth(x4, y4, z4, 2);
	
	//Los coeficientes sumados y restados acá solo son para centrar
	// la grafica en la simulacion, no son parte de los metodos.
	point1xy.setPosition(x1-55, y1+5);
	point2xy.setPosition(x2-55, y2+5);
	point3xy.setPosition(x3-55, y3+5);
	point4xy.setPosition(x4-55, y4+5);
	
	point1xz.setPosition(z1-30, x1);
	point2xz.setPosition(z2-30, x2);
	point3xz.setPosition(z3-30, x3);
	point4xz.setPosition(z4-30, x4);
	
	point1yz.setPosition(y1+55, z1-25);
	point2yz.setPosition(y2+55, z2-25);
	point3yz.setPosition(y3+55, z3-25);
	point4yz.setPosition(y4+55, z4-25);
}
