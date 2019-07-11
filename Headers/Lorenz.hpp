#ifndef LORENZ_H
#define LORENZ_H
#include <bits/stdc++.h>
#include <SFML/Graphics.hpp>

using namespace std;
using namespace sf;

typedef unsigned long long ull;
typedef long long           ll;
typedef long double         ld;

class Lorenz
{
	private:
		
		ld x1, y1, z1, dt;
		ld x2, y2, z2;
		ld x3, y3, z3;
		ld x4, y4, z4;
		double sigma, rho, beta;
		
		ld dx(ld, ld);
		ld dy(ld, ld, ld);
		ld dz(ld, ld, ld);
		
		void Euler(ld&, ld&, ld&);
		void Runge_Kutta(ld&, ld&, ld&);
		void Kutta_Merson(ld&, ld&, ld&);
		void Adams_Moulton(ld&, ld&, ld&, int);
		void Adams_Bashforth(ld&, ld&, ld&, int);
		
	public:
		// Puntos dibujables.
		// para cada metodo numerico
		CircleShape point1xy, point2xy, point3xy, point4xy,
								point1xz, point2xz, point3xz, point4xz,
								point1yz, point2yz, point3yz, point4yz;
			
		Lorenz();
		void run();
};
#endif
