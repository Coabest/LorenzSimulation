#include "Headers/Lorenz.hpp"

using namespace std;
using namespace sf;

int main(){
  RenderWindow window(sf::VideoMode(1800, 600), "Atractores de Lorenz :D");
  window.setFramerateLimit(120);
  View view(Vector2f(-30, 4), Vector2f(1800, 600));
  view.zoom(0.13);
  
  Lorenz lorenz;
  Font font;
  font.loadFromFile("Sources/Pixeled.ttf");
  Text Euler, RK, AM, KM, AB;
  
  Euler.setFont(font);
	Euler.setCharacterSize(5);
	Euler.setString("EULER");
	Euler.setPosition(-145, -15);
	Euler.setFillColor(Color(255,255,0,250));
	
	RK.setFont(font);
	RK.setCharacterSize(5);
	RK.setString("RUNGE-KUTTA");
	RK.setPosition(-145, -5);
	RK.setFillColor(Color(255,0,255,250));
	
	KM.setFont(font);
	KM.setCharacterSize(5);
	KM.setString("KUTTA-MERSON");
	KM.setPosition(-145, 5);
	KM.setFillColor(Color(255,0,0,250));
	
	AB.setFont(font);
	AB.setCharacterSize(5);
	AB.setString("ADAMS-BASHFORTH");
	AB.setPosition(-145, 15);
	AB.setFillColor(Color(0,255,255,250));
  
  int R = 0, G = 0, B = 0;

  while (window.isOpen()){
      Event event;
      while (window.pollEvent(event)){
          if (event.type == sf::Event::Closed)  window.close();
      
			}
			
			lorenz.run();
			
  		R += 6;
			if(R > 255){R = 0; G += 5;}
			if(G > 255){G = 0; B += 5;}
			if(B > 255)	B = 0;
			
			//limpiar la pantalla cada cierto tiempo
			if(G == 0 or G == 127) window.clear();
      window.setView(view);
      
      window.draw(Euler);
      window.draw(RK);
      window.draw(KM);
      window.draw(AB);
      
      //Point 1: Euler
      //Point 2: Runge-Kutta
      //Point 3: Kutta-Merson
      //Point 4: Adams-Bashforth
      //Comentar cualquiera para que no se muestre
      
      window.draw(lorenz.point1xy);
      window.draw(lorenz.point2xy);
      window.draw(lorenz.point3xy);
      window.draw(lorenz.point4xy);
      
      window.draw(lorenz.point1xz);
      window.draw(lorenz.point2xz);
      window.draw(lorenz.point3xz);
      window.draw(lorenz.point4xz);
      
      window.draw(lorenz.point1yz);
      window.draw(lorenz.point2yz);
      window.draw(lorenz.point3yz);
      window.draw(lorenz.point4yz);
      window.display();
  }

  return 0;
}
