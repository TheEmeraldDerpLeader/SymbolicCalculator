#include "SFMLRenderer.hpp"

#include <SFML/Graphics.hpp>

class SFMLRenderer
{
public:
    sf::RenderWindow window;
    sf::Texture sfTex;
    sf::Sprite drawRectangle;

    SFMLRenderer()
    {
        window.create(sf::VideoMode(800, 800), "Display");
        drawRectangle.setTexture(sfTex, true);
        window.setFramerateLimit(60);
    }

    int Loop()
    {
        if (window.isOpen())
        {
            sf::Event event;
            while (window.pollEvent(event))
            {
                if (event.type == sf::Event::Closed)
                    window.close();
            }

            window.clear();
            window.draw(drawRectangle);
            window.display();
            return 0;
        }
        else
            return 1;
    }

    void SetTex(unsigned char* tex, int width, int height)
    {
        sf::Image sfImg; sfImg.create(width, height, tex);
        sfTex.loadFromImage(sfImg);
        drawRectangle.setTexture(sfTex, true);
        drawRectangle.setScale(800/drawRectangle.getTextureRect().width, 800/drawRectangle.getTextureRect().height);
    }
};

SFMLRendererAPI::SFMLRendererAPI()
{
	sfr = new SFMLRenderer();
}
int SFMLRendererAPI::RunLoop()
{
    return sfr->Loop();
}
void SFMLRendererAPI::SetTex(unsigned char* tex, int width, int height)
{
    sfr->SetTex(tex, width, height);
}
SFMLRendererAPI::~SFMLRendererAPI()
{
	delete sfr;
}
