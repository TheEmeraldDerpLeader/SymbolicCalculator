## Makefile in the works, not working and tested yet. Should also add -M stuff for when headers change


Default C++11 settings should work, may have to remove some warnings. See Makefile.mkf for a reference on what settings to use.

# SFML Dependency
A small, disjointed portion of this project depends on SFML, and currently the main file (Program.cpp) has a test that utilizes this. However, it is easy to remove this dependency: remove all references to SFMLRendererAPI in Program.cpp, then remove SFMLRenderer.cpp from the build. SFMLRendererAPI.hpp provides the few method definitions needed for this test while separating the dependency to SFML from the rest of the project. If you would like to build this project with SFMLRendererAPI.cpp, see https://www.sfml-dev.org/tutorials/2.6/ under "Getting Started"