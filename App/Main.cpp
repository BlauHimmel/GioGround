#include "GioApplication.hpp"

GioApplication * g_pApp = nullptr;

int main(int argc, char ** argv)
{
	g_pApp = new GioApplication(argc, argv);
	g_pApp->start();
	delete g_pApp;
	system("PAUSE");
    return 0;
}

