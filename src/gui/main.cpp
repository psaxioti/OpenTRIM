#include "mainui.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    MainUI w;
    w.show();
    return app.exec();
}
