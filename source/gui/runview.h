#ifndef RUNVIEW_H
#define RUNVIEW_H

#include <QWidget>

class QLineEdit;
class QLabel;
class QSpinBox;
class QProgressBar;
class QToolButton;
class QLabel;
class QGroupBox;
class QTimer;

class MyDataWidgetMapper;
class MainUI;

class RunView : public QWidget
{
    Q_OBJECT

public:

    explicit RunView(MainUI* iui, QWidget *parent = nullptr);

public slots:
    void revert();
    void onTallyUpdate();
    void onDriverStatusChanged();

signals:

private:

    MainUI* ionsui;
    MyDataWidgetMapper* mapper;
    QLineEdit* simTitle;
    QGroupBox* box3;

    std::vector<QLineEdit*> simTotals;

};



#endif // RUNVIEW_H
