#ifndef SIMULATIONOPTIONSVIEW_H
#define SIMULATIONOPTIONSVIEW_H

#include <QWidget>


class QPushButton;
class QLabel;
class QTreeView;
class QTextBrowser;
class QTabWidget;

class MyDataWidgetMapper;
class MaterialsDefView;
class TargetGeometryView;
class OptionsModel;
class IonsUI;

class SimulationOptionsView : public QWidget
{
    Q_OBJECT

    Q_PROPERTY(bool modified READ modified WRITE setModified NOTIFY modifiedChanged FINAL)

public:

    QAction* whatsThisAction;
    QPushButton* helpButton;
    QLabel* ionLabel;
    QTreeView* treeView;
    QTextBrowser* jsonView;
    QTabWidget* tabWidget;
    MyDataWidgetMapper* mapper;
    MaterialsDefView* materialsView;
    TargetGeometryView* targetView;

    SimulationOptionsView(IonsUI *iui, QWidget *parent = nullptr);

    void setOptions(const QJsonDocument& json);

    bool modified() const { return modified_; }

public slots:
    void setModified(bool b = true) {
        modified_ = b;
        applyRules();
        emit modifiedChanged(b);
    }
    void setModified2(const QModelIndex&, const QModelIndex&, const QVector<int>&) {
        setModified();
    }

signals:
    void modifiedChanged(bool);
    void optionsChanged();

public slots:
    void revert();
    void submit();
    void help();
    void selectIonZ();
    void validateOptions();

private:
    bool modified_{false};
    void applyRules();


    IonsUI* ionsui;
};

#endif // SIMULATIONOPTIONSVIEW_H
