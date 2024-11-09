#include "welcomeview.h"

#include "ionsui.h"

#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QButtonGroup>
#include <QFontMetrics>
#include <QListWidget>
#include <QListWidgetItem>
#include <QTextBrowser>
#include <QStackedWidget>
#include <QLabel>
#include <QDir>
#include <QDialogButtonBox>
#include <QFile>
#include <QFileDialog>
#include <QMessageBox>

WelcomeView::WelcomeView(IonsUI* iui, QWidget *parent)
    : QWidget{parent}, ionsui(iui)
{
    QFontMetrics fm = this->fontMetrics();
    QRect rect = fm.boundingRect("OOOOpen SimulationOOO");
    int h = rect.height()*2;
    int w = rect.width();
    btNew =     createButton("New Simulation",w,h,0);
    btOpen =    createButton("Open Simulation",w,h,0);
    btSave =    createButton("Save",w,h,0);
    btSaveAs =  createButton("Save As",w,h,0);
    btRecent =  createButton("Recent",w,h,1);
    btExamples = createButton("Examples",w,h,1);
    btGettingStarted = createButton("Getting Started",w,h,1);
    btAbout =   createButton("About",w,h,1);
    btRecent->setChecked(true);

    QVBoxLayout* vbox = new QVBoxLayout;
    vbox->addWidget(btNew);
    vbox->addWidget(btOpen);
    vbox->addWidget(btSave);
    vbox->addWidget(btSaveAs);
    vbox->addStretch();
    vbox->addWidget(btRecent);
    vbox->addWidget(btExamples);
    vbox->addWidget(btGettingStarted);
    vbox->addWidget(btAbout);
    buttonGrp = new QButtonGroup(this);
    buttonGrp->addButton(btRecent,0);
    buttonGrp->addButton(btExamples,1);
    buttonGrp->addButton(btGettingStarted,2);
    buttonGrp->addButton(btAbout,3);

    QHBoxLayout* viewLayout = new QHBoxLayout;
    viewLayout->addLayout(vbox,0);
    //viewLayout->addWidget(new QWidget,2);
    viewLayout->setContentsMargins(0,0,0,0);
    setLayout(viewLayout);

    /* create right stacked pane */
    stackedWidget = new QStackedWidget;
    viewLayout->addWidget(stackedWidget,2);

    pushCenterWidget("Recent simulations", new QWidget);

    /* create example view */
    QWidget* examplePage = new QWidget;
    {
        exampleList = new QListWidget;
        foreach( const QString &ex, QDir(":/examples").entryList() )
        {
            new QListWidgetItem(
                QIcon(":/icons/assets/ionicons/document-text-outline.svg"),
                ex,
                exampleList);
        }
        jsonView = new QTextBrowser;

        QDialogButtonBox* buttonBox = new QDialogButtonBox(Qt::Horizontal);
        QPushButton* btOpenExample = new QPushButton(QIcon(":/icons/assets/ionicons/checkmark-done-outline.svg"),
                                                  "Open example");
        buttonBox->addButton(btOpenExample, QDialogButtonBox::AcceptRole);
        connect(btOpenExample, &QPushButton::clicked,
                this, &WelcomeView::onOpenExample);

        QVBoxLayout* vbox = new QVBoxLayout;
        vbox->addWidget(exampleList);
        vbox->addSpacing(20);
        vbox->addWidget(jsonView);
        vbox->addWidget(buttonBox);
        vbox->setContentsMargins(0,0,0,0);
        examplePage->setLayout(vbox);
    }
    pushCenterWidget("Example simulations",examplePage);

    pushCenterWidget("Getting Started", new QWidget);
    pushCenterWidget("About", new QWidget);

    /* signal/slot connections */
    connect(buttonGrp, &QButtonGroup::idClicked,
            this, &WelcomeView::changeCenterWidget);
    connect(exampleList, &QListWidget::currentItemChanged,
            this, &WelcomeView::exampleSelected);

    connect(btOpen, &QPushButton::clicked,
            this, &WelcomeView::onOpenJson);
}

void WelcomeView::changeCenterWidget(int id)
{
    stackedWidget->setCurrentIndex(id);
}

void WelcomeView::onOpenExample()
{
    QListWidgetItem* i = exampleList->currentItem();
    if (!i) return;
    QFile f(QString(":/examples/%1").arg(i->text()));
    if (f.open(QIODevice::ReadOnly | QIODevice::Text)) {
        f.close();
        ionsui->ions_driver->loadJson(f.fileName());
        ionsui->setCurrentPage(IonsUI::idConfigPage);
    }
}

void WelcomeView::onOpenJson()
{
    QString fileName = QFileDialog::getOpenFileName(this,
                                            tr("Open Image"), QString(),
                                            tr("Json Files [*.json](*.json);;All Files [*.*](*.*)"));
    if (fileName.isNull()) return; // cancelled by user
    QFile f(fileName);
    if (f.open(QIODevice::ReadOnly | QIODevice::Text)) {
        f.close();
        ionsui->ions_driver->loadJson(fileName);
        ionsui->setCurrentPage(IonsUI::idConfigPage);
    }
    else QMessageBox::warning(this, "Open Failed",
                             QString("Could not open file:\n%1").arg(fileName),
                             QMessageBox::Ok);

}

QPushButton* WelcomeView::createButton(const QString& txt, int w, int h, int ch)
{
    QPushButton* bt = new QPushButton(txt);
    bt->setMinimumHeight(h);
    bt->setMinimumWidth(w);
    if (ch) {
        bt->setCheckable(true);
        bt->setAutoExclusive(true);
    }
    return bt;
}

void WelcomeView::pushCenterWidget(const QString &title, QWidget *page)
{
    QWidget* w = new QWidget;
    QVBoxLayout* vbox = new QVBoxLayout;
    QLabel* lbl = new QLabel(title);
    lbl->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::Fixed);
    lbl->setStyleSheet("font-size : 14pt; font-weight : bold;");
    vbox->addWidget(lbl);
    vbox->addWidget(page);
    w->setLayout(vbox);
    stackedWidget->addWidget(w);
}

void WelcomeView::exampleSelected()
{
    QListWidgetItem* i = exampleList->currentItem();
    QFile f(QString(":/examples/%1").arg(i->text()));
    if (f.open(QIODevice::ReadOnly | QIODevice::Text))
        jsonView->setText(f.readAll());
}
