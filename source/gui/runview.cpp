#include "runview.h"
#include "simulationoptionsview.h"
#include "mainui.h"
#include "mydatawidgetmapper.h"
#include "optionsmodel.h"
#include "mcdriverobj.h"

#include "tally.h"

#include "error_fmt.h"

#include <QFile>
#include <QTimer>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QFormLayout>
#include <QGroupBox>
#include <QProgressBar>
#include <QToolButton>
#include <QLineEdit>
#include <QSpinBox>
#include <QLabel>
#include <QButtonGroup>
#include <QMessageBox>

RunView::RunView(MainUI *iui, QWidget *parent) : QWidget{ parent }, ionsui(iui)
{
    /* Create & Map widgets to OptionsModel */

    OptionsModel *model = ionsui->optionsModel;
    mapper = new MyDataWidgetMapper(model, this);

    // main title widget
    QLabel *simTitleLabel = new QLabel("Simulation title:");
    {
        QModelIndex idxOut = model->index("Output", 0);
        QModelIndex idxTitle = model->index("title", 0, idxOut);
        OptionsItem *item = model->getItem(idxTitle);
        simTitle = (QLineEdit *)item->createEditor(this);
        simTitle->setReadOnly(true);
        simTitleLabel->setToolTip(simTitle->toolTip());
        simTitleLabel->setWhatsThis(simTitle->whatsThis());
        mapper->addMapping(simTitle, idxTitle, item->editorSignal());
        simTitleLabel->setStyleSheet("font-size : 14pt; font-weight : bold;");
        simTitle->setStyleSheet("font-size : 14pt");
    }

    /* Create Info items */

    int ntotals = tally::std_tallies - 1;
    QStringList totalLabels;
    QColor clr = palette().color(QPalette::Window);
    QString styleSheet = QString("background: %1").arg(clr.name());
    for (int i = 0; i < ntotals; ++i) {
        QLineEdit *edt = new QLineEdit;
        edt->setReadOnly(true);
        edt->setStyleSheet(styleSheet);
        simTotals.push_back(edt);
        totalLabels << tally::arrayName(i + 1);
    }

    /* Layout items */

    box3 = new QGroupBox("Simulation Results");
    {
        int k = 0;
        int Nrow[] = { tally::cL, tally::eLost - tally::cL, tally::isCollision - tally::eLost };
        const char *titles[] = { "Point Defects", "Energy Loss", "Damage" };
        const char *units[] = { "(counts/ion)", "(eV/ion)", "(eV/ion)" };
        QHBoxLayout *hbox = new QHBoxLayout;
        for (int i = 0; i < 3; ++i) {
            QFormLayout *flayout = new QFormLayout;
            flayout->addRow(new QLabel(QString("%1 %2").arg(titles[i]).arg(units[i])));
            for (int j = 0; j < Nrow[i]; ++j) {
                flayout->addRow(totalLabels.at(k), simTotals[k]);
                k++;
            }
            hbox->addLayout(flayout);
            if (i != 2)
                hbox->addSpacing(10);
        }
        box3->setLayout(hbox);
        box3->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
    }

    QVBoxLayout *vbox = new QVBoxLayout;
    {
        QHBoxLayout *hbox = new QHBoxLayout;
        hbox->addWidget(simTitleLabel);
        hbox->addWidget(simTitle);
        hbox->addStretch();
        vbox->addLayout(hbox);
    }

    vbox->addWidget(box3);
    vbox->addStretch();
    setLayout(vbox);
    vbox->setContentsMargins(0, 0, 0, 0);

    /* Connect Signals */

    bool ret = connect(ionsui->driverObj(), &McDriverObj::statusChanged, this,
                       &RunView::onDriverStatusChanged, Qt::QueuedConnection);
    assert(ret);

    connect(ionsui->driverObj(), &McDriverObj::tallyUpdate, this, &RunView::onTallyUpdate,
            Qt::QueuedConnection);

    connect(ionsui->driverObj(), &McDriverObj::configChanged, this, &RunView::revert);
}

void RunView::revert()
{
    mapper->revert();
}

void RunView::onTallyUpdate()
{
    auto T = ionsui->driverObj()->totals();
    auto dT = ionsui->driverObj()->dtotals();
    if (!T.isNull()) {
        for (int i = 1; i < T.size(); ++i)
            simTotals[i - 1]->setText(QString::fromStdString(print_with_err(T[i], dT[i], 'g', 1)));
    }
}

void RunView::onDriverStatusChanged()
{
    const McDriverObj *D = ionsui->driverObj();
    McDriverObj::DriverStatus st = D->status();

    if (st == McDriverObj::mcReset)
        for (auto edt : simTotals)
            edt->setText("");
}
