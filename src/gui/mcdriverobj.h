#ifndef MCDRIVEROBJ_H
#define MCDRIVEROBJ_H

#include <chrono>

#include <QObject>
#include <QTimer>

#include "mcdriver.h"

class McDriverObj : public QObject
{
    Q_OBJECT

public:

    explicit McDriverObj();
    virtual ~McDriverObj() override;

    const mcdriver::options& options() const;
    void setOptions(const mcdriver::options& opt);
    std::string json() const;

    bool isModified() const { return modified_; }
    void setModified(bool b);
    bool isTemplate() const { return template_; }
    void setTemplate(bool b);

    QString fileName() const { return fileName_; }
    QString filePath() const { return filePath_; }
    void setFileName(const QString& path);
    void setTemplateName(const QString& name);

    QString title() const;

    bool validateOptions(QString* msg = nullptr) const;

    enum DriverStatus {
        mcReset = 0,
        mcIdle = 1,
        mcRunning = 2,
        mcMax = 3
    };

    DriverStatus status() const { return (DriverStatus)int(status_); };

    // load a json example or the default if path = QString()
    void loadJsonTemplate(const QString& path = QString());
    // load from a file on disk
    void loadJsonFile(const QString& path);
    bool loadH5File(const QString& path);
    QString ioErrorMsg() const;

    void saveJson(const QString& fname);
    bool saveH5(const QString& fname);

    void start(bool b);
    void reset();

    void init_run_data();
    void update_run_data();
    int progress() const { return progress_; }
    double elapsed() const { return total_elapsed_ + elapsed_; }
    double ips() const { return ips_; }
    double eta() const { return eta_; }
    size_t nions() const { return ncurr_; }
    const ArrayNDd& totals() const { return totals_; }
    const ArrayNDd& dtotals() const { return dtotals_; }

    const tally& getTally() const;
    const mccore* getSim() const;

    size_t maxIons() const { return max_ions_; }
    int nThreads() const { return nThreads_; }
    int seed() const { return seed_; }
    int updInterval() const { return updInterval_; }

public slots:
    void setMaxIons(int n) { max_ions_ = n; };
    void setNThreads(int n) { nThreads_ = n; };
    void setSeed(int n) { seed_ = n; }
    void setUpdInterval(int n) { updInterval_ = n; }

private slots:
    void start_();
    void onLoadH5_(const QString& path);
    void onSaveH5_();

signals:
    void modificationChanged(bool b);
    void configChanged();
    void contentsChanged();
    void fileNameChanged();
    void simulationCreated();
    void simulationDestroyed();
    void startSignal();
    void simulationStarted(bool b);
    void loadH5_(const QString& path);
    void saveH5_();

    // these are sent from within the worker/simulation thread
    // must be connected with Qt::QueuedConn.
    void statusChanged();
    void tallyUpdate();

private:
    mcdriver* driver_;
    mcdriver* test_driver_;
    mcdriver::options options_;
    bool modified_;
    bool template_;
    QString fileName_;
    QString filePath_;
    std::atomic_bool is_running_;
    std::atomic_int status_;
    std::atomic_bool io_op_active_;
    int io_ret_;
    std::string io_err_;

    size_t max_ions_;
    int nThreads_;
    int seed_;
    int updInterval_;

    // state
    typedef std::chrono::high_resolution_clock my_clock_t;
    typedef my_clock_t::time_point time_point;

    time_point tstart_; // time when sim started
    size_t nstart_; // # of ions when starting
    size_t ncurr_; // # of ions current
    size_t ntarget_; // max # of ions
    int progress_; // 0..100 in percent
    double elapsed_;
    double total_elapsed_;
    double eta_; // s
    double ips_; // ions per second
    ArrayNDd totals_, dtotals_;

    void setStatus(DriverStatus s);

    static void mc_callback_(const mcdriver& d, void* p);

    void update_tally_totals_();
};

#endif // MCDRIVEROBJ_H
