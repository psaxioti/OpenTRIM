#ifndef EVENT_RECORDER_H
#define EVENT_RECORDER_H

#include "simulation.h"

#include <vector>

class event_recorder
{
    const simulation* sim_;

    std::vector<float> buffer_;

    static const int event_size_ = 8;

public:

    event_recorder(const simulation* s);

    void record(simulation::simulation_event_t ev, const ion* i);

    void prepare(int nEvents);

    void flush();

};

#endif // EVENT_RECORDER_H
