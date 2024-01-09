#include "event_recorder.h"

event_recorder::event_recorder(const simulation* s) :
    sim_(s)
{

}


void event_recorder::record(simulation::simulation_event_t ev, const ion* i)
{

}

void event_recorder::prepare(int nEvents)
{
    buffer_.reserve(nEvents);
}

void event_recorder::flush()
{

}
