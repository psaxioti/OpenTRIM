#ifndef ION_QUEUES_H
#define ION_QUEUES_H

#include <queue>

template<class ion_t>
class ion_queues {

public:

    // FIFO ion buffer
    typedef ion_t* ion_ptr_t;
    typedef std::queue< ion_t* > ion_queue_t;

private:
    ion_queue_t ion_buffer_; // buffer of allocated ion objects
    ion_queue_t recoil_queue_; // queue of generated recoils
    ion_queue_t pka_queue_; // queue of generated PKAs

    // pop an ion from the respective queue
    static ion_ptr_t pop_one_(ion_queue_t& Q)
    {
        if (Q.empty()) return nullptr;
        ion_ptr_t i = Q.front();
        Q.pop();
        return i;
    }

public:

    // return finished ion to buffer
    void free_ion(ion_ptr_t i) { ion_buffer_.push(i); }
    // delete all allocated ions
    void clear() {
        while (!ion_buffer_.empty()) {
            ion_ptr_t i = ion_buffer_.front();
            ion_buffer_.pop();
            delete i;
        }
    }
    ion_ptr_t new_ion(const ion_t& p) {
        ion_ptr_t i;
        if (ion_buffer_.empty()) i = new ion_t(p);
        else {
            i = ion_buffer_.front();
            ion_buffer_.pop();
            *i = p;
        }
        return i;
    }
    // ...
    ion_ptr_t pop_pka() { return pop_one_(pka_queue_); }
    // ...
    ion_ptr_t pop_recoil() { return pop_one_(recoil_queue_); }
    // push an ion. if recoil_id>0 goes to recoil queue
    void push_pka(ion_ptr_t i) { pka_queue_.push(i); }
    void push_recoil(ion_ptr_t i) { recoil_queue_.push(i); }

};

#endif // ION_QUEUES_H
