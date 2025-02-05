#ifndef CASCADE_QUEUE_H
#define CASCADE_QUEUE_H

#include <queue>
#include <list>
#include <iostream>
#include <cassert>

#include "geometry.h"
#include "target.h"
#include "ion.h"
#include "tally.h"

class cascade_queue {

    enum defect_type_t {
        Vacancy,
        Interstitial
    };

    struct defect {
        defect_type_t type;
        double t;
        vector3 pos;
        vector3 dir;
        int icell;
        const atom* a;
        size_t pair_id;
    };

    struct interstitial : public defect {};
    struct vacancy : public defect {};

    // defect comparison operator
    // returns true of the time of lhs is greater than rhs
    struct defect_time_cmp {
        bool operator()(const defect* lhs, const defect* rhs)
        {
            return lhs->t > rhs->t;
        }
    };

    // type of a defect containers
    // std::list is used for optimized removal of items (regardless of position)
    typedef std::list<const interstitial*> interstitial_list_t;
    typedef std::list<const vacancy*> vacancy_list_t;
    // type for a container of recombination pairs
    typedef std::vector<std::pair<const interstitial*, const vacancy*>> pair_buff_t;
    // type for defect object buffers
    // std::vector selected for optimized allocation & pop_back/push_back at the back
    typedef std::vector<defect*> defect_buff_t;
    defect_buff_t all_buff_, free_buff_;

    // time ordered defect queue
    typedef std::priority_queue<defect*,
                                std::vector<defect*>,
                                defect_time_cmp> defect_queue_t;

    defect* new_defect(defect_type_t atype,
                       const double& at, const vector3& p, const vector3& n,
                       int ic, const atom* a, size_t pid) {
        defect* d;
        if (free_buff_.empty()) {
            d = new defect;
            all_buff_.push_back(d);
        } else {
            d = free_buff_.back();
            free_buff_.pop_back();
        }
        d->type = atype;
        d->t = at;
        d->pos = p;
        d->dir = n;
        d->icell = ic;
        d->a = a;
        d->pair_id = pid;
        return d;
    }



    defect_queue_t defect_queue;
    vacancy_list_t v_;
    interstitial_list_t i_;
    pair_buff_t riv_;
    size_t nv{0}, ni{0};
    ion pka;

    bool allow_correlated_recombination_;
    float i_rc_boost_;

    bool can_recombine(const defect* d1, const defect* d2)
    {
        return (d1->a == d2->a) && // same atom type
               ((!allow_correlated_recombination_ && (d1->pair_id != d2->pair_id)) || allow_correlated_recombination_);
    }

    bool recombine(const vacancy* d1, const grid3D& g)
    {
        float rc = d1->a->Rc(); // recombination radius
        auto i = i_.begin();
        for(; i!=i_.end(); ++i) {
            auto d2 = *i;
            if (can_recombine(d1,d2) &&
                g.distance(d1->pos,d2->pos) < rc)
            {
                riv_.emplace_back(d2,d1);
                i_.erase(i);
                return true;
            }
        }
        v_.push_back(d1);
        return false;
    }
    bool recombine(const interstitial* d1, const grid3D& g)
    {
        float rc = d1->a->Rc(); // recombination radius
        vector3 x = d1->pos;
        // special I recombination
        // x -= d1->dir*rc/2;
        rc *= i_rc_boost_;

        auto i = v_.begin();
        for(; i!=v_.end(); ++i) {
            auto d2 = *i;
            if (can_recombine(d1,d2) &&
                g.distance(x,d2->pos) < rc)
            {
                riv_.emplace_back(d1,d2);
                v_.erase(i);
                return true;
            }
        }
        i_.push_back(d1);
        return false;
    }

    static void print(std::ostream& os, const defect& d)
    {
        os << d.t << '\t'
           << d.a->id() << '\t'
           << ((d.type==Interstitial) ? '1' : '0') << '\t'
           << d.pos.x() << '\t'
           << d.pos.y() << '\t'
           << d.pos.z() << '\t'
           << d.icell << '\t'
           << d.pair_id << std::endl;
    }

public:
    cascade_queue(bool allow_correlated_recombination, float i_rc_boost) :
        allow_correlated_recombination_(allow_correlated_recombination),
        i_rc_boost_(i_rc_boost)
    {}
    ~cascade_queue() {
        for(defect* d : all_buff_) delete d;
    }

    void push_i(const ion* i)
    {

        defect_queue.push(new_defect(Interstitial, i->t(),
                                     i->pos(), i->dir(), i->cellid(),
                                     i->myAtom(), i->uid()));
        ni++;

    }
    void push_v(const ion* i)
    {

        defect_queue.push(new_defect(Vacancy, i->t(),
                                     i->pos(), i->dir(), i->cellid(),
                                     i->myAtom(), i->uid()));
        nv++;

    }

    void intra_cascade_recombination(const grid3D& g, tally& t)
    {
        //if (defect_queue.size()<=2) return;

#ifndef NDEBUG

        bool dbg = defect_queue.size() > 10;
        if (dbg) {
            std::cout << "PKA E=" << pka.erg() << std::endl;
            std::cout << "R=" << pka.pos() << std::endl;
            std::cout << "defect_queue (N=" << defect_queue.size() << "):" << std::endl;
        }

#endif

        for (; !defect_queue.empty(); defect_queue.pop()) {
            defect* d = defect_queue.top();

#ifndef NDEBUG
            if (dbg) print(std::cout, *d);
#endif

            if (d->type == Interstitial) recombine((interstitial*)d,g);
            else recombine((vacancy*)d,g);
        }

#ifndef NDEBUG
        if (dbg) {
            std::cout << "END defect_queue\n\n";
            std::cout << "Recombinations (" << riv_.size() << "):\n";
            for(auto& p : riv_) {
                print(std::cout, *p.first);
                print(std::cout, *p.second);

            }
        }
#endif

        // update tally due to I-V recombinations
        // get relevant tally arrays
        auto& AI = t.at(tally::cI);
        auto& AV = t.at(tally::cV);
        auto& As = t.at(tally::eStored);
        auto& Al = t.at(tally::eLattice);
        // update tally
        for(auto& p : riv_) {
            // interstitial
            auto& i = p.first;
            int iid = i->a->id();
            float El = i->a->El();
            int cid =  i->icell;
            AI(iid,cid)--;
            As(iid,cid) -= El/2;
            Al(iid,cid) += El/2;
            // vacancy
            auto& v = p.second;
            cid =  v->icell;
            AV(iid,cid)--;
            As(iid,cid) -= El/2;
            Al(iid,cid) += El/2;
        }
    }

    void clear() {

        // clear all buffers
        if (!defect_queue.empty()) {
            defect_queue_t empty;
            std::swap(defect_queue, empty);
        }
        i_.clear();
        v_.clear();
        riv_.clear();
        ni = nv = 0;
        free_buff_.clear();

        // if there are allocated defect objects copy them to the "free" buffer
        if (!all_buff_.empty()) {
            free_buff_.resize(all_buff_.size());
            std::copy(all_buff_.begin(),all_buff_.end(),free_buff_.begin());
        }

    }

    void init(const ion* i)
    {
        clear();
        pka = *i;
        defect_queue.push(new_defect(Vacancy, i->t(),
                                     i->pos0(), i->dir(),i->cellid0(),
                                     i->myAtom(), i->uid()));
        nv++;
    }

    void count_riv(float* s) const
    {
        for(auto& p : riv_) {
            auto& d = p.first; // interstitial
            int i = d->a->id()-1;
            s[i]++;
        }
    }

};





#endif // CASCADE_QUEUE_H
