#ifndef PERIODIC_TABLE_H
#define PERIODIC_TABLE_H

#include <string>
#include <vector>
#include <map>

struct periodic_table
{
    enum class phase_t { gas, liquid, solid };

    struct isotope
    {
        int A;
        std::string symbol;
        double mass;
        double abundance;
        double half_life;
    };

    struct element
    {
        int Z;
        std::string symbol;
        std::string name;
        double mass;
        double density;
        phase_t phase;
        std::vector<isotope> isotopes;
        bool is_valid() const { return Z >= 0 && Z < periodic_table::size(); }
    };

    static int size() { return elements_.size(); }

    static const std::vector<element> &elements() { return elements_; }

    static const element &at(int Z)
    {
        return (Z >= 0 && Z < size()) ? elements_[Z] : invalid_element_;
    }

    static const element &at(const std::string &symbol)
    {
        auto i = symbol2z_.find(symbol);
        return (i != symbol2z_.end()) ? elements_[i->second] : invalid_element_;
    }

private:
    const static std::map<std::string, int> symbol2z_;
    const static std::vector<element> elements_;
    inline const static element invalid_element_{ -1 };
};

#endif // PERIODIC_TABLE_H
