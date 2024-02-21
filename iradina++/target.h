#ifndef _TARGET_H_
#define _TARGET_H_

#include "geometry.h"
#include "arrays.h"
#include "elements.h"

#include <string>
#include <vector>

class target;
class material;

class atom {

public:

    struct parameters {
        int Z; // atomic number
        float M; // mass
        float Ed; // Displacement threshold energy (eV)
        float El; // Lattice energy (eV)
        float Es; // Surface binding energy (eV)
        float Er; // Replacement threshold energy (eV)
    };

private:

    int id_;
    material* mat_;
    target* target_;
    parameters p_;
    float X_; // proxy concentration

    // LSS & NRT
    float lss_Kd_;
    float lss_Efact_;
    float nrt_L_;

public:

    int id() const { return id_; }
    const material* mat() const { return mat_; }

    const char* name() const { return elements::name(p_.Z); }
    int Z() const { return p_.Z; }
    float M() const { return p_.M; }
    float X() const { return X_; }
    float Ed() const { return p_.Ed; }
    float El() const { return p_.El; }
    float Es() const { return p_.Es; }
    float Er() const { return p_.Er; }

    bool operator==(const atom& other) const
    { return (p_.Z == other.p_.Z) && (p_.M == other.p_.M); }

    float LSS_Tdam(float recoilE) const;
    float NRT(float Tdam) const;

private:
    atom();
    atom(target* t, material* m, int id);

    friend class target;
    friend class material;
};


class material {

    int id_;
    std::string name_;
    float massDensity_; // gr/cm^3

    target* target_;
    std::vector<atom*> atoms_;
    std::vector<float> X_;
    std::vector<float> cumX_;

    float atomicDistance_; // nm
    float sqrtAtomicDistance_;
    float layerDistance_;

    // float atomicDensity_; // at/cm^3 ??
    float atomicDensityNM_; // at/nm^3
    float sqrtRecFlDensity_;
    float meanZ_;
    float meanM_;

    // SRIM like stuff
    float meanF_;
    float meanA_;
    float meanMinRedTransfer_;
    float meanImpactPar_;

    // LSS & NRT
    float lss_Kd_;
    float lss_Efact_;
    float nrt_Ed_;
    float nrt_L_;

public:

    struct material_desc_t {
        float density{1.f};
        bool isMassDensity{true};
        std::vector<int> Z;
        std::vector<float> M;
        std::vector<float> X;
        std::vector<float> Ed;
        std::vector<float> El;
        std::vector<float> Es;
        std::vector<float> Er;
    };

    /**
     * @brief Set atomic density of the material
     *
     * Invalidates other previous density setting
     *
     * @param v density in [at/nm^3]
     */
    void setAtomicDensity(float v) { atomicDensityNM_ = v; massDensity_ = -1; }
    /**
     * @brief Set mass density of the material
     *
     * Invalidates other previous density setting
     *
     * @param v density in [g/cm^3]
     */
    void setMassDensity(float v) { massDensity_ = v; atomicDensityNM_ = -1; }

    const std::string& name() const { return name_; }
    float atomicDensity() const { return atomicDensityNM_; }
    float massDensity() const { return massDensity_; }
    float atomicDistance() const { return atomicDistance_; }
    float layerDistance() const { return layerDistance_; }

    float meanF() const { return meanF_; }
    float meanA() const { return meanA_; }
    float meanMinRedTransfer() const { return meanMinRedTransfer_; }
    float meanImpactPar() const { return meanImpactPar_; }
    float sqrtRecFlDensity() const { return sqrtRecFlDensity_; }

    material_desc_t getDescription() const;

    atom* addAtom(const atom::parameters& p, float x);

    void init();

    int id() const { return id_; }
    const std::vector<float>& X();
    const std::vector<atom*>& atoms() const { return atoms_; }

    template<class _U>
    const atom* selectAtom(_U& g) const
    {
        if (atoms_.size()==1) return atoms_.front();
        float u = g.u01ropen();
        int i=0;
        while((i < atoms_.size()-1) && (u > cumX_[i])) i++;
        return atoms_[i];
    }

    float LSS_Tdam(float recoilE) const;
    float NRT(float Tdam) const;

private:
    material();
    material(class target* t, const char* name, int id);

    friend class target;
};

class target
{
protected:

    std::vector<atom*> atoms_;
    std::vector<material*> materials_;

    // grid points
    grid3D grid_;

    // cells
    Array3D<const material*> cells_;

    friend class material;

public:

    struct target_desc_t {
        std::vector< std::string > materials;
        std::vector< std::string > regions;
        ivector3 cell_count{1, 1, 1};
        ivector3 periodic_bc{0, 1, 1};
        vector3 cell_size{100.f, 100.f, 100.f};
    };

    struct region_desc_t {
        std::string material_id;
        vector3 min;
        vector3 max;
    };

protected:
    std::vector< region_desc_t > regions_;

public:
    target();
    target(const target& t);
    ~target();

    grid3D& grid() { return grid_; }
    const grid3D& grid() const { return grid_; }

    void fill(const box3D& box, const material* m);

    const std::vector< region_desc_t >& regions() const { return regions_; }
    target_desc_t getDescription() const;
    void getMaterialDescriptions(std::vector< material::material_desc_t >& mds);

    const material* cell(const ivector3& i) const
    { return cells_[i.x()][i.y()][i.z()]; }

    const material* cell(int i) const
    { return cells_.data()[i]; }

    void setProjectile(int Z, float M);
    const atom* projectile() const { return atoms_.front(); }
    material* addMaterial(const char* name);

    const std::vector<atom*>& atoms() const { return atoms_; }
    const std::vector<material*>& materials() const { return materials_; }

    void init();
};

#endif // TARGET_H
