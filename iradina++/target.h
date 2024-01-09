#ifndef _TARGET_H_
#define _TARGET_H_

#include "geometry.h"
#include "arrays.h"
#include "urbg.h"

#include <string>
#include <vector>

class inventory;
class material;
class scatteringXS;
class reducedXS;

typedef enum {
    NoStraggling = 0,
    BohrStraggling,
    ChuStraggling,
    YangStraggling
} StragglingModel;

class atom {
    int id_;
    material* mat_;
    inventory* inv_;
    int Z_; // atomic number
    float M_; // mass
    float X_; // fraction, concentration
    float Ed_; // Displacement threshold energy (eV)
    float El_; // Lattice energy (eV)
    float Es_; // Surface binding energy (eV)
    float Er_; // Replacement threshold energy (eV)

public:

    int id() const { return id_; }
    const material* mat() const { return mat_; }

    int Z() const { return Z_; }
    float M() const { return M_; }
    float X() const { return X_; }
    float Ed() const { return Ed_; }
    float El() const { return El_; }
    float Es() const { return Es_; }
    float Er() const { return Er_; }

    bool operator==(const atom& other) const
    { return (Z_ == other.Z_) && (M_ == other.M_); }


private:
    atom();
    atom(class inventory* i, class material* m, int id);

    friend class inventory;
    friend class material;
};


class material {

    int id_;
    std::string name_;
    float massDensity_; // gr/cm^3

    inventory* inv_;
    std::vector<atom*> atoms_;
    std::vector<float> cumX_;

    float atomicDistance_; // nm
    float sqrtAtomicDistance_;
    float layerDistance_;

    float atomicDensity_; // at/cm^3 ??
    float atomicDensityNM_; // at/nm^3
    float sqrtRecFlDensity_;
    float meanZ_;
    float meanM_;

    // SRIM like stuff
    float meanF_;
    float meanA_;
    float meanMinRedTransfer_;
    float meanImpactPar_;

public:

    const std::string& name() const { return name_; }
    float atomicDensity() const { return atomicDensityNM_; }
    float atomicDistance() const { return atomicDistance_; }
    float layerDistance() const { return layerDistance_; }

    float meanF() const { return meanF_; }
    float meanA() const { return meanA_; }
    float meanMinRedTransfer() const { return meanMinRedTransfer_; }
    float meanImpactPar() const { return meanImpactPar_; }
    float sqrtRecFlDensity() const { return sqrtRecFlDensity_; }

    atom* addAtom(int Z, float M, float X,
                  float Ed, float El, float Es, float Er);


    void init(const reducedXS &xs);

    int id() const { return id_; }
    const std::vector<float>& X();
    const std::vector<atom*>& atoms() const { return atoms_; }

    const atom* selectAtom(URBG& g) const
    {
        if (atoms_.size()==1) return atoms_.front();
        float u = g.u01();
        int i=0;
        while((i < atoms_.size()-1) && (u > cumX_[i])) i++;
        return atoms_[i];
    }

private:
    material();
    material(class inventory* i, const char* name, const float& density, int id);

    friend class inventory;
};



class inventory {

    std::vector<atom*> atoms_;
    std::vector<material*> materials_;
    //std::vector<scatteringXS*> scattering_matrix_;
    Array2D<scatteringXS*> scattering_matrix_;
    Array3Df dedx_; // stopping data (atoms x materials x energy)
    Array2Df dedx1; // proton stopping (materials x energy)
    Array3Df de_strag_; // straggling data (atoms x materials x energy)

public:

    inventory();
    ~inventory();

    void setProjectile(int Z, float M);
    const atom* projectile() const { return atoms_.front(); }
    material* addMaterial(const char* name, const float &density);

    const std::vector<atom*>& atoms() { return atoms_; }
    const std::vector<material*>& materials() { return materials_; }

    void init(const reducedXS& xs, StragglingModel smodel = YangStraggling);

    const scatteringXS* getScatteringXS(const atom* z1, const atom* z2) const;
    int getDE(const atom* z1, const material* m, const float& erg,
                float &dedx, float &de_stragg) const;
    int getDEtables(const atom* z1, const material* m,
              const float *&dedx, const float *&de_stragg) const;

    friend class material;

};

class target
{
protected:

    // grid points
    grid3D grid_;

    // cells
    Array3D<const material*> cells_;

public:
    target();

    grid3D& grid() { return grid_; }
    const grid3D& grid() const { return grid_; }

    void fill(const box3D& box, const material* m);

    const material* cell(const ivector3& i) const
    { return cells_[i.x()][i.y()][i.z()]; }

};

#endif // TARGET_H
