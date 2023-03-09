#include <dolfin.h>
#include "heat.h"

using namespace dolfin;

class Source : public Expression
{
  public:
  double t = 0;
  double t0 = 3;
  double z0 = -45;

  void eval(Array<double>& values, const Array<double>& x) const
  {
    
    if (t >= t0 && x[0] < 25){
        values[0] = 300 * exp(-pow(((x[2] - z0) - 10 * (t- t0)), 2));
    }else{
        values[0] = 0;
    }
  }
};

class HandleBoundaryFunction: public Expression
{
    public:
    double t = 0;
    void eval(Array<double>& values, const Array<double>& x) const
    {
        values[0] = 310;
    }
};
class HandleBoundary: public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return on_boundary && (x[0] > 50 - DOLFIN_EPS);
    }
};
class HeadBoundaryFunction: public Expression
{
    public:
    double t = 0;
    void eval(Array<double>& values, const Array<double>& x) const
    {
        values[0] = 270;
    }
};
class HeadBoundary: public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return on_boundary && (x[0] < 25 - DOLFIN_EPS);
    }
};

int main(){

    unsigned mu = 5;
    double T = 20;
    unsigned num_steps = 60;
    double dt = T / num_steps;
    auto Dt = std::make_shared<Constant>(dt);
    auto Mu = std::make_shared<Constant>(mu);

    auto mesh = std::make_shared<Mesh>();
    auto mesh_file = std::make_shared<XDMFFile>(MPI_COMM_WORLD, "mesh/mjolnir.xdmf");
    mesh_file->read(*mesh);

    auto V = std::make_shared<heat::FunctionSpace>(mesh);

    auto handle_func = std::make_shared<HandleBoundaryFunction>();
    auto handle_boundary = std::make_shared<HandleBoundary>();

    auto head_func = std::make_shared<HeadBoundaryFunction>();
    auto head_boundary = std::make_shared<HeadBoundary>();

    DirichletBC handle_bc(V, handle_func, handle_boundary);
    DirichletBC head_bc(V, head_func, head_boundary);
    std::vector<const DirichletBC*> bcs = {{&handle_bc, &head_bc}};
    
    auto u_n = std::make_shared<Function>(V);
    u_n->interpolate(*head_func);

    auto f = std::make_shared<Source>();

    heat::BilinearForm a(V, V);
    heat::LinearForm L(V);
    a.Dt = Dt;
    a.mu = Mu;
    L.Dt = Dt;
    L.u_n = u_n;
    L.f = f;
    
    Function u(V);

    double t = 0;
    for (int i = 0; i < num_steps; ++i){
        t += dt;
        f->t = t;
        solve(a == L, u, bcs);

        *u_n = u;
        File file("heat-step-" + std::to_string(i) + ".pvd");
        file << u;
    }
    return 0;
}
