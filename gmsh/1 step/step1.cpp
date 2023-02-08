#include <set>
#include <cmath>
#include <gmsh.h>
#include <cstdlib>
#include <vector>

int tor(double r1, double r2, double lc, int num, int i){
    gmsh::model::geo::addPoint(0, r1, 0, lc, i+1);
    gmsh::model::geo::addPoint(r1, 0, 0, lc, i+2);
    gmsh::model::geo::addPoint(0, -r1, 0, lc, i+3);
    gmsh::model::geo::addPoint(-r1, 0, 0, lc, i+4);//main centres

    gmsh::model::geo::addPoint(0, r1, r2, lc, i+5);
    gmsh::model::geo::addPoint(r1, 0, r2, lc, i+6);
    gmsh::model::geo::addPoint(0, -r1, r2, lc, i+7);
    gmsh::model::geo::addPoint(-r1, 0, r2, lc, i+8);//upper centres

    gmsh::model::geo::addPoint(0, r1, -r2, lc, i+9);
    gmsh::model::geo::addPoint(r1, 0, -r2, lc, i+10);
    gmsh::model::geo::addPoint(0, -r1, -r2, lc, i+11);
    gmsh::model::geo::addPoint(-r1, 0, -r2, lc, i+12);//down centres

    gmsh::model::geo::addPoint(0, r1 - r2, 0, lc, i+13);
    gmsh::model::geo::addPoint(r1 - r2, 0, 0, lc, i+14);
    gmsh::model::geo::addPoint(0, -r1 + r2, 0, lc, i+15);
    gmsh::model::geo::addPoint(-r1 + r2, 0, 0, lc, i+16);//inner centres

    gmsh::model::geo::addPoint(0, r1 + r2, 0, lc, i+17);
    gmsh::model::geo::addPoint(r1 + r2, 0, 0, lc, i+18);
    gmsh::model::geo::addPoint(0, -r1 - r2, 0, lc, i+19);
    gmsh::model::geo::addPoint(-r1 - r2, 0, 0, lc, i+20);//outer centres

    gmsh::model::geo::addPoint(0, 0, 0, lc, i+21);
    gmsh::model::geo::addPoint(0, 0, r2, lc, i+22);
    gmsh::model::geo::addPoint(0, 0, -r2, lc, i+23);//z-centres

    gmsh::model::geo::addCircleArc(i+11, i+23, i+12, i+29);
    gmsh::model::geo::addCircleArc(i+20, i+4, i+12, i+30);
    gmsh::model::geo::addCircleArc(i+16, i+4, i+12, i+31);
    gmsh::model::geo::addCircleArc(i+12, i+23, i+9, i+32);


    gmsh::model::geo::addCircleArc(i+7, i+22, i+8, i+25);
    gmsh::model::geo::addCircleArc(i+8, i+4, i+16, i+26);
    gmsh::model::geo::addCircleArc(i+20, i+4, i+8, i+27);
    gmsh::model::geo::addCircleArc(i+8, i+22, i+5, i+28);


    gmsh::model::geo::addCircleArc(i+14, i+2, i+10, i+21);
    gmsh::model::geo::addCircleArc(i+10, i+23, i+9, i+22);
    gmsh::model::geo::addCircleArc(i+10, i+2, i+18, i+23);
    gmsh::model::geo::addCircleArc(i+11, i+23, i+10, i+24);



    gmsh::model::geo::addCircleArc(i+14, i+2, i+6, i+17);
    gmsh::model::geo::addCircleArc(i+6, i+22, i+5, i+18);
    gmsh::model::geo::addCircleArc(i+6, i+2, i+18, i+19);
    gmsh::model::geo::addCircleArc(i+7, i+22, i+6, i+20);



    gmsh::model::geo::addCircleArc(i+13, i+21, i+16, i+13);
    gmsh::model::geo::addCircleArc(i+16, i+21, i+15, i+14);
    gmsh::model::geo::addCircleArc(i+19, i+21, i+20, i+15);
    gmsh::model::geo::addCircleArc(i+20, i+21, i+17, i+16);

    gmsh::model::geo::addCircleArc(i+15, i+3, i+11, i+9);
    gmsh::model::geo::addCircleArc(i+11, i+3, i+19, i+10);
    gmsh::model::geo::addCircleArc(i+17, i+1, i+9, i+11);
    gmsh::model::geo::addCircleArc(i+9, i+1, i+13, i+12);

    gmsh::model::geo::addCircleArc(i+13, i+21, i+14, i+1);
    gmsh::model::geo::addCircleArc(i+14, i+21, i+15, i+2);
    gmsh::model::geo::addCircleArc(i+15, i+3, i+7, i+3);
    gmsh::model::geo::addCircleArc(i+7, i+3, i+19, i+4);
    gmsh::model::geo::addCircleArc(i+19, i+21, i+18, i+5);
    gmsh::model::geo::addCircleArc(i+18, i+21, i+17, i+6);
    gmsh::model::geo::addCircleArc(i+17, i+1, i+5, i+7);
    gmsh::model::geo::addCircleArc(i+5, i+1, i+13, i+8);
    
    
    gmsh::model::geo::addCurveLoop({ i+10, i+15, i+30, -(i+29) }, i+2);
    gmsh::model::geo::addCurveLoop({ i+13, i+31, i+32, i+12 }, i+3);
    gmsh::model::geo::addCurveLoop({ i+16, i+11, -(i+32), -(i+30)}, i+4);
    gmsh::model::geo::addCurveLoop({ i+14, i+3, i+25, i+26 }, i+5);
    gmsh::model::geo::addCurveLoop({ i+13, -(i+26), i+28, i+8 }, i+6);
    gmsh::model::geo::addCurveLoop({ i+4, i+15, i+27, -(i+25) }, i+7);
    gmsh::model::geo::addCurveLoop({ i+16, i+7, -(i+28), -(i+27) }, i+8);
    gmsh::model::geo::addCurveLoop({ i+1, i+21, i+22, i+12 }, i+9);
    gmsh::model::geo::addCurveLoop({ i+6, i+11, -(i+22), i+23 }, i+10);
    gmsh::model::geo::addCurveLoop({ i+2, i+9, i+24, -(i+21) }, i+11);
    gmsh::model::geo::addCurveLoop({ i+10, i+5, -(i+23), -(i+24) }, i+12);
    gmsh::model::geo::addCurveLoop({ i+1, i+17, i+18, i+8 }, i+13);
    gmsh::model::geo::addCurveLoop({ i+6, i+7, -(i+18), i+19 }, i+14);
    gmsh::model::geo::addCurveLoop({ i+2, i+3, i+20, -(i+17) }, i+15);
    gmsh::model::geo::addCurveLoop({ i+4, i+5, -(i+19), -(i+20) }, i+16);
    gmsh::model::geo::addCurveLoop({ i+14, i+9, i+29, -(i+31) }, i+17);

    gmsh::model::geo::addSurfaceFilling({i+2}, i+2);
    gmsh::model::geo::addSurfaceFilling({i+3}, i+3);
    gmsh::model::geo::addSurfaceFilling({i+4}, i+4);
    gmsh::model::geo::addSurfaceFilling({i+5}, i+5);
    gmsh::model::geo::addSurfaceFilling({i+6}, i+6);
    gmsh::model::geo::addSurfaceFilling({i+7}, i+7);
    gmsh::model::geo::addSurfaceFilling({i+8}, i+8);
    gmsh::model::geo::addSurfaceFilling({i+9}, i+9);
    gmsh::model::geo::addSurfaceFilling({i+10}, i+10);
    gmsh::model::geo::addSurfaceFilling({i+11}, i+11);
    gmsh::model::geo::addSurfaceFilling({i+12}, i+12);
    gmsh::model::geo::addSurfaceFilling({i+13}, i+13);
    gmsh::model::geo::addSurfaceFilling({i+14}, i+14);
    gmsh::model::geo::addSurfaceFilling({i+15}, i+15);
    gmsh::model::geo::addSurfaceFilling({i+16}, i+16);
    gmsh::model::geo::addSurfaceFilling({i+17}, i+17);

    int ans = gmsh::model::geo::addSurfaceLoop({i+12, i+2, i+3, i+4, i+5, i+6, i+7, i+8, i+9, i+10, i+11, i+13, i+14, i+15, i+16, i+17 }, num);

    return ans;
}

int main(int argc, char** argv) {
    gmsh::initialize();

    gmsh::model::add("tor");

    double lc = 1e-1;
    int i = 0;
    tor(2, 0.5, lc, 1, i);
    i += 40;
    tor(2, 1, lc, 2, i);


    gmsh::model::geo::addVolume({ 2, -1 }, 1);

        
    gmsh::model::geo::synchronize();

    gmsh::model::mesh::generate(3);

    gmsh::write("tor.msh");

    std::set<std::string> args(argv, argv + argc);
    if (!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();

    return 0;

}
