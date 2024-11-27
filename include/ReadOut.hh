#include <G4Types.hh>
#include <array>
#include <cmath>
class ReadOut
{
    constexpr static size_t XSIZE = 32;
    constexpr static size_t YSIZE = 20;

    constexpr static float X_SIZE = 64*CLHEP::mm;
    constexpr static float Y_SIZE = 20*CLHEP::mm;

    std::array<std::array<G4float,XSIZE>, YSIZE> charge;

    public:

    ReadOut(void){
        std::fill(reinterpret_cast<float*>(charge.data()),reinterpret_cast<float*>(charge.data()) + XSIZE*YSIZE, 0.0);
    }


    G4float & operator() (float x, float y) 
    {
        int j=floor(x*(XSIZE/X_SIZE))+0.5*XSIZE;
        int i=floor(y*(YSIZE/Y_SIZE))+YSIZE*0.5;
        return charge[i][j];
    }
    void AddCharge(float x, float y, float q)
    {
        if (std::abs(x) > 0.5*X_SIZE || std::abs(y) > 0.5*Y_SIZE) return;
        (*this)(x,y) += q;
    }

};


